"""Decision-impact pilot. Usage: python -m src.pilot --config configs/pilot.yaml"""
from __future__ import annotations

import argparse
import copy
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import simpy
import yaml

from . import metrics, plots, report, scenarios
from .model import EDModel, build_model, calibrate
from .network import Network
from .policies import PlainBP
from .shadow import ShadowLogger, ThetaSpec, policy_key
from .sim import EDSimulation


def load_config(path: str | Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def setup(cfg: dict) -> tuple[EDModel, Network]:
    model = build_model(cfg["model"])
    calibrate(model, float(cfg["model"]["target_utilisation"]))
    return model, Network.from_model(model)


def run_rep(model: EDModel, network: Network, rep: int, seed: int, exp: dict,
            kappas: list[float], thetas: list[ThetaSpec], hop: dict | None = None):
    hop = hop or {}
    logger = ShadowLogger(rep, seed, kappas, thetas, hop_kappa=hop.get("kappa"),
                          hop_thetas=hop.get("thetas", []), katz_ref=hop.get("katz_ref"))
    sim = EDSimulation(model, network, PlainBP(), seed,
                       warmup=exp["warmup_minutes"], run_length=exp["run_minutes"], observer=logger)
    result = sim.run()
    return pd.DataFrame(logger.rows), result, logger


def scaled_thetas(model, network, cfg) -> tuple[float, str, list[ThetaSpec]]:
    """theta = eta x median top-2 BP gap over contested, post-warm-up epochs of replication 1."""
    exp = cfg["experiment"]
    df, _, _ = run_rep(model, network, 1, exp["seeds"][0], exp, kappas=[], thetas=[])
    gaps = df.loc[df["contested"], "bp_top2_gap"] if len(df) else pd.Series(dtype=float)
    if gaps.empty:
        raise RuntimeError("no contested epochs in replication 1; cannot pressure-scale theta")
    median, note = float(gaps.median()), ""
    if median == 0:
        median, note = float(gaps[gaps > 0].median()), "median was 0 (ties); used median of non-zero gaps"
    specs = [ThetaSpec(f"eta{eta:g}", float(eta) * median, "eta") for eta in cfg["katz"]["theta_eta"]]
    return median, note, specs


def candidate_set(kcfg: dict) -> tuple[float, list[float]]:
    """The named candidate (kappa, theta) set; its kappa must be in the grid."""
    cand = kcfg.get("candidate_set") or {}
    if not cand:
        return float(kcfg["kappa_default"]), []
    kappa = float(cand["kappa"])
    if kappa not in [float(k) for k in kcfg["kappa_grid"]]:
        raise ValueError(f"candidate_set kappa {kappa} is not in kappa_grid")
    return kappa, [float(t) for t in cand["theta"]]


def one_hop_setup(cfg: dict, kappas: list[float], thetas: list[ThetaSpec]) -> dict:
    """OneHopBP grid and the KatzBP reference it is compared with (empty if not configured)."""
    ocfg = cfg.get("one_hop") or {}
    if not ocfg:
        return {}
    ref = ocfg["katz_reference"]
    ref_kappa, ref_theta = float(ref["kappa"]), float(ref["theta"])
    ref_spec = next((t for t in thetas if t.kind == "abs" and t.value == ref_theta), None)
    if ref_kappa not in kappas or ref_spec is None:
        raise ValueError(f"one_hop.katz_reference (κ={ref_kappa:g}, θ={ref_theta:g}) is not in the Katz grid")
    return {"kappa": float(ocfg["kappa"]),
            "thetas": [ThetaSpec(f"h{t:g}", float(t), "hop") for t in ocfg["theta_h"]],
            "katz_kappa": ref_kappa, "katz_spec": ref_spec,
            "katz_ref": policy_key(ref_kappa, ref_spec)}


def run_pilot(cfg: dict, out_dir: str | Path, make_plots: bool = True, verbose: bool = True) -> dict:
    t_start = time.perf_counter()
    out = Path(out_dir)
    (out / "epochs").mkdir(parents=True, exist_ok=True)
    say = print if verbose else (lambda *a, **k: None)

    exp, kcfg = cfg["experiment"], cfg["katz"]
    reps = int(exp["replications"])
    seeds = [int(s) for s in exp["seeds"]]
    if len(seeds) < reps:
        raise ValueError("need one seed per replication")
    seeds = seeds[:reps]

    model, network = setup(cfg)
    kappas = [float(k) for k in kcfg["kappa_grid"]]
    target = model.utilisation
    bottleneck = max(target, key=target.get)
    say(f"λ = {model.arrival_rate:.4f}/min ({model.arrival_rate * 1440:.1f}/day); analytic utilisation:")
    for n, u in target.items():
        say(f"  {n:<13} {u:.3f}" + ("  <- bottleneck" if n == bottleneck else ""))
    policy_nodes = [model.nodes[i] for i in model.policy_nodes]
    say("FCFS nodes (no policy, no decision epochs): "
        + ", ".join(model.nodes[i] for i in sorted(model.fcfs_nodes)) if model.fcfs_nodes else "FCFS nodes: none")
    say(f"ρ(A) = {network.rho:.4f}; α = " + ", ".join(f"{network.alpha(k):.4f} (κ={k:g})" for k in kappas))

    median_gap, gap_note, eta_specs = scaled_thetas(model, network, cfg)
    cand_kappa, cand_thetas = candidate_set(kcfg)
    theta_abs = [float(t) for t in kcfg["theta_abs"]]
    theta_abs += [t for t in cand_thetas if t not in theta_abs]   # candidate thetas join the absolute grid
    thetas = [ThetaSpec(f"abs{t:g}", t, "abs") for t in theta_abs] + eta_specs
    cand_specs = [ThetaSpec(f"abs{t:g}", t, "abs") for t in cand_thetas]
    hop = one_hop_setup(cfg, kappas, thetas)
    say(f"Median top-2 BP gap (rep 1) = {median_gap:.4g}; scaled θ = "
        + ", ".join(f"{t.value:.4g}" for t in eta_specs))

    records, all_epochs = [], []
    for r, seed in enumerate(seeds, start=1):
        df, res, logger = run_rep(model, network, r, seed, exp, kappas, thetas, hop)
        if res.arrivals != res.departures + res.in_system:
            raise AssertionError(f"rep {r}: patient conservation violated")
        df.to_csv(out / "epochs" / f"epochs_rep{r}.csv.gz", index=False)
        records += metrics.rep_metrics(r, df, model.nodes, policy_nodes, kappas, thetas, res.utilisation, target,
                                       logger.katz_us_per_call, logger.bp_us_per_call)
        if hop:
            records += metrics.hop_metrics(r, df, hop["kappa"], hop["thetas"])
        all_epochs.append(df)
        say(f"rep {r} (seed {seed}): {res.arrivals} arrivals, {len(df)} epochs, "
            f"{int(df['contested'].sum()) if len(df) else 0} contested")

    summary = metrics.summarise(records)
    summary.to_csv(out / "pilot_summary.csv", index=False)
    verdict_lines = metrics.verdict(summary, kappas, thetas, cfg["verdict"])
    if hop:
        hop["matched"], hop["katz_dom"], hop["matched_dom"] = metrics.matched_theta_h(
            summary, hop["kappa"], hop["thetas"], hop["katz_kappa"], hop["katz_spec"])
        verdict_lines += metrics.one_hop_verdict(summary, hop["kappa"], hop["matched"], cfg["verdict"])

    scen = None
    if cfg.get("scenarios"):
        say("\nScenario comparison:")
        scen = scenarios.run_scenarios(cfg, out, say)
        scen["summary"].to_csv(out / "scenario_summary.csv", index=False)
        scen["baseline_check"] = _baseline_check(summary, scen)
        verdict_lines += [f"SCENARIO {k} ({scen['labels'].get(k.split('-')[0], k)}"
                          f"{', window' if k.endswith('-window') else ''}): {v}" for k, v in scen["verdicts"].items()]

    pd.DataFrame([{"class_id": b.id, "name": model.buffer_name(b.id), "level": b.level,
                   "pathway": b.pathway, "stage": b.stage, "node": model.nodes[b.node],
                   "next_node": model.node_name(b.next_node), "mean_service_min": b.mean,
                   "class_prob": model.classes[b.class_id].prob} for b in model.buffers]
                 ).to_csv(out / "classes.csv", index=False)

    resolved = copy.deepcopy(cfg)
    resolved["resolved"] = {
        "arrival_rate_per_min": float(model.arrival_rate),
        "arrivals_per_day": float(model.arrival_rate * 1440),
        "analytic_utilisation": {n: round(u, 6) for n, u in target.items()},
        "bottleneck": bottleneck,
        "spectral_radius_A": network.rho,
        "alpha": {float(k): network.alpha(k) for k in kappas},
        "adjacency_A": {"nodes": model.nodes, "matrix": network.A.round(6).tolist()},
        "median_top2_bp_gap": median_gap,
        "median_gap_note": gap_note,
        "theta_scaled": {float(t.label[3:]): t.value for t in eta_specs},
        "theta_abs_evaluated": theta_abs,
        "candidate_set": {"kappa": cand_kappa, "theta": cand_thetas},
        "fcfs_nodes": [model.nodes[i] for i in sorted(model.fcfs_nodes)],
        "policy_nodes": policy_nodes,
        "one_hop": ({"kappa": hop["kappa"], "theta_h": [t.value for t in hop["thetas"]],
                     "katz_reference": {"kappa": hop["katz_kappa"], "theta": hop["katz_spec"].value},
                     "matched_theta_h": hop["matched"].value,
                     "katz_dominance_ratio": float(hop["katz_dom"]),
                     "matched_dominance_ratio_1h": float(hop["matched_dom"])} if hop else None),
        "scenarios": scen["setup"] if scen else None,
        "seeds_used": seeds,
        "class_mix": {model.classes[c].name: round(model.classes[c].prob, 6) for c in range(len(model.classes))},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                     "pandas": pd.__version__, "simpy": simpy.__version__},
    }
    with open(out / "resolved_config.yaml", "w") as f:
        yaml.safe_dump(resolved, f, sort_keys=False, allow_unicode=True)

    ctx = {"kappas": kappas, "thetas": thetas, "nodes": model.nodes, "policy_nodes": policy_nodes,
           "fcfs_nodes": [model.nodes[i] for i in sorted(model.fcfs_nodes)],
           "cand_kappa": cand_kappa, "cand_specs": cand_specs, "vcfg": cfg["verdict"], "hop": hop, "scen": scen,
           "kappa_default": float(kcfg["kappa_default"]),
           "reps": reps, "seeds": seeds, "warmup_h": exp["warmup_minutes"] / 60, "run_days": exp["run_minutes"] / 1440,
           "lambda": model.arrival_rate, "bottleneck": bottleneck, "target": model.target_utilisation,
           "rho_A": network.rho, "alphas": {k: network.alpha(k) for k in kappas},
           "median_gap": median_gap, "median_gap_note": gap_note,
           "servers": dict(zip(model.nodes, model.servers.tolist()))}
    report.write_report(out / "pilot_report.md", summary, verdict_lines, ctx)
    if make_plots:
        plots.r_change_vs_theta(summary, kappas, thetas, cfg["verdict"], out / "r_change_vs_theta.png")
        plots.bp_gap_hist(pd.concat(all_epochs, ignore_index=True), out / "bp_top2_gap_hist.png")
        plots.contested_by_node(summary, policy_nodes, out / "contested_share_by_node.png")

    runtime = time.perf_counter() - t_start
    say("\n" + report.console_summary(summary, verdict_lines, ctx, runtime))
    return {"summary": summary, "verdict": verdict_lines, "thetas": thetas, "runtime": runtime,
            "scenarios": scen}


def _baseline_check(summary: pd.DataFrame, scen: dict) -> tuple[float, float] | None:
    """Scenario A shadows the same trajectory as the main pilot, so Katz vs 1H must match exactly."""
    (kk, kt), (hk, ht) = scen["katz"], scen["hop"]
    row = scenarios.lookup(scen["summary"], "A", "r_katz_1h_nontie")
    main = summary[(summary["metric"] == "r_katz_vs_1h_nontie") & (summary["kappa"] == hk)
                   & (summary["theta_label"] == f"h{ht:g}") & (summary["node"] == metrics.ALL)]
    if row is None or len(main) != 1:
        return None
    return float(main["mean"].iloc[0]), float(row["mean"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="Katz-BP decision-impact pilot (shadow evaluation)")
    ap.add_argument("--config", required=True, help="path to the pilot YAML config")
    ap.add_argument("--out", default=None, help="output directory (default: experiment.output_dir)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    run_pilot(cfg, args.out or cfg["experiment"].get("output_dir", "results"))


if __name__ == "__main__":
    main()
