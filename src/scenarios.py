"""Scenario comparison: KatzBP vs OneHopBP under baseline, high load, a lab shock and blocking.

Each scenario runs the pilot's replications and seeds under PlainBP (FCFS at
the FCFS nodes) and shadow-evaluates one KatzBP and one OneHopBP. The key
metric is how often Katz and 1H choose different classes (Katz vs 1H).
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import metrics
from .model import build_model, calibrate
from .network import Network
from .policies import PlainBP
from .shadow import ShadowLogger, ThetaSpec, policy_key
from .sim import EDSimulation, Shock

COMPARISONS = {                     # metric stem -> epoch-log column (filled in per run)
    "r_katz_bp": "changed_{katz}",
    "r_1h_bp": "changed_oh_{hop}",
    "r_katz_1h": "katzvs1h_{hop}",
}


@dataclass
class Scenario:
    key: str
    label: str
    cfg: dict                 # full config with the scenario's model overrides merged in
    shock: dict | None


def deep_merge(base: dict, over: dict) -> dict:
    """Dicts merge key by key; any other value (including lists) replaces the base value."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        out[k] = deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def build_scenarios(cfg: dict) -> list[Scenario]:
    out = []
    for key, s in cfg["scenarios"]["list"].items():
        s = s or {}
        scfg = copy.deepcopy(cfg)
        scfg["model"] = deep_merge(cfg["model"], s.get("model"))
        out.append(Scenario(key, s.get("label", key), scfg, s.get("shock")))
    return out


def _shock(sc: Scenario, model, exp: dict) -> tuple[Shock | None, tuple[float, float] | None]:
    if not sc.shock:
        return None, None
    sh = sc.shock
    if sh["node"] not in model.nodes:
        raise ValueError(f"scenario {sc.key}: shock node {sh['node']} is not a model node")
    t0 = exp["warmup_minutes"] + float(sh["start_after_warmup_minutes"])
    window = (t0, t0 + float(sh["window_minutes"]))
    shock = Shock(model.nodes.index(sh["node"]), t0, t0 + float(sh["duration_minutes"]),
                  float(sh["rate_multiplier"]))
    if max(shock.end, window[1]) > exp["warmup_minutes"] + exp["run_minutes"]:
        raise ValueError(f"scenario {sc.key}: shock window ends after the run")
    return shock, window


def scenario_metrics(rep: int, df: pd.DataFrame, policy_nodes: list[str], katz: str, hop: str) -> list[dict]:
    """Contested share, tie share and the three disagreement rates, including and excluding ties."""
    out = []
    groups = [(metrics.ALL, df)] + [(n, df[df["node"] == n]) for n in policy_nodes]
    for node, d in groups:
        con = d[d["contested"]] if len(d) else d
        pos = con[con["bp_top2_gap"] > 0] if len(con) else con
        out.append(metrics._rec(rep, "epochs", len(d), node=node))
        out.append(metrics._rec(rep, "contested_epochs", len(con), node=node))
        out.append(metrics._rec(rep, "contested_share", len(con) / len(d) if len(d) else np.nan, node=node))
        out.append(metrics._rec(rep, "bp_tie_share", 1 - len(pos) / len(con) if len(con) else np.nan, node=node))
        for stem, col in COMPARISONS.items():
            col = col.format(katz=katz, hop=hop)
            out.append(metrics._rec(rep, stem, con[col].mean() if len(con) else np.nan, node=node))
            out.append(metrics._rec(rep, f"{stem}_nontie", pos[col].mean() if len(pos) else np.nan, node=node))
    return out


def run_scenarios(cfg: dict, out_dir: Path, say=print) -> dict:
    exp, scfg = cfg["experiment"], cfg["scenarios"]
    reps = int(exp["replications"])
    seeds = [int(s) for s in exp["seeds"]][:reps]
    kk, kt = float(scfg["katz"]["kappa"]), float(scfg["katz"]["theta"])
    hk, ht = float(scfg["one_hop"]["kappa"]), float(scfg["one_hop"]["theta_h"])
    katz_spec, hop_spec = ThetaSpec(f"abs{kt:g}", kt, "abs"), ThetaSpec(f"h{ht:g}", ht, "hop")
    katz_key, hop_key = policy_key(kk, katz_spec), policy_key(hk, hop_spec)

    rows, setup = {}, []          # rows: row key -> list of metric records
    for sc in build_scenarios(cfg):
        model = build_model(sc.cfg["model"])
        calibrate(model, float(sc.cfg["model"]["target_utilisation"]))
        network = Network.from_model(model)
        shock, window = _shock(sc, model, exp)
        policy_nodes = [model.nodes[i] for i in model.policy_nodes]
        (out_dir / "scenarios" / sc.key).mkdir(parents=True, exist_ok=True)
        whole, win, blocking, util, epochs = [], [], [], [], []
        for r, seed in enumerate(seeds, start=1):
            logger = ShadowLogger(r, seed, [kk], [katz_spec], hop_kappa=hk, hop_thetas=[hop_spec], katz_ref=katz_key)
            sim = EDSimulation(model, network, PlainBP(), seed, exp["warmup_minutes"], exp["run_minutes"],
                               observer=logger, shock=shock)
            res = sim.run()
            if res.arrivals != res.departures + res.in_system:
                raise AssertionError(f"scenario {sc.key} rep {r}: patient conservation violated")
            df = pd.DataFrame(logger.rows)
            df.to_csv(out_dir / "scenarios" / sc.key / f"epochs_rep{r}.csv.gz", index=False)
            whole += scenario_metrics(r, df, policy_nodes, katz_key, hop_key)
            if window:
                w = df[(df["time"] >= window[0]) & (df["time"] < window[1])]
                win += scenario_metrics(r, w, policy_nodes, katz_key, hop_key)
            blocking.append(res.blocking_events)
            util.append(res.utilisation)
            epochs.append(len(df))
        rows[sc.key] = whole
        if window:
            rows[f"{sc.key}-window"] = win
        bottleneck = max(model.utilisation, key=model.utilisation.get)
        setup.append({
            "scenario": sc.key, "label": sc.label, "lambda_per_min": float(model.arrival_rate),
            "arrivals_per_day": float(model.arrival_rate * 1440), "bottleneck": bottleneck,
            "bottleneck_rho": float(model.utilisation[bottleneck]), "rho_A": float(network.rho),
            "alpha": float(network.alpha(kk)), "blocking": model.blocking, "shock": dict(sc.shock) if sc.shock else None,
            "shock_window_minutes": list(window) if window else None,
            "analytic_utilisation": {n: round(u, 4) for n, u in model.utilisation.items()},
            "achieved_utilisation": {n: round(float(np.mean([u[n] for u in util])), 4) for n in model.nodes},
            "epochs_per_rep": float(np.mean(epochs)), "blocking_events_per_rep": float(np.mean(blocking)),
            "policy_nodes": policy_nodes,
        })
        say(f"scenario {sc.key} ({sc.label}): λ = {model.arrival_rate * 1440:.1f}/day, bottleneck {bottleneck} "
            f"ρ = {model.utilisation[bottleneck]:.2f}, ρ(A) = {network.rho:.4f}, "
            f"{np.mean(epochs):.0f} epochs/rep, {np.mean(blocking):.0f} blocking events/rep")

    summary = pd.concat([metrics.summarise(recs).assign(scenario=key) for key, recs in rows.items()],
                        ignore_index=True)
    summary = summary[["scenario"] + [c for c in summary.columns if c != "scenario"]]
    lo = cfg["verdict"]["r_change_min"]
    verdicts = {}
    for key in rows:
        r = summary[(summary["scenario"] == key) & (summary["metric"] == "r_katz_1h_nontie")
                    & (summary["node"] == metrics.ALL)]["mean"].iloc[0]
        verdicts[key] = ("n/a (no non-tied contested epochs)" if np.isnan(r)
                         else "KATZ DISTINCT" if r >= lo else "KATZ ≈ ONE-HOP")
    return {"summary": summary, "setup": setup, "verdicts": verdicts, "rows": list(rows),
            "katz": (kk, kt), "hop": (hk, ht), "labels": {s["scenario"]: s["label"] for s in setup}}


def lookup(summary: pd.DataFrame, row: str, metric: str, node: str = metrics.ALL) -> pd.Series | None:
    m = summary[(summary["scenario"] == row) & (summary["metric"] == metric) & (summary["node"] == node)]
    return m.iloc[0] if len(m) == 1 else None
