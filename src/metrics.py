"""Per-replication pilot metrics and their mean +- 95% CI across replications."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .shadow import ThetaSpec, policy_key

ALL = "ALL"
NA = ""   # placeholder key for metrics without a kappa / theta / node dimension


def _rec(rep, metric, value, kappa=NA, theta: ThetaSpec | None = None, node=ALL):
    return {"rep": rep, "metric": metric, "kappa": kappa,
            "theta_label": theta.label if theta else NA,
            "theta": theta.value if theta else np.nan,
            "node": node, "value": float(value)}


def _iqr(x: pd.Series) -> float:
    return float(x.quantile(0.75) - x.quantile(0.25)) if len(x) else np.nan


def _median(x: pd.Series) -> float:
    return float(x.median()) if len(x) else np.nan


def rep_metrics(rep: int, df: pd.DataFrame, nodes: list[str], policy_nodes: list[str],
                kappas: list[float], thetas: list[ThetaSpec], utilisation: dict[str, float],
                target: dict[str, float], katz_us: float, bp_us: float) -> list[dict]:
    """Metrics for one replication. ``df`` holds its post-warm-up epoch log.

    Per-node decision metrics cover ``policy_nodes`` only (FCFS nodes have no epochs).
    """
    out = []
    groups = [(ALL, df)] + [(n, df[df["node"] == n]) for n in policy_nodes]
    for node, d in groups:
        n_ep = len(d)
        n_con = int(d["contested"].sum())
        con_gap = d.loc[d["contested"], "bp_top2_gap"]
        out.append(_rec(rep, "epochs", n_ep, node=node))
        out.append(_rec(rep, "contested_epochs", n_con, node=node))
        out.append(_rec(rep, "contested_share", n_con / n_ep if n_ep else np.nan, node=node))
        out.append(_rec(rep, "bp_tie_share", (con_gap == 0).mean() if n_con else np.nan, node=node))
        for k in kappas:
            for t in thetas:
                changed = d[f"changed_{policy_key(k, t)}"]
                r = (changed & d["contested"]).sum() / n_con if n_con else np.nan
                out.append(_rec(rep, "r_change", r, kappa=k, theta=t, node=node))

    con = df[df["contested"]]
    gap = con["bp_top2_gap"]
    out.append(_rec(rep, "bp_top2_gap_median", _median(gap)))
    out.append(_rec(rep, "bp_top2_gap_iqr", _iqr(gap)))
    positive = con[con["bp_top2_gap"] > 0]
    for k in kappas:
        for t in thetas:
            # Changes that override a strictly positive BP gap (not just a tie-break); verdict basis.
            r = positive[f"changed_{policy_key(k, t)}"].mean() if len(positive) else np.nan
            out.append(_rec(rep, "r_change_nontie", r, kappa=k, theta=t))
        cg = con[f"cent_gap_top2_k{k:g}"]
        out.append(_rec(rep, "cent_gap_top2_median", _median(cg), kappa=k))
        out.append(_rec(rep, "cent_gap_top2_iqr", _iqr(cg), kappa=k))
        for t in thetas:
            ratio = t.value * positive[f"cent_gap_top2_k{k:g}"] / positive["bp_top2_gap"]
            out.append(_rec(rep, "dominance_ratio", _median(ratio), kappa=k, theta=t))
            # Sanity: a change outside a contested epoch is impossible by construction.
            bad = (df[f"changed_{policy_key(k, t)}"] & ~df["contested"]).sum()
            out.append(_rec(rep, "changed_uncontested", bad, kappa=k, theta=t))

    for n in nodes:
        out.append(_rec(rep, "utilisation_achieved", utilisation[n], node=n))
        out.append(_rec(rep, "utilisation_target", target[n], node=n))
    out.append(_rec(rep, "decision_time_katz_us", katz_us))
    out.append(_rec(rep, "decision_time_bp_us", bp_us))
    return out


def summarise(records: list[dict]) -> pd.DataFrame:
    """Mean and t-based 95% CI across replications for every metric key."""
    long = pd.DataFrame(records)
    keys = ["metric", "kappa", "theta_label", "theta", "node"]
    rows = []
    for key, g in long.groupby(keys, dropna=False, sort=False):
        v = g["value"].dropna().to_numpy()
        n = len(v)
        mean = v.mean() if n else np.nan
        hw = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
        rows.append(dict(zip(keys, key)) | {"mean": mean, "ci_halfwidth": hw,
                                            "ci_low": mean - hw, "ci_high": mean + hw, "n_reps": n})
    return pd.DataFrame(rows)


def lookup(summary: pd.DataFrame, metric: str, kappa=NA, theta_label=NA, node=ALL) -> pd.Series:
    m = summary[(summary["metric"] == metric) & (summary["kappa"] == kappa)
                & (summary["theta_label"] == theta_label) & (summary["node"] == node)]
    if len(m) != 1:
        raise KeyError((metric, kappa, theta_label, node))
    return m.iloc[0]


def verdict(summary: pd.DataFrame, kappas: list[float], thetas: list[ThetaSpec], vcfg: dict) -> list[str]:
    """Automatic verdict lines, evaluated on the across-replication means."""
    lines = []
    contested = lookup(summary, "contested_share")["mean"]
    if contested < vcfg["contested_min"]:
        lines.append("MODEL: insufficient branching/competition — revise ED topology or class mix.")

    lo = vcfg["r_change_min"]
    r = {(k, t.label): lookup(summary, "r_change", k, t.label)["mean"] for k in kappas for t in thetas}
    rnt = {(k, t.label): lookup(summary, "r_change_nontie", k, t.label)["mean"] for k in kappas for t in thetas}
    dom = {(k, t.label): lookup(summary, "dominance_ratio", k, t.label)["mean"] for k in kappas for t in thetas}
    if all(not (v >= lo) for v in r.values()):
        lines.append("BONUS TOO WEAK — rescale θ using the pressure gap.")

    # Active range: R_change EXCLUDING exact BP ties >= 5% and dominance ratio < 1.
    for k in kappas:
        active = sorted({t.value for t in thetas
                        if rnt[(k, t.label)] >= lo and dom[(k, t.label)] < vcfg["dominance_max"]})
        if active:
            lines.append(f"ACTIVE RANGE (excluding ties): θ ∈ [{active[0]:.4g}, {active[-1]:.4g}], κ = {k:g}")
    for k in kappas:
        tie_only = sorted({t.value for t in thetas if r[(k, t.label)] >= lo and not (rnt[(k, t.label)] >= lo)})
        if tie_only:
            vals = ", ".join(f"{v:.4g}" for v in tie_only)
            lines.append(f"TIE-BREAK ONLY κ = {k:g}: θ ∈ {{{vals}}} (R_change ≥ {lo:.0%} with ties, "
                         f"< {lo:.0%} without)")
    for k in kappas:
        strong = sorted({t.value for t in thetas if r[(k, t.label)] > vcfg["r_change_max"]})
        if strong:
            vals = ", ".join(f"{v:.4g}" for v in strong)
            lines.append(f"FLAG κ = {k:g}: θ ∈ {{{vals}}} have R_change > {vcfg['r_change_max']:.0%} "
                         "— centrality overrides BP — too strong.")
    if not any(l.startswith(("MODEL", "BONUS", "ACTIVE")) for l in lines):
        lines.append(f"NO ACTIVE RANGE: no θ has R_change excluding ties ≥ {vcfg['r_change_min']:.0%} "
                     f"with dominance ratio < {vcfg['dominance_max']:g}.")
    return lines


def theta_status(summary: pd.DataFrame, kappa: float, theta: ThetaSpec, vcfg: dict) -> str:
    """Per-(kappa, theta) label using the same rules as ``verdict``."""
    lo = vcfg["r_change_min"]
    r = lookup(summary, "r_change", kappa, theta.label)["mean"]
    rnt = lookup(summary, "r_change_nontie", kappa, theta.label)["mean"]
    dom = lookup(summary, "dominance_ratio", kappa, theta.label)["mean"]
    if rnt >= lo and dom < vcfg["dominance_max"]:
        label = "ACTIVE"
    elif rnt >= lo:
        label = "BP not dominant (ratio ≥ 1)"
    elif r >= lo:
        label = "TIE-BREAK ONLY"
    else:
        label = f"below {lo:.0%}"
    if r > vcfg["r_change_max"]:
        label += "; overrides BP (> 50%)"
    return label


# ---------- Katz vs one-hop ----------

def hop_metrics(rep: int, df: pd.DataFrame, hop_kappa: float, hop_thetas: list[ThetaSpec]) -> list[dict]:
    """1H-BP vs PlainBP and KatzBP(reference) vs 1H-BP, including and excluding exact BP ties."""
    out = []
    con = df[df["contested"]]
    pos = con[con["bp_top2_gap"] > 0]
    gap_col = f"hop_gap_top2_k{hop_kappa:g}"
    for t in hop_thetas:
        key = policy_key(hop_kappa, t)
        for metric, col in (("r_change_1h", f"changed_oh_{key}"), ("r_katz_vs_1h", f"katzvs1h_{key}")):
            out.append(_rec(rep, metric, con[col].mean() if len(con) else np.nan, kappa=hop_kappa, theta=t))
            out.append(_rec(rep, f"{metric}_nontie", pos[col].mean() if len(pos) else np.nan,
                            kappa=hop_kappa, theta=t))
        ratio = t.value * pos[gap_col] / pos["bp_top2_gap"]
        out.append(_rec(rep, "dominance_ratio_1h", _median(ratio), kappa=hop_kappa, theta=t))
        # Neither policy can differ from BP (or from each other) outside a contested epoch.
        bad = ((df[f"changed_oh_{key}"] | df[f"katzvs1h_{key}"]) & ~df["contested"]).sum()
        out.append(_rec(rep, "changed_uncontested", bad, kappa=hop_kappa, theta=t))
    return out


def matched_theta_h(summary: pd.DataFrame, hop_kappa: float, hop_thetas: list[ThetaSpec],
                    katz_kappa: float, katz_theta: ThetaSpec) -> tuple[ThetaSpec, float, float]:
    """theta_h whose mean 1H dominance ratio is closest (absolute difference) to the Katz reference's."""
    target = lookup(summary, "dominance_ratio", katz_kappa, katz_theta.label)["mean"]
    dom = {t.label: lookup(summary, "dominance_ratio_1h", hop_kappa, t.label)["mean"] for t in hop_thetas}
    best = min(hop_thetas, key=lambda t: (abs(dom[t.label] - target), t.value))
    return best, target, dom[best.label]


def one_hop_verdict(summary: pd.DataFrame, hop_kappa: float, matched: ThetaSpec, vcfg: dict) -> list[str]:
    r = lookup(summary, "r_katz_vs_1h_nontie", hop_kappa, matched.label)["mean"]
    if r < vcfg["r_change_min"]:
        return ["KATZ ≈ ONE-HOP: full-network signal adds few distinct decisions."]
    return [f"KATZ ≠ ONE-HOP: Katz and 1H-BP disagree on {r:.1%} of non-tied contested decisions "
            f"at matched θ_h = {matched.value:g}."]
