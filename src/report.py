"""Markdown report and console summary for the pilot."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .metrics import lookup, theta_status
from .shadow import ThetaSpec


def fmt(row: pd.Series, pct: bool = False, digits: int = 3) -> str:
    m, h = row["mean"], row["ci_halfwidth"]
    if np.isnan(m):
        return "n/a"
    if pct:
        return f"{m:.1%} ± {h:.1%}" if not np.isnan(h) else f"{m:.1%}"
    return f"{m:.{digits}g} ± {h:.{digits}g}" if not np.isnan(h) else f"{m:.{digits}g}"


def theta_header(t: ThetaSpec) -> str:
    return f"{t.value:g}" if t.kind == "abs" else f"η={t.label[3:]} ({t.value:.3g})"


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def _kappa_theta_table(summary, metric, kappas, order, pct=True) -> str:
    return _table(["κ \\ θ"] + [theta_header(t) for t in order],
                  [[f"{k:g}"] + [fmt(lookup(summary, metric, k, t.label), pct=pct) for t in order]
                   for k in kappas])


def write_report(path, summary: pd.DataFrame, verdict_lines: list[str], ctx: dict) -> None:
    kappas, thetas, nodes = ctx["kappas"], ctx["thetas"], ctx["nodes"]
    pnodes, kc, cand = ctx["policy_nodes"], ctx["cand_kappa"], ctx["cand_specs"]
    k0 = ctx["kappa_default"]
    order = sorted(thetas, key=lambda t: t.value)
    L = []
    L.append("# Decision-impact pilot: Katz-BP vs plain BP (shadow evaluation)\n")
    L.append("## Verdict\n")
    L += [f"- **{v}**" for v in verdict_lines]
    L.append("\nThe ACTIVE RANGE uses R_change **excluding exact BP ties** (bp_top2_gap > 0) ≥ 5% together "
             "with a dominance ratio < 1. TIE-BREAK ONLY marks θ whose R_change reaches 5% only because of "
             "tied BP pressures, where the centrality term replaces the head-of-line tie-break.\n")
    L.append("> Model parameters are PROVISIONAL placeholders. This pilot asks only whether the "
             "centrality term changes decisions; it says nothing about congestion.\n")

    if cand:
        L.append(f"## Candidate set (κ = {kc:g})\n")
        L.append(_table(["θ", "R_change incl. ties", "R_change excl. ties", "dominance ratio", "status"],
                        [[f"{t.value:g}", fmt(lookup(summary, 'r_change', kc, t.label), pct=True),
                          fmt(lookup(summary, 'r_change_nontie', kc, t.label), pct=True),
                          fmt(lookup(summary, 'dominance_ratio', kc, t.label)),
                          theta_status(summary, kc, t, ctx["vcfg"])] for t in cand]))
        L.append("")

    L.append("## Set-up\n")
    L.append(f"- Replications: {ctx['reps']} (seeds {ctx['seeds']}); warm-up {ctx['warmup_h']:g} h discarded; "
             f"{ctx['run_days']:g} observed days per replication.")
    L.append(f"- Arrival rate λ = {ctx['lambda']:.4f} /min ({ctx['lambda'] * 1440:.1f} patients/day), "
             f"calibrated so the bottleneck ({ctx['bottleneck']}) has ρ = {ctx['target']:.2f}.")
    L.append(f"- Served FCFS (no policy, no decision epochs): {', '.join(ctx['fcfs_nodes']) or 'none'}. "
             f"Policy nodes: {', '.join(pnodes)}.")
    L.append(f"- Spectral radius ρ(A) = {ctx['rho_A']:.4f}; α = κ/ρ(A): "
             + ", ".join(f"κ={k:g} → α={ctx['alphas'][k]:.4f}" for k in kappas) + ".")
    L.append(f"- Median top-2 BP gap (rep 1, contested, post warm-up) = {ctx['median_gap']:.4g} /min"
             f"{' (' + ctx['median_gap_note'] + ')' if ctx['median_gap_note'] else ''}; "
             "pressure-scaled θ: " + ", ".join(f"η={t.label[3:]} → θ={t.value:.4g}"
                                                for t in thetas if t.kind == "eta") + ".")
    L.append(f"- Post-warm-up decision epochs: {fmt(lookup(summary, 'epochs'), digits=4)} per replication.\n")

    L.append("### Node utilisation (analytic target vs achieved)\n")
    L.append(_table(["node", "servers", "analytic", "achieved (mean ± 95% CI)"],
                    [[n, str(ctx['servers'][n]), f"{lookup(summary, 'utilisation_target', node=n)['mean']:.3f}",
                      fmt(lookup(summary, 'utilisation_achieved', node=n))] for n in nodes]))
    L.append("")

    L.append("## Contested share\n")
    L.append(f"Overall: **{fmt(lookup(summary, 'contested_share'), pct=True)}** of decision epochs.\n")
    L.append(_table(["node", "epochs / rep", "contested share"],
                    [[n, fmt(lookup(summary, 'epochs', node=n), digits=4),
                      fmt(lookup(summary, 'contested_share', node=n), pct=True)] for n in pnodes]))
    L.append("\n![contested share by node](contested_share_by_node.png)\n")

    L.append("## R_change(κ, θ) including ties — changed ÷ contested epochs (overall)\n")
    L.append(_kappa_theta_table(summary, "r_change", kappas, order))
    L.append("\n![R_change vs theta](r_change_vs_theta.png)\n")

    L.append("## R_change(κ, θ) excluding exact BP ties — verdict basis\n")
    L.append("Changed ÷ contested epochs whose top-2 BP gap is > 0, i.e. decisions where the centrality "
             "term overrides a real pressure difference rather than breaking a tie.\n")
    L.append(_kappa_theta_table(summary, "r_change_nontie", kappas, order))
    L.append("")

    L.append(f"### R_change including ties, by node (κ = {k0:g})\n")
    L.append(_table(["node"] + [theta_header(t) for t in order],
                    [[n] + [fmt(lookup(summary, 'r_change', k0, t.label, node=n), pct=True) for t in order]
                     for n in pnodes]))
    L.append("\nAll κ and nodes are in `pilot_summary.csv` (metric `r_change`). "
             "`n/a` means the node had no contested epochs.\n")

    L.append("## Pressure and centrality gaps (contested epochs)\n")
    rows = [["bp_top2_gap (1/min)", fmt(lookup(summary, 'bp_top2_gap_median')),
             fmt(lookup(summary, 'bp_top2_gap_iqr'))]]
    rows += [[f"cent_gap_top2, κ={k:g}", fmt(lookup(summary, 'cent_gap_top2_median', k)),
              fmt(lookup(summary, 'cent_gap_top2_iqr', k))] for k in kappas]
    L.append(_table(["quantity", "median", "IQR"], rows))
    L.append("")
    L.append("### Exact BP ties among contested epochs, by node\n")
    L.append("At a tie (gap = 0) any θ > 0 with a non-zero centrality gap decides the outcome instead of "
             "the head-of-line tie-break.\n")
    L.append(_table(["node", "contested epochs / rep", "exact-tie share"],
                    [[n, fmt(lookup(summary, 'contested_epochs', node=n), digits=4),
                      fmt(lookup(summary, 'bp_tie_share', node=n), pct=True)] for n in pnodes + ["ALL"]]))
    L.append("\n![BP gap histogram](bp_top2_gap_hist.png)\n")

    L.append("## Dominance ratio — median of θ·cent_gap_top2 ÷ bp_top2_gap\n")
    L.append("Below 1 means BP pressure still dominates the typical contested decision "
             "(ties with gap = 0 excluded).\n")
    L.append(_kappa_theta_table(summary, "dominance_ratio", kappas, order, pct=False))
    L.append("")

    if ctx.get("hop"):
        L += _one_hop_section(summary, ctx["hop"])
    if ctx.get("scen"):
        L += _scenario_section(ctx["scen"], ctx["vcfg"])

    L.append("## Decision computation time\n")
    L.append(f"- KatzBP (Katz solve + argmax), per call: **{fmt(lookup(summary, 'decision_time_katz_us'))} µs**")
    L.append(f"- PlainBP (argmax only), per call: {fmt(lookup(summary, 'decision_time_bp_us'))} µs")
    L.append("- Wall-clock timings vary between machines and runs; they are excluded from determinism checks.\n")

    bad = summary[summary["metric"] == "changed_uncontested"]["mean"].sum()
    L.append("## Sanity checks\n")
    L.append(f"- Changed decisions outside contested epochs (must be 0): {bad:g}")
    L.append(f"- Shadow PlainBP matched the acting policy at every epoch (asserted during the run).\n")
    L.append("Metric definitions are in the README. Per-epoch logs are in `epochs/`.")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def _one_hop_section(summary: pd.DataFrame, hop: dict) -> list[str]:
    hk, ref, kk, m = hop["kappa"], hop["katz_spec"], hop["katz_kappa"], hop["matched"]
    L = ["## Katz vs one-hop\n",
         f"1H-BP uses h = q + α Aᵀq, the first two terms of the Katz series, with the same α = κ/ρ(A) "
         f"(κ = {hk:g}): ĥ = h / (1 + max|h|) and W_c = μ_i^c (Q_i^c − Q_j^c) + θ_h (ĥ_i − ĥ_j). "
         f"It is shadow-evaluated on the same PlainBP trajectory. \"Katz vs 1H\" is the share of contested "
         f"epochs where KatzBP (κ = {kk:g}, θ = {ref.value:g}) and 1H-BP choose different classes. "
         "\"Excl. ties\" restricts to contested epochs with bp_top2_gap > 0.\n"]
    rows = []
    for t in hop["thetas"]:
        rows.append([f"{t.value:g}" + (" (matched)" if t.label == m.label else ""),
                     fmt(lookup(summary, "r_change_1h", hk, t.label), pct=True),
                     fmt(lookup(summary, "r_change_1h_nontie", hk, t.label), pct=True),
                     fmt(lookup(summary, "r_katz_vs_1h", hk, t.label), pct=True),
                     fmt(lookup(summary, "r_katz_vs_1h_nontie", hk, t.label), pct=True),
                     fmt(lookup(summary, "dominance_ratio_1h", hk, t.label))])
    L.append(_table(["θ_h", "1H vs BP incl. ties", "1H vs BP excl. ties", "Katz vs 1H incl. ties",
                     "Katz vs 1H excl. ties", "dominance ratio (1H)"], rows))
    L.append(f"\n**Matched θ_h = {m.value:g}**: its 1H dominance ratio ({hop['matched_dom']:.3g}) is the "
             f"closest in the grid to KatzBP θ = {ref.value:g} ({hop['katz_dom']:.3g}). At this value KatzBP "
             f"and 1H-BP disagree on {fmt(lookup(summary, 'r_katz_vs_1h', hk, m.label), pct=True)} of contested "
             f"epochs including ties and {fmt(lookup(summary, 'r_katz_vs_1h_nontie', hk, m.label), pct=True)} "
             "excluding ties. For comparison, KatzBP itself changes "
             f"{fmt(lookup(summary, 'r_change_nontie', kk, ref.label), pct=True)} of non-tied contested "
             "decisions relative to PlainBP.\n")
    return L


def _pct(row) -> str:
    return "n/a" if row is None else fmt(row, pct=True)


def _scenario_section(scen: dict, vcfg: dict) -> list[str]:
    from .scenarios import lookup as slook
    s = scen["summary"]
    (kk, kt), (hk, ht) = scen["katz"], scen["hop"]
    L = ["## Scenario comparison\n",
         f"Each scenario re-runs the {len(scen['setup'])} set-ups below with the same seeds and replications, "
         f"under PlainBP with FCFS at registration and triage, and shadow-evaluates KatzBP (κ = {kk:g}, θ = {kt:g}) "
         f"and OneHopBP (κ = {hk:g}, θ_h = {ht:g}, the matched value). \"Katz vs 1H\" is the share of contested "
         "epochs where the two choose different classes; \"excl.\" restricts to contested epochs with "
         f"bp_top2_gap > 0. Verdict: KATZ DISTINCT if Katz vs 1H excl. ties ≥ {vcfg['r_change_min']:.0%}, "
         "else KATZ ≈ ONE-HOP.\n"]
    L.append(_table(["scenario", "λ (per day)", "bottleneck ρ (analytic)", "ρ(A)", "α", "epochs / rep",
                     "blocking events / rep", "notes"],
                    [[f"{d['scenario']} {d['label']}", f"{d['arrivals_per_day']:.1f}",
                      f"{d['bottleneck']} {d['bottleneck_rho']:.2f}", f"{d['rho_A']:.4f}", f"{d['alpha']:.4f}",
                      f"{d['epochs_per_rep']:.0f}", f"{d['blocking_events_per_rep']:.0f}",
                      _scenario_note(d)] for d in scen["setup"]]))
    L.append("")
    rows = []
    for key in scen["rows"]:
        label = key if key.endswith("-window") else f"{key} {scen['labels'][key]}"
        rows.append([label, f"{slook(s, key, 'contested_epochs')['mean']:.0f}", _pct(slook(s, key, "contested_share"))]
                    + [_pct(slook(s, key, m)) for stem in ("r_katz_bp", "r_1h_bp", "r_katz_1h")
                       for m in (stem, f"{stem}_nontie")]
                    + [f"**{scen['verdicts'][key]}**"])
    L.append(_table(["scenario", "contested epochs / rep", "contested share", "Katz vs BP incl.", "Katz vs BP excl.", "1H vs BP incl.",
                     "1H vs BP excl.", "Katz vs 1H incl.", "Katz vs 1H excl.", "verdict"], rows))
    L.append("")
    L.append("The C-window row covers only the 8 h after t0, so it rests on few contested epochs per "
             "replication; read its rates together with that count.\n")
    if any(d["blocking"] for d in scen["setup"]):
        L.append("In D, A is recomputed with the admission node. ρ(A) is unchanged because admission is not on a "
                 "cycle, so it adds only a zero eigenvalue. Achieved utilisation in D counts blocked time as busy "
                 "(see `resolved_config.yaml`).\n")
    L.append("### Katz vs 1H by node (incl. ties / excl. ties, mean across replications)\n")
    nodes = []
    for d in scen["setup"]:
        nodes += [n for n in d["policy_nodes"] if n not in nodes]
    def cell(key, n):
        a, b = slook(s, key, "r_katz_1h", n), slook(s, key, "r_katz_1h_nontie", n)
        if a is None or np.isnan(a["mean"]):
            return "n/a"
        return f"{a['mean']:.1%} / " + ("n/a" if np.isnan(b["mean"]) else f"{b['mean']:.1%}")
    L.append(_table(["node"] + scen["rows"], [[n] + [cell(k, n) for k in scen["rows"]] for n in nodes]))
    L.append("\n`n/a` means no contested epochs at that node (or the node does not exist in that scenario). "
             "Confidence intervals for every cell are in `scenario_summary.csv`.\n")
    chk = scen.get("baseline_check")
    if chk:
        ok = "matches" if abs(chk[0] - chk[1]) < 1e-12 else "DOES NOT match"
        L.append(f"Consistency: scenario A's Katz vs 1H excl. ties ({chk[1]:.2%}) {ok} the baseline pilot's "
                 f"value at θ_h = {ht:g} ({chk[0]:.2%}); both shadow the same PlainBP trajectory.\n")
    return L


def _scenario_note(d: dict) -> str:
    notes = []
    if d["shock"]:
        sh = d["shock"]
        notes.append(f"{sh['node']} rate ×{sh['rate_multiplier']:g} for {sh['duration_minutes'] / 60:g} h from "
                     f"t0 = warm-up + {sh['start_after_warmup_minutes'] / 1440:g} d; window "
                     f"[t0, t0 + {sh['window_minutes'] / 60:g} h]")
    if d["blocking"]:
        notes.append("blocking-after-service with admission node")
    return "; ".join(notes) or "—"


def console_summary(summary: pd.DataFrame, verdict_lines: list[str], ctx: dict, runtime: float) -> str:
    kappas, thetas, pnodes = ctx["kappas"], ctx["thetas"], ctx["policy_nodes"]
    kc, cand = ctx["cand_kappa"], ctx["cand_specs"]
    order = sorted(thetas, key=lambda t: t.value)
    per_node = ", ".join(f"{n} {lookup(summary, 'contested_share', node=n)['mean']:.1%}" for n in pnodes
                         if lookup(summary, 'contested_share', node=n)['mean'] > 0)
    heads = [f"{t.value:.3g}" + ("*" if t.kind == "eta" else "") for t in order]
    out = [f"Contested share: {fmt(lookup(summary, 'contested_share'), pct=True)}  ({per_node})",
           "R_change %  θ→ " + " ".join(f"{h:>6}" for h in heads) + "   (incl. ties; * pressure-scaled)"]
    for k in kappas:
        vals = [lookup(summary, "r_change", k, t.label)["mean"] * 100 for t in order]
        out.append(f"  κ={k:<8g}  " + " ".join(f"{v:6.1f}" for v in vals))
    if cand:
        out.append(f"Candidate set κ={kc:g}   θ: " + " ".join(f"{t.value:>6g}" for t in cand))
        for metric, name in (("r_change", "incl. ties"), ("r_change_nontie", "excl. ties")):
            vals = [lookup(summary, metric, kc, t.label)["mean"] * 100 for t in cand]
            out.append(f"  R_change % {name}: " + " ".join(f"{v:6.1f}" for v in vals))
    hop = ctx.get("hop")
    if hop:
        m, hk = hop["matched"], hop["kappa"]
        out.append(f"Katz(θ={hop['katz_spec'].value:g}) vs 1H at matched θ_h={m.value:g}: "
                   f"{lookup(summary, 'r_katz_vs_1h', hk, m.label)['mean']:.1%} incl. ties, "
                   f"{lookup(summary, 'r_katz_vs_1h_nontie', hk, m.label)['mean']:.1%} excl. ties "
                   f"(dominance {hop['katz_dom']:.2f} vs {hop['matched_dom']:.2f})")
    out.append(f"Median BP gap {lookup(summary, 'bp_top2_gap_median')['mean']:.3g}/min, ties "
               f"{lookup(summary, 'bp_tie_share')['mean']:.1%}; Katz decision "
               f"{lookup(summary, 'decision_time_katz_us')['mean']:.1f} µs")
    scen = ctx.get("scen")
    if scen:
        from .scenarios import lookup as slook
        for key in scen["rows"]:
            a = slook(scen["summary"], key, "r_katz_1h")["mean"]
            b = slook(scen["summary"], key, "r_katz_1h_nontie")["mean"]
            out.append(f"Scenario {key:<9} Katz vs 1H {a:6.1%} incl. / {b:6.1%} excl. ties")
    out += [f"VERDICT: {v}" for v in verdict_lines]
    out.append(f"Total runtime: {runtime:.1f} s")
    return "\n".join(out)
