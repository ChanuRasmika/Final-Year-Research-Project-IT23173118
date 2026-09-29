import pandas as pd

from src.metrics import summarise, theta_status, verdict
from src.shadow import ThetaSpec

VCFG = {"contested_min": 0.02, "r_change_min": 0.05, "r_change_max": 0.50, "dominance_max": 1.0}
THETAS = [ThetaSpec("abs0.01", 0.01, "abs"), ThetaSpec("abs0.25", 0.25, "abs"), ThetaSpec("abs5", 5.0, "abs")]


def _summary(r, rnt, dom):
    recs = [{"rep": 1, "metric": "contested_share", "kappa": "", "theta_label": "", "theta": float("nan"),
             "node": "ALL", "value": 0.3}]
    for t, a, b, c in zip(THETAS, r, rnt, dom):
        for metric, v in (("r_change", a), ("r_change_nontie", b), ("dominance_ratio", c)):
            recs.append({"rep": 1, "metric": metric, "kappa": 0.5, "theta_label": t.label,
                         "theta": t.value, "node": "ALL", "value": v})
    return summarise(recs)


def test_active_range_excludes_ties_and_tie_break_only_is_labelled():
    s = _summary(r=[0.22, 0.29, 0.60], rnt=[0.0, 0.12, 0.40], dom=[0.04, 0.8, 20.0])
    lines = verdict(s, [0.5], THETAS, VCFG)
    assert "ACTIVE RANGE (excluding ties): θ ∈ [0.25, 0.25], κ = 0.5" in lines
    assert any(l.startswith("TIE-BREAK ONLY κ = 0.5: θ ∈ {0.01}") for l in lines)
    assert any("centrality overrides BP" in l and "{5}" in l for l in lines)
    assert theta_status(s, 0.5, THETAS[0], VCFG) == "TIE-BREAK ONLY"
    assert theta_status(s, 0.5, THETAS[1], VCFG) == "ACTIVE"


def test_no_active_range_when_only_ties_change():
    s = _summary(r=[0.22, 0.25, 0.30], rnt=[0.0, 0.01, 0.02], dom=[0.04, 0.8, 0.9])
    lines = verdict(s, [0.5], THETAS, VCFG)
    assert not any(l.startswith("ACTIVE") for l in lines)
    assert any(l.startswith("NO ACTIVE RANGE") for l in lines)


def test_matched_theta_h_and_one_hop_verdict():
    from src.metrics import matched_theta_h, one_hop_verdict
    ref = ThetaSpec("abs0.1", 0.1, "abs")
    hops = [ThetaSpec("h0.05", 0.05, "hop"), ThetaSpec("h0.1", 0.1, "hop"), ThetaSpec("h0.25", 0.25, "hop")]
    recs = [{"rep": 1, "metric": "dominance_ratio", "kappa": 0.5, "theta_label": "abs0.1", "theta": 0.1,
             "node": "ALL", "value": 0.49}]
    for t, dom, dis in zip(hops, [0.3, 0.52, 1.3], [0.01, 0.03, 0.2]):
        recs += [{"rep": 1, "metric": m, "kappa": 0.5, "theta_label": t.label, "theta": t.value,
                  "node": "ALL", "value": v} for m, v in (("dominance_ratio_1h", dom), ("r_katz_vs_1h_nontie", dis))]
    s = summarise(recs)
    matched, katz_dom, hop_dom = matched_theta_h(s, 0.5, hops, 0.5, ref)
    assert matched.label == "h0.1" and katz_dom == 0.49 and hop_dom == 0.52
    assert one_hop_verdict(s, 0.5, matched, VCFG) == ["KATZ ≈ ONE-HOP: full-network signal adds few distinct decisions."]
    assert one_hop_verdict(s, 0.5, hops[2], VCFG)[0].startswith("KATZ ≠ ONE-HOP")
