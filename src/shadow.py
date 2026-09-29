"""Shadow evaluation: at each post-warm-up decision epoch of the PlainBP
trajectory, record the BP choice and the KatzBP choice for every (kappa, theta)
computed on the same state, without acting on it. Optionally also the
OneHopBP choice for each theta_h, and whether it disagrees with a reference
KatzBP policy.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from .policies import DecisionState, KatzBP, OneHopBP, PlainBP, bp_pressure, feasible_set, rank
from .sim import STREAM_SHADOW, EDSimulation, make_rng


@dataclass(frozen=True)
class ThetaSpec:
    label: str      # e.g. "abs0.1", "eta1" or "h0.1"
    value: float
    kind: str       # "abs", "eta" or "hop" (a OneHopBP theta_h)


def policy_key(kappa: float, theta: ThetaSpec) -> str:
    return f"k{kappa:g}_{theta.label}"


class ShadowLogger:
    def __init__(self, rep: int, seed: int, kappas: list[float], thetas: list[ThetaSpec],
                 hop_kappa: float | None = None, hop_thetas: list[ThetaSpec] = (),
                 katz_ref: str | None = None):
        """``katz_ref`` is the policy_key of the KatzBP policy that OneHopBP is compared with."""
        self.rep = rep
        self.kappas = list(kappas)
        self.thetas = list(thetas)
        self.bp = PlainBP()
        self.shadow = {policy_key(k, t): KatzBP(t.value, k) for k in self.kappas for t in self.thetas}
        self.hop_kappa = hop_kappa
        self.hop = {policy_key(hop_kappa, t): OneHopBP(t.value, hop_kappa) for t in hop_thetas}
        if self.hop and katz_ref not in self.shadow:
            raise ValueError(f"Katz reference {katz_ref!r} is not among the shadow KatzBP policies")
        self.katz_ref = katz_ref
        self.rng = make_rng(seed, STREAM_SHADOW)
        self.rows: list[dict] = []
        self.katz_ns = 0
        self.katz_calls = 0
        self.bp_ns = 0
        self.bp_calls = 0

    def __call__(self, sim: EDSimulation, state: DecisionState, choice: int) -> None:
        if sim.env.now < sim.warmup:
            return
        node = state.node
        model = sim.model

        t0 = time.perf_counter_ns()
        bp_choice = self.bp.select(node, state, self.rng)
        self.bp_ns += time.perf_counter_ns() - t0
        self.bp_calls += 1
        if bp_choice != choice:
            raise AssertionError("shadow PlainBP disagrees with the acting policy")

        F = feasible_set(state)
        n_distinct = len({c.next_node for c in F})
        row = {
            "rep": self.rep,
            "time": state.time,
            "node": model.nodes[node],
            "triage_level": F[0].level,
            "n_feasible": len(F),
            "n_distinct_next": n_distinct,
            "contested": len(F) >= 2 and n_distinct >= 2,
            "bp_choice": bp_choice,
            "bp_top2_gap": np.nan,
        }
        top2 = None
        if len(F) >= 2:
            P = [bp_pressure(c) for c in F]
            order = rank(F, P)
            row["bp_top2_gap"] = abs(P[order[0]] - P[order[1]])
            top2 = (F[order[0]].next_node, F[order[1]].next_node)
        for k in self.kappas:
            gap = np.nan
            if top2 is not None:
                c_hat = state.network.c_hat(state.q, k)
                gap = abs(c_hat[top2[0]] - c_hat[top2[1]])
            row[f"cent_gap_top2_k{k:g}"] = gap
        if self.hop:
            gap = np.nan
            if top2 is not None:
                h_hat = state.network.h_hat(state.q, self.hop_kappa)
                gap = abs(h_hat[top2[0]] - h_hat[top2[1]])
            row[f"hop_gap_top2_k{self.hop_kappa:g}"] = gap
        row["q_total"] = float(state.q.sum())
        for i, name in enumerate(model.nodes):
            row[f"q_{name}"] = int(state.q[i])

        for key, pol in self.shadow.items():
            t0 = time.perf_counter_ns()
            kc = pol.select(node, state, self.rng)
            self.katz_ns += time.perf_counter_ns() - t0
            self.katz_calls += 1
            row[f"katz_{key}"] = kc
            row[f"changed_{key}"] = kc != bp_choice
        for key, pol in self.hop.items():
            hc = pol.select(node, state, self.rng)
            row[f"oh_{key}"] = hc
            row[f"changed_oh_{key}"] = hc != bp_choice
            row[f"katzvs1h_{key}"] = hc != row[f"katz_{self.katz_ref}"]
        self.rows.append(row)

    @property
    def katz_us_per_call(self) -> float:
        return self.katz_ns / self.katz_calls / 1e3 if self.katz_calls else float("nan")

    @property
    def bp_us_per_call(self) -> float:
        return self.bp_ns / self.bp_calls / 1e3 if self.bp_calls else float("nan")
