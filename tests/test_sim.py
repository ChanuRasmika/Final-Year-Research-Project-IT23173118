import copy

import pandas as pd
import pytest

from src.model import build_model, calibrate
from src.network import Network
from src.pilot import run_pilot
from src.policies import KatzBP, PlainBP
from src.sim import EDSimulation


@pytest.mark.parametrize("policy", [PlainBP(), KatzBP(5.0, 0.7)], ids=["BP", "KatzBP"])
def test_triage_never_violated(ed, policy):
    """At policy nodes, no epoch serves level k while a more urgent class waits."""
    model, network = ed
    violations, seen_nodes = [], set()

    def check(sim, state, choice):
        seen_nodes.add(state.node)
        most_urgent = min(c.level for c in state.candidates)
        if model.buffers[choice].level != most_urgent:
            violations.append((state.time, state.node, choice))

    res = EDSimulation(model, network, policy, seed=7, warmup=0, run_length=3 * 1440, observer=check).run()
    assert res.epochs > 500
    assert violations == []
    assert seen_nodes and not seen_nodes & model.fcfs_nodes   # FCFS nodes are never decision epochs


class _StartRecorder(EDSimulation):
    """Records (node, queue-join time) of every patient as service starts."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.starts = []

    def _start(self, node, b):
        self.starts.append((node, self.queues[b][0].t_join))
        super()._start(node, b)


def test_fcfs_nodes_serve_in_arrival_order(ed):
    model, network = ed
    assert {model.nodes[i] for i in model.fcfs_nodes} == {"registration", "triage"}
    sim = _StartRecorder(model, network, KatzBP(5.0, 0.7), seed=9, warmup=0, run_length=3 * 1440)
    sim.run()
    for node in model.fcfs_nodes:
        joins = [t for n, t in sim.starts if n == node]
        assert len(joins) > 500
        assert joins == sorted(joins)


def test_patient_conservation(ed):
    model, network = ed
    res = EDSimulation(model, network, PlainBP(), seed=8, warmup=1440, run_length=2 * 1440).run()
    assert res.arrivals > 0
    assert res.arrivals == res.departures + res.in_system


def test_mm1_tandem_matches_jackson():
    """Single class, 2 nodes, exponential service, rho = 0.7: L = rho / (1 - rho) per node."""
    cfg = {
        "nodes": ["A", "B"],
        "servers": {"A": 1, "B": 1},
        "service": {"A": {"mean": 1.0, "dist": "exponential"}, "B": {"mean": 1.0, "dist": "exponential"}},
        "pathways": {"P1": ["A", "B", "sink"]},
        "class_mix": {"level_probs": {3: 1.0}, "pathway_given_level": {3: {"P1": 1.0}}},
    }
    model = build_model(cfg)
    calibrate(model, 0.7)
    assert model.arrival_rate == pytest.approx(0.7)
    # Tandem is acyclic (rho(A) = 0), so Katz is undefined; BP does not need the network.
    res = EDSimulation(model, network=None, policy=PlainBP(), seed=123,
                       warmup=1_000, run_length=300_000).run()
    L = 0.7 / 0.3
    for node in ("A", "B"):
        assert res.mean_in_system[node] == pytest.approx(L, rel=0.05)
        assert res.utilisation[node] == pytest.approx(0.7, rel=0.02)


def test_same_seed_same_summary(cfg, tmp_path):
    small = copy.deepcopy(cfg)
    small["experiment"].update(replications=2, warmup_minutes=360, run_minutes=1440)
    small["katz"].update(kappa_grid=[0.5], theta_abs=[0.1, 1.0], theta_eta=[1])
    small["scenarios"]["list"]["C"]["shock"].update(start_after_warmup_minutes=600, window_minutes=480)
    ra = run_pilot(small, tmp_path / "a", make_plots=False, verbose=False)
    rb = run_pilot(small, tmp_path / "b", make_plots=False, verbose=False)
    a, b = ra["summary"], rb["summary"]
    timing = a["metric"].str.startswith("decision_time")
    pd.testing.assert_frame_equal(a[~timing].reset_index(drop=True), b[~timing].reset_index(drop=True))
    pd.testing.assert_frame_equal(ra["scenarios"]["summary"], rb["scenarios"]["summary"])
    assert (tmp_path / "a" / "pilot_summary.csv").exists()
