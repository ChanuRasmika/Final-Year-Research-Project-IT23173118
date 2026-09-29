import copy

import numpy as np
import pytest

from src.model import build_model, calibrate
from src.network import Network
from src.policies import KatzBP, PlainBP
from src.scenarios import build_scenarios, deep_merge
from src.sim import EDSimulation, Shock


def _scenario_model(cfg, key, **model_overrides):
    sc = next(s for s in build_scenarios(cfg) if s.key == key)
    mcfg = deep_merge(sc.cfg["model"], model_overrides)
    model = build_model(mcfg)
    calibrate(model, float(mcfg["target_utilisation"]))
    return model, Network.from_model(model)


class _CapacityChecker(EDSimulation):
    """_touch runs before every change to a node's counts and once at the end, so it sees every state."""

    def _touch(self, node):
        for i, cap in enumerate(self.capacity):
            assert self.waiting[i] <= cap, f"{self.model.nodes[i]}: {self.waiting[i]} waiting > buffer {cap}"
        super()._touch(node)


def test_scenario_d_adds_admission_to_network(cfg):
    model, network = _scenario_model(cfg, "D")
    adm = model.nodes.index("admission")
    assert model.blocking and model.capacity[adm] == 4
    assert network.A.shape == (8, 8) and network.A[model.nodes.index("treatment"), adm] > 0
    assert max(model.utilisation, key=model.utilisation.get) == "consultation"


@pytest.mark.parametrize("tight", [False, True], ids=["config-buffers", "buffers-of-1"])
def test_blocking_never_exceeds_buffer_and_conserves_patients(cfg, tight):
    over = {"capacity": {"lab": 1, "imaging": 1, "treatment": 1, "admission": 1}} if tight else {}
    model, network = _scenario_model(cfg, "D", **over)
    for policy in (PlainBP(), KatzBP(5.0, 0.5)):
        sim = _CapacityChecker(model, network, policy, seed=21, warmup=0, run_length=5 * 1440)
        res = sim.run()
        assert res.blocking_events > 0
        assert res.arrivals == res.departures + res.in_system
        assert sum(len(q) for q in sim.blocked) <= sum(sim.busy)   # blocked patients still hold a server


class _ServiceRecorder(EDSimulation):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.records = []

    def _base_service_time(self, b, node):
        base = super()._base_service_time(b, node)
        self._last_base = base
        return base

    def _service_time(self, b, node):
        final = super()._service_time(b, node)
        self.records.append((node, self.env.now, self._last_base, final))
        return final


def test_shock_changes_only_lab_service_within_window(ed):
    model, network = ed
    lab = model.nodes.index("lab")
    t0, t1 = 1440 + 2 * 1440, 1440 + 2 * 1440 + 240
    shock = Shock(lab, t0, t1, 0.5)
    shocked = _ServiceRecorder(model, network, PlainBP(), seed=31, warmup=1440, run_length=4 * 1440, shock=shock)
    shocked.run()
    in_window = 0
    for node, t, base, final in shocked.records:
        if node == lab and t0 <= t < t1:
            assert final == pytest.approx(2 * base)
            in_window += 1
        else:
            assert final == base
    assert in_window > 5
    # Same seed without the shock: identical service draws up to t0.
    plain = _ServiceRecorder(model, network, PlainBP(), seed=31, warmup=1440, run_length=4 * 1440)
    plain.run()
    before = [r for r in shocked.records if r[1] < t0]
    assert before == plain.records[:len(before)]
    assert shocked.shock_factor(lab, t0) == 2.0 and shocked.shock_factor(lab, t1) == 1.0
    assert all(shocked.shock_factor(i, t0 + 1) == 1.0 for i in range(model.n_nodes) if i != lab)


def test_deep_merge_replaces_lists_and_merges_dicts():
    base = {"a": {"x": 1, "y": 2}, "l": [1, 2]}
    assert deep_merge(base, {"a": {"y": 3}, "l": [9]}) == {"a": {"x": 1, "y": 3}, "l": [9]}
    assert base == {"a": {"x": 1, "y": 2}, "l": [1, 2]}
