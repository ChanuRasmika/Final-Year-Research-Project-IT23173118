import numpy as np

from src.policies import FCFS, KatzBP, PlainBP, feasible_set

from .conftest import random_state


def test_theta_zero_matches_plain_bp(ed):
    model, network = ed
    rng = np.random.default_rng(0)
    bp = PlainBP()
    katz = [KatzBP(0.0, k) for k in (0.3, 0.5, 0.7)]
    for _ in range(1000):
        s = random_state(model, network, rng)
        expected = bp.select(s.node, s, rng)
        for pol in katz:
            assert pol.select(s.node, s, rng) == expected


def test_choice_is_most_urgent_level(ed):
    model, network = ed
    rng = np.random.default_rng(1)
    policies = [PlainBP(), KatzBP(5.0, 0.7), KatzBP(100.0, 0.3)]
    for _ in range(1000):
        s = random_state(model, network, rng)
        r_star = min(c.level for c in s.candidates)
        for pol in policies:
            chosen = model.buffers[pol.select(s.node, s, rng)]
            assert chosen.level == r_star and chosen.node == s.node


def test_fcfs_serves_earliest_join_regardless_of_level(ed):
    model, network = ed
    rng = np.random.default_rng(5)
    for _ in range(1000):
        s = random_state(model, network, rng)
        earliest = min(s.candidates, key=lambda c: (c.hol_time, c.class_id))
        assert FCFS().select(s.node, s, rng) == earliest.class_id


def test_uncontested_decisions_never_change(ed):
    """If all feasible classes share one next node, the centrality term is a constant shift."""
    model, network = ed
    rng = np.random.default_rng(2)
    bp = PlainBP()
    for _ in range(1000):
        s = random_state(model, network, rng)
        F = feasible_set(s)
        if len(F) >= 2 and len({c.next_node for c in F}) >= 2:
            continue
        for theta in (0.1, 5.0, 100.0):
            assert KatzBP(theta, 0.5).select(s.node, s, rng) == bp.select(s.node, s, rng)


def test_one_hop_theta_zero_matches_plain_bp_and_respects_triage(ed):
    from src.policies import OneHopBP
    model, network = ed
    rng = np.random.default_rng(7)
    bp = PlainBP()
    for _ in range(1000):
        s = random_state(model, network, rng)
        assert OneHopBP(0.0, 0.5).select(s.node, s, rng) == bp.select(s.node, s, rng)
        r_star = min(c.level for c in s.candidates)
        assert model.buffers[OneHopBP(5.0, 0.5).select(s.node, s, rng)].level == r_star
