import numpy as np
import pytest

from src.model import load_per_arrival
from src.network import Network


def test_c_hat_bounded(ed):
    model, network = ed
    rng = np.random.default_rng(3)
    for _ in range(1000):
        q = rng.integers(0, 200, size=network.n) * (rng.random(network.n) < 0.7)
        for k in (0.3, 0.5, 0.7):
            c_hat = network.c_hat(q, k)
            assert c_hat.shape == (network.n + 1,)
            assert np.all(c_hat >= 0) and np.all(c_hat < 1)
            assert c_hat[model.sink] == 0


def test_katz_solve_matches_truncated_series(ed):
    _, network = ed
    rng = np.random.default_rng(4)
    for _ in range(100):
        q = rng.uniform(0, 50, size=network.n)
        for k in (0.3, 0.5, 0.7):
            alpha, At = network.alpha(k), network.A.T
            series, term = np.zeros(network.n), q.copy()
            for _ in range(500):
                series += term
                term = alpha * At @ term
            assert np.allclose(network.katz_raw(q, k), series, rtol=0, atol=1e-8)


def test_zero_spectral_radius_fails_loudly():
    A = np.array([[0.0, 1.0], [0.0, 0.0]])  # acyclic: rho(A) = 0
    with pytest.raises(ValueError, match="spectral radius"):
        Network(A, ["a", "b"])


def test_adjacency_is_substochastic(ed):
    model, network = ed
    rows = network.A.sum(axis=1)
    assert np.all(rows <= 1 + 1e-12) and np.all(network.A >= 0)
    reg, tri = model.nodes.index("registration"), model.nodes.index("triage")
    assert network.A[reg, tri] == pytest.approx(1.0)


def test_calibration_hits_target(ed, cfg):
    model, _ = ed
    util = model.arrival_rate * load_per_arrival(model)
    assert util.max() == pytest.approx(cfg["model"]["target_utilisation"])


def test_one_hop_equals_first_two_katz_terms(ed):
    """h = q + alpha A^T q is the Katz series truncated after n = 1."""
    _, network = ed
    rng = np.random.default_rng(6)
    for _ in range(100):
        q = rng.uniform(0, 50, size=network.n)
        for k in (0.3, 0.5, 0.7):
            alpha, At = network.alpha(k), network.A.T
            two_terms = sum(alpha ** n * np.linalg.matrix_power(At, n) @ q for n in range(2))
            assert np.allclose(network.one_hop_raw(q, k), two_terms, rtol=0, atol=1e-12)
            h_hat = network.h_hat(q, k)
            assert np.all(h_hat >= 0) and np.all(h_hat < 1) and h_hat[-1] == 0
