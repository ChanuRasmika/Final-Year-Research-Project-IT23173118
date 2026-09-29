from pathlib import Path

import numpy as np
import pytest

from src.pilot import load_config, setup
from src.policies import Candidate, DecisionState

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def cfg():
    return load_config(ROOT / "configs" / "pilot.yaml")


@pytest.fixture(scope="session")
def ed(cfg):
    return setup(cfg)


def random_state(model, network, rng) -> DecisionState:
    """A random but structurally valid decision state at a random node."""
    while True:
        node = int(rng.integers(model.n_nodes))
        cands = []
        for b in model.node_buffers[node]:
            q = int(rng.integers(0, 6)) if rng.random() < 0.6 else 0
            if q:
                buf = model.buffers[b]
                nq = 0 if buf.next_buffer is None else int(rng.integers(0, 6))
                cands.append(Candidate(b, buf.level, q, buf.next_node, nq, buf.mu, float(rng.uniform(0, 100))))
        if cands:
            q_vec = rng.integers(0, 30, size=model.n_nodes).astype(float)
            return DecisionState(0.0, node, tuple(cands), q_vec, network)
