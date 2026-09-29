"""Scheduling policies behind a single interface: ``select(node, state, rng) -> class_id``.

Every BP-type policy first restricts to the most urgent triage level waiting at
the node (never compares across levels), then ranks the feasible buffers.
Ties break on earliest head-of-line join time, then on class (buffer) id.
``FCFS`` is used at the nodes listed in ``model.fcfs_nodes``.
``OneHopBP`` truncates the Katz series after one hop.
Future policies (APQ) subclass ``Policy`` the same way.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .network import Network


@dataclass(frozen=True)
class Candidate:
    """A non-empty buffer waiting at the decision node."""
    class_id: int       # buffer id
    level: int          # triage level (1 = most urgent)
    queue: int          # Q_i^c (waiting, excluding those in service)
    next_node: int      # j = n(i, c); the sink index when exiting
    next_queue: int     # Q_j^c (0 for the sink)
    mu: float           # 1 / E[S_i^c]
    hol_time: float     # time the head-of-line patient joined this buffer


@dataclass(frozen=True)
class DecisionState:
    time: float
    node: int
    candidates: tuple[Candidate, ...]
    q: np.ndarray       # aggregate waiting queue per real node, q_i = sum_c Q_i^c
    network: Network


def feasible_set(state: DecisionState) -> list[Candidate]:
    """F = waiting buffers at the most urgent triage level present."""
    if not state.candidates:
        raise ValueError("decision requested at a node with no waiting patients")
    r_star = min(c.level for c in state.candidates)
    return [c for c in state.candidates if c.level == r_star]


def bp_pressure(c: Candidate) -> float:
    return c.mu * (c.queue - c.next_queue)


def rank(F: list[Candidate], weights: list[float]) -> list[int]:
    """Indices of F from best to worst: highest weight, earliest HOL, lowest id."""
    return sorted(range(len(F)), key=lambda k: (-weights[k], F[k].hol_time, F[k].class_id))


class Policy(ABC):
    name = "policy"

    @abstractmethod
    def select(self, node: int, state: DecisionState, rng: np.random.Generator) -> int:
        """Return the class (buffer) id to serve next at ``node``."""


class FCFS(Policy):
    """Serve the patient who joined the node's queue earliest, ignoring triage level."""
    name = "FCFS"

    def select(self, node, state, rng):
        if not state.candidates:
            raise ValueError("decision requested at a node with no waiting patients")
        return min(state.candidates, key=lambda c: (c.hol_time, c.class_id)).class_id


class WeightedPolicy(Policy):
    """Triage filter + argmax of a per-candidate weight."""

    def weights(self, node: int, state: DecisionState, F: list[Candidate]) -> list[float]:
        raise NotImplementedError

    def select(self, node, state, rng):
        F = feasible_set(state)
        if len(F) == 1:
            return F[0].class_id
        return F[rank(F, self.weights(node, state, F))[0]].class_id


class PlainBP(WeightedPolicy):
    name = "BP"

    def weights(self, node, state, F):
        return [bp_pressure(c) for c in F]


class KatzBP(WeightedPolicy):
    """W_c = mu_i^c (Q_i^c - Q_j^c) + theta (c_hat_i - c_hat_j)."""

    def __init__(self, theta: float, kappa: float):
        if theta < 0:
            raise ValueError("theta must be non-negative")
        self.theta = float(theta)
        self.kappa = float(kappa)
        self.name = f"KatzBP(theta={theta:g}, kappa={kappa:g})"

    def weights(self, node, state, F):
        c_hat = state.network.c_hat(state.q, self.kappa)
        return [bp_pressure(c) + self.theta * (c_hat[node] - c_hat[c.next_node]) for c in F]


class OneHopBP(WeightedPolicy):
    """1H-BP: W_c = mu_i^c (Q_i^c - Q_j^c) + theta_h (h_hat_i - h_hat_j), h = q + alpha A^T q."""

    def __init__(self, theta_h: float, kappa: float):
        if theta_h < 0:
            raise ValueError("theta_h must be non-negative")
        self.theta = float(theta_h)
        self.kappa = float(kappa)
        self.name = f"OneHopBP(theta_h={theta_h:g}, kappa={kappa:g})"

    def weights(self, node, state, F):
        h_hat = state.network.h_hat(state.q, self.kappa)
        return [bp_pressure(c) + self.theta * (h_hat[node] - h_hat[c.next_node]) for c in F]
