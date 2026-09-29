"""Routing network and the bounded Katz congestion signal.

A[i, j] = Pr(next = j | at i), from the class mix and visit counts; the sink
column is dropped, so rows of nodes that discharge patients sum to < 1.
The Katz score solves (I - alpha A^T) c = q with alpha = kappa / rho(A); node
i's score accumulates its own queue plus the discounted queues upstream of it.
"""
from __future__ import annotations

import numpy as np

from .model import EDModel


def transition_matrix(model: EDModel) -> np.ndarray:
    """Weighted adjacency among real nodes (sink excluded)."""
    n = model.n_nodes
    visits = np.zeros(n)
    trans = np.zeros((n, n + 1))
    probs = model.class_probs
    for b in model.buffers:
        p = probs[b.class_id]
        visits[b.node] += p
        trans[b.node, b.next_node] += p
    A = np.zeros((n, n))
    seen = visits > 0
    A[seen] = trans[seen, :n] / visits[seen, None]
    return A


class Network:
    def __init__(self, A: np.ndarray, node_names: list[str]):
        A = np.asarray(A, dtype=float)
        if A.shape != (len(node_names), len(node_names)):
            raise ValueError("adjacency shape does not match node list")
        self.A = A
        self.nodes = list(node_names)
        self.n = len(node_names)
        self.rho = float(np.max(np.abs(np.linalg.eigvals(A))))
        if self.rho <= 1e-12:
            raise ValueError("spectral radius of A is 0 (acyclic routing): alpha = kappa / rho(A) "
                             "is undefined. Add re-review loops or change the Katz scaling.")
        self._M: dict[float, np.ndarray] = {}

    @classmethod
    def from_model(cls, model: EDModel) -> "Network":
        return cls(transition_matrix(model), model.nodes)

    def alpha(self, kappa: float) -> float:
        if not 0 < kappa < 1:
            raise ValueError("kappa must lie in (0, 1) so the Katz series converges")
        return kappa / self.rho

    def _system(self, kappa: float) -> np.ndarray:
        M = self._M.get(kappa)
        if M is None:
            M = np.eye(self.n) - self.alpha(kappa) * self.A.T
            self._M[kappa] = M
        return M

    def katz_raw(self, q: np.ndarray, kappa: float) -> np.ndarray:
        """Solve (I - alpha A^T) c = q (no explicit inverse)."""
        return np.linalg.solve(self._system(kappa), np.asarray(q, dtype=float))

    def one_hop_raw(self, q: np.ndarray, kappa: float) -> np.ndarray:
        """h = q + alpha A^T q: the first two terms of the Katz series (direct upstream only)."""
        q = np.asarray(q, dtype=float)
        return q + self.alpha(kappa) * (self.A.T @ q)

    def c_hat(self, q: np.ndarray, kappa: float) -> np.ndarray:
        """Bounded Katz score c_i / (1 + max_k |c_k|), with the sink appended as 0.

        Returned array has length n + 1 so ``c_hat[model.sink] == 0``.
        """
        return self._bounded(self.katz_raw(q, kappa))

    def h_hat(self, q: np.ndarray, kappa: float) -> np.ndarray:
        """Bounded one-hop score h_i / (1 + max_k |h_k|), with the sink appended as 0."""
        return self._bounded(self.one_hop_raw(q, kappa))

    def _bounded(self, c: np.ndarray) -> np.ndarray:
        c_hat = c / (1.0 + np.max(np.abs(c)))
        if c_hat.min() < -1e-12 or c_hat.max() >= 1.0:
            raise AssertionError(f"bounded score outside [0, 1): {c_hat}")
        out = np.empty(self.n + 1)
        out[:-1] = np.clip(c_hat, 0.0, None)   # remove round-off negatives only
        out[-1] = 0.0
        return out
