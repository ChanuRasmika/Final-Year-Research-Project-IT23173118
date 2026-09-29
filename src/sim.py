"""Event-driven SimPy simulation of the ED network.

Servers are not SimPy Resources: each node keeps a busy-server count and one
FIFO deque per buffer, and the acting policy picks the buffer at every
dispatch. Nodes in ``model.fcfs_nodes`` are served FCFS instead of by the
policy. A *decision epoch* is a server release at a policy node that finds
patients waiting; the optional observer sees each epoch (used for shadow
evaluation). The released patient is routed downstream before the decision.

Blocking-after-service (``model.blocking``): a node's waiting room holds at
most ``model.capacity[j]`` patients. A patient who finishes at i while the
next node j has no idle server and a full waiting room stays in i's server
(blocking it) until a place frees at j; blocked patients enter j in the order
they became blocked. A server is released, and a decision epoch occurs, only
when its patient actually leaves.

Shock (``shock``): services at one node that *start* inside [start, end) take
1 / rate_multiplier times as long (same random draw, scaled).
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Callable, NamedTuple

import numpy as np
import simpy

from .model import EDModel
from .network import Network
from .policies import FCFS, Candidate, DecisionState, Policy

# Stable stream ids keep common random numbers aligned across later experiments.
STREAM_ARRIVALS, STREAM_CLASS, STREAM_POLICY, STREAM_SHADOW = 0, 1, 2, 3
STREAM_SERVICE_BASE = 100


def make_rng(seed: int, stream: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([int(seed), int(stream)]))


class Patient:
    __slots__ = ("pid", "class_id", "t_arrival", "t_join")

    def __init__(self, pid: int, class_id: int, t_arrival: float):
        self.pid = pid
        self.class_id = class_id
        self.t_arrival = t_arrival
        self.t_join = t_arrival


class Shock(NamedTuple):
    node: int
    start: float
    end: float
    rate_multiplier: float   # e.g. 0.5 halves the service rate (doubles service times)


@dataclass
class SimResult:
    arrivals: int
    departures: int
    in_system: int
    epochs: int                       # post-warm-up decision epochs
    utilisation: dict[str, float]     # achieved, post-warm-up
    mean_in_system: dict[str, float]  # time-average waiting + in service, post-warm-up
    mean_waiting: dict[str, float]
    blocking_events: int = 0          # post-warm-up


Observer = Callable[["EDSimulation", DecisionState, int], None]


class EDSimulation:
    def __init__(self, model: EDModel, network: Network, policy: Policy, seed: int,
                 warmup: float, run_length: float, observer: Observer | None = None,
                 shock: Shock | None = None):
        if model.arrival_rate is None:
            raise ValueError("model has no arrival rate; call calibrate() first")
        if shock is not None and not (0 < shock.rate_multiplier and shock.start < shock.end):
            raise ValueError(f"invalid shock {shock}")
        self.model = model
        self.shock = shock
        self.network = network
        self.policy = policy
        self.fcfs = FCFS()
        self.observer = observer
        self.warmup = float(warmup)
        self.horizon = float(warmup + run_length)
        self.env = simpy.Environment()

        self.rng_arrivals = make_rng(seed, STREAM_ARRIVALS)
        self.rng_class = make_rng(seed, STREAM_CLASS)
        self.rng_policy = make_rng(seed, STREAM_POLICY)
        self.rng_service = [make_rng(seed, STREAM_SERVICE_BASE + i) for i in range(model.n_nodes)]

        n = model.n_nodes
        self.queues: list[deque[Patient]] = [deque() for _ in model.buffers]
        self.waiting = [0] * n
        self.busy = [0] * n
        self.arrivals = 0
        self.departures = 0
        self.epochs = 0
        self.blocking_events = 0
        self.capacity = list(model.capacity) if model.blocking else [np.inf] * n
        self.blocked: list[deque] = [deque() for _ in range(n)]   # per target node: (node, buffer, patient)

        self._last = [0.0] * n
        self._area_busy = [0.0] * n
        self._area_wait = [0.0] * n

        # Lognormal parameters from mean and CV: sigma^2 = ln(1 + cv^2), m = ln(mean) - sigma^2/2.
        self._svc = []
        for b in model.buffers:
            s2 = np.log1p(b.cv ** 2)
            self._svc.append((b.dist, b.mean, np.log(b.mean) - s2 / 2, np.sqrt(s2)))

    # ---------- bookkeeping ----------
    def _touch(self, node: int) -> None:
        """Accumulate time-weighted busy/waiting areas since the last change (post-warm-up only)."""
        now = self.env.now
        t0 = max(self._last[node], self.warmup)
        if now > t0:
            dt = now - t0
            self._area_busy[node] += self.busy[node] * dt
            self._area_wait[node] += self.waiting[node] * dt
        self._last[node] = now

    def _base_service_time(self, b: int, node: int) -> float:
        dist, mean, m, s = self._svc[b]
        rng = self.rng_service[node]
        return rng.exponential(mean) if dist == "exponential" else rng.lognormal(m, s)

    def shock_factor(self, node: int, t: float) -> float:
        """Service-time multiplier for a service starting at ``t`` (1 outside the shock)."""
        sh = self.shock
        if sh is not None and node == sh.node and sh.start <= t < sh.end:
            return 1.0 / sh.rate_multiplier
        return 1.0

    def _service_time(self, b: int, node: int) -> float:
        base = self._base_service_time(b, node)
        factor = self.shock_factor(node, self.env.now)
        return base if factor == 1.0 else base * factor

    def _is_full(self, node: int) -> bool:
        return self.busy[node] >= self.model.servers[node] and self.waiting[node] >= self.capacity[node]

    # ---------- events ----------
    def _arrival_process(self):
        probs = self.model.class_probs
        classes = self.model.classes
        pid = 0
        while True:
            yield self.env.timeout(self.rng_arrivals.exponential(1.0 / self.model.arrival_rate))
            cid = int(self.rng_class.choice(len(probs), p=probs))
            self.arrivals += 1
            pid += 1
            self._join(Patient(pid, cid, self.env.now), classes[cid].first_buffer)

    def _join(self, patient: Patient, b: int) -> None:
        node = self.model.buffers[b].node
        self._touch(node)
        patient.t_join = self.env.now
        self.queues[b].append(patient)
        self.waiting[node] += 1
        if self.busy[node] < self.model.servers[node]:
            # Idle server: the newcomer is the only one waiting, so there is no real choice.
            pol = self.fcfs if node in self.model.fcfs_nodes else self.policy
            self._start(node, pol.select(node, self.state(node), self.rng_policy))

    def _start(self, node: int, b: int) -> None:
        self._touch(node)
        patient = self.queues[b].popleft()
        self.waiting[node] -= 1
        self.busy[node] += 1
        done = self.env.timeout(self._service_time(b, node))
        done.callbacks.append(lambda _ev: self._complete(node, b, patient))
        if self.blocked[node]:
            self._unblock(node)   # a waiting place just freed at this node

    def _unblock(self, node: int) -> None:
        while self.blocked[node] and not self._is_full(node):
            i, b, patient = self.blocked[node].popleft()
            self._leave(i, b, patient)

    def _complete(self, node: int, b: int, patient: Patient) -> None:
        nxt = self.model.buffers[b].next_buffer
        if nxt is not None and self._is_full(self.model.buffers[nxt].node):
            # Blocking-after-service: the patient keeps this server until space frees downstream.
            self.blocked[self.model.buffers[nxt].node].append((node, b, patient))
            if self.env.now >= self.warmup:
                self.blocking_events += 1
            return
        self._leave(node, b, patient)

    def _leave(self, node: int, b: int, patient: Patient) -> None:
        """The patient leaves node's server: route downstream, then fill the freed server."""
        self._touch(node)
        self.busy[node] -= 1
        nxt = self.model.buffers[b].next_buffer
        if nxt is None:
            self.departures += 1
        else:
            self._join(patient, nxt)
        if self.waiting[node] > 0 and node in self.model.fcfs_nodes:
            self._start(node, self.fcfs.select(node, self.state(node), self.rng_policy))
        elif self.waiting[node] > 0:
            state = self.state(node)
            choice = self.policy.select(node, state, self.rng_policy)
            if self.env.now >= self.warmup:
                self.epochs += 1
            if self.observer is not None:
                self.observer(self, state, choice)
            self._start(node, choice)

    # ---------- state for policies ----------
    def state(self, node: int) -> DecisionState:
        buffers = self.model.buffers
        cands = []
        for b in self.model.node_buffers[node]:
            dq = self.queues[b]
            if dq:
                buf = buffers[b]
                nq = len(self.queues[buf.next_buffer]) if buf.next_buffer is not None else 0
                cands.append(Candidate(b, buf.level, len(dq), buf.next_node, nq, buf.mu, dq[0].t_join))
        return DecisionState(self.env.now, node, tuple(cands),
                             np.array(self.waiting, dtype=float), self.network)

    # ---------- run ----------
    def run(self) -> SimResult:
        self.env.process(self._arrival_process())
        self.env.run(until=self.horizon)
        for i in range(self.model.n_nodes):
            self._touch(i)
        span = self.horizon - self.warmup
        nodes = self.model.nodes
        servers = self.model.servers
        return SimResult(
            arrivals=self.arrivals,
            departures=self.departures,
            in_system=sum(self.waiting) + sum(self.busy),
            epochs=self.epochs,
            utilisation={n: self._area_busy[i] / (servers[i] * span) for i, n in enumerate(nodes)},
            mean_in_system={n: (self._area_busy[i] + self._area_wait[i]) / span for i, n in enumerate(nodes)},
            mean_waiting={n: self._area_wait[i] / span for i, n in enumerate(nodes)},
            blocking_events=self.blocking_events,
        )
