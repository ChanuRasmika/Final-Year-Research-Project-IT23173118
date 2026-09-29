"""ED model: nodes, patient classes, per-stage buffers, and lambda calibration.

A patient class is c = (triage level k, pathway p). Because some pathways
revisit a node (consultation re-review), a class's queue at a node is split by
visit stage. Each (k, p, stage) queue is a *buffer*; it lives at exactly one
node and has exactly one next node, so n(i, c) is well defined. Policies choose
between buffers, and the buffer id is the ``class_id`` returned by ``select``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SINK = "sink"
DISTRIBUTIONS = ("lognormal", "exponential")


@dataclass(frozen=True)
class PatientClass:
    id: int
    level: int
    pathway: str
    prob: float
    first_buffer: int

    @property
    def name(self) -> str:
        return f"L{self.level}-{self.pathway}"


@dataclass(frozen=True)
class Buffer:
    id: int
    class_id: int
    level: int
    pathway: str
    stage: int              # position in the pathway (0 = first node)
    node: int
    next_node: int          # index into nodes, or model.sink for exit
    next_buffer: int | None  # None when the patient exits after this service
    mean: float             # E[S_i^c] in minutes
    cv: float
    dist: str

    @property
    def mu(self) -> float:
        return 1.0 / self.mean


@dataclass
class EDModel:
    nodes: list[str]
    servers: np.ndarray
    classes: list[PatientClass]
    buffers: list[Buffer]
    node_buffers: list[list[int]]
    blocking: bool = False
    fcfs_nodes: frozenset[int] = frozenset()   # served FCFS, outside the scheduling policy
    capacity: np.ndarray | None = None         # max waiting per node when blocking (inf = unbounded)
    arrival_rate: float | None = None
    target_utilisation: float | None = None
    utilisation: dict[str, float] = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def policy_nodes(self) -> list[int]:
        """Nodes where the scheduling policy decides (and decision epochs are counted)."""
        return [i for i in range(self.n_nodes) if i not in self.fcfs_nodes]

    @property
    def sink(self) -> int:
        """Index of the virtual sink (one past the last real node)."""
        return len(self.nodes)

    def node_name(self, idx: int) -> str:
        return SINK if idx == self.sink else self.nodes[idx]

    def buffer_name(self, b: int) -> str:
        buf = self.buffers[b]
        return f"L{buf.level}-{buf.pathway}-{self.nodes[buf.node]}#{buf.stage}"

    @property
    def class_probs(self) -> np.ndarray:
        return np.array([c.prob for c in self.classes])


def build_model(cfg: dict) -> EDModel:
    """Build the model from the ``model`` section of the config (lambda not yet set)."""
    nodes = list(cfg["nodes"])
    if SINK in nodes:
        raise ValueError("'sink' is virtual and must not be listed in model.nodes")
    idx = {n: i for i, n in enumerate(nodes)}
    sink = len(nodes)
    servers = np.array([int(cfg["servers"][n]) for n in nodes])
    if (servers < 1).any():
        raise ValueError("every node needs at least one server")

    level_mult = {int(k): float(v) for k, v in cfg.get("level_multiplier", {}).items()}
    revisit_mult = float(cfg.get("revisit_multiplier", 1.0))

    pathways = {}
    for name, route in cfg["pathways"].items():
        route = list(route)
        if route[-1] != SINK or SINK in route[:-1]:
            raise ValueError(f"pathway {name} must end with '{SINK}' exactly once")
        unknown = set(route[:-1]) - set(nodes)
        if unknown:
            raise ValueError(f"pathway {name} uses unknown nodes {unknown}")
        pathways[name] = route[:-1]

    mix = _class_mix(cfg["class_mix"], pathways)

    classes: list[PatientClass] = []
    buffers: list[Buffer] = []
    for cid, (level, pw, prob) in enumerate(mix):
        route = pathways[pw]
        first = len(buffers)
        classes.append(PatientClass(cid, level, pw, prob, first))
        for s, node in enumerate(route):
            svc = cfg["service"][node]
            dist = svc.get("dist", "lognormal")
            if dist not in DISTRIBUTIONS:
                raise ValueError(f"unknown service distribution {dist!r} at {node}")
            mean = float(svc["mean"]) * level_mult.get(level, 1.0)
            if node in route[:s]:
                mean *= revisit_mult
            last = s == len(route) - 1
            buffers.append(Buffer(
                id=len(buffers), class_id=cid, level=level, pathway=pw, stage=s,
                node=idx[node],
                next_node=sink if last else idx[route[s + 1]],
                next_buffer=None if last else len(buffers) + 1,
                mean=mean, cv=float(svc.get("cv", 1.0)), dist=dist,
            ))

    fcfs = cfg.get("fcfs_nodes", []) or []
    unknown = set(fcfs) - set(nodes)
    if unknown:
        raise ValueError(f"fcfs_nodes lists unknown nodes {unknown}")
    # Finite waiting rooms (max patients waiting, excluding those in service); used only when blocking.
    capacity = np.full(len(nodes), np.inf)
    for n, c in (cfg.get("capacity") or {}).items():
        if n not in idx:
            raise ValueError(f"capacity lists unknown node {n}")
        if int(c) < 0:
            raise ValueError(f"capacity at {n} must be >= 0")
        capacity[idx[n]] = int(c)
    node_buffers = [[b.id for b in buffers if b.node == i] for i in range(len(nodes))]
    return EDModel(nodes=nodes, servers=servers, classes=classes, buffers=buffers,
                   node_buffers=node_buffers, blocking=bool(cfg.get("blocking", False)),
                   fcfs_nodes=frozenset(idx[n] for n in fcfs), capacity=capacity)


def _class_mix(mix_cfg: dict, pathways: dict) -> list[tuple[int, str, float]]:
    """Joint Pr(k, p) = Pr(k) Pr(p | k), sorted by (level, pathway) for stable class ids."""
    level_probs = {int(k): float(v) for k, v in mix_cfg["level_probs"].items()}
    cond = {int(k): v for k, v in mix_cfg["pathway_given_level"].items()}
    if not np.isclose(sum(level_probs.values()), 1.0):
        raise ValueError("class_mix.level_probs must sum to 1")
    out = []
    for level in sorted(level_probs):
        if level not in range(1, 6):
            raise ValueError(f"triage level {level} outside 1..5")
        pw_probs = cond[level]
        if not np.isclose(sum(pw_probs.values()), 1.0):
            raise ValueError(f"class_mix.pathway_given_level[{level}] must sum to 1")
        for pw in sorted(pw_probs):
            if pw not in pathways:
                raise ValueError(f"class mix uses unknown pathway {pw}")
            p = level_probs[level] * float(pw_probs[pw])
            if p > 0:
                out.append((level, pw, p))
    return out


def load_per_arrival(model: EDModel) -> np.ndarray:
    """Per node: sum_c Pr(c) x visits_c(i) x E[S_i^c] / servers_i (utilisation per unit lambda)."""
    work = np.zeros(model.n_nodes)
    probs = model.class_probs
    for b in model.buffers:
        work[b.node] += probs[b.class_id] * b.mean
    return work / model.servers


def calibrate(model: EDModel, target: float) -> dict[str, float]:
    """Set lambda so the bottleneck node's analytic utilisation equals ``target``."""
    if not 0 < target < 1:
        raise ValueError("target utilisation must lie in (0, 1)")
    load = load_per_arrival(model)
    model.arrival_rate = target / load.max()
    model.target_utilisation = target
    model.utilisation = {n: float(model.arrival_rate * load[i]) for i, n in enumerate(model.nodes)}
    return model.utilisation
