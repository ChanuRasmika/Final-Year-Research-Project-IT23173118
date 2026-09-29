"""Pilot figures (static PNGs, light surface)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .metrics import lookup  # noqa: E402
from .shadow import ThetaSpec  # noqa: E402

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]   # fixed categorical order (blue, orange, aqua)


def _style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=TEXT, fontsize=12, pad=10)
    ax.set_xlabel(xlabel, color=TEXT_2)
    ax.set_ylabel(ylabel, color=TEXT_2)
    ax.tick_params(colors=TEXT_2, length=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def r_change_vs_theta(summary: pd.DataFrame, kappas, thetas: list[ThetaSpec], vcfg: dict, path: Path):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    order = sorted(thetas, key=lambda t: t.value)
    x = np.array([t.value for t in order])
    for i, k in enumerate(kappas):
        rows = [lookup(summary, "r_change", k, t.label) for t in order]
        y = np.array([r["mean"] for r in rows]) * 100
        err = np.nan_to_num(np.array([r["ci_halfwidth"] for r in rows]) * 100)
        color = SERIES[i % len(SERIES)]
        ax.errorbar(x, y, yerr=err, color=color, linewidth=2, marker="o", markersize=5,
                    capsize=3, label=f"κ = {k:g}")
    for t in order:
        if t.kind == "eta":
            ax.axvline(t.value, color=GRID, linestyle=":", linewidth=1.2, zorder=0)
            ax.annotate(f"η={t.label[3:]}", (t.value, 0.0), xycoords=("data", "axes fraction"),
                        xytext=(2, 4), textcoords="offset points", color=TEXT_2, fontsize=8)
    for lvl, name in ((vcfg["r_change_min"], "5% floor"), (vcfg["r_change_max"], "50% override")):
        ax.axhline(lvl * 100, color=TEXT_2, linestyle="--", linewidth=0.9)
        ax.annotate(name, (x[-1], lvl * 100), xytext=(0, 3), textcoords="offset points",
                    ha="right", color=TEXT_2, fontsize=8)
    ax.set_xscale("log")
    _style(ax, "Share of contested decisions changed by the centrality term",
           "θ (log scale; dotted = pressure-scaled θ)", "R_change (%)  mean ± 95% CI")
    ax.set_ylim(0, max(55, ax.get_ylim()[1]))
    ax.legend(frameon=False, labelcolor=TEXT_2, loc="lower right", bbox_to_anchor=(1, 0.12))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def bp_gap_hist(epochs: pd.DataFrame, path: Path):
    gap = epochs.loc[epochs["contested"], "bp_top2_gap"].dropna()
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    if len(gap):
        hi = gap.quantile(0.995)
        ax.hist(gap[gap <= hi], bins=40, color=SERIES[0], edgecolor=SURFACE, linewidth=1)
        med = gap.median()
        ax.axvline(med, color=TEXT, linestyle="--", linewidth=1)
        ax.annotate(f"median {med:.3g}", (med, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(4, -14), textcoords="offset points", color=TEXT, fontsize=9)
        ties = (gap == 0).mean()
        ax.annotate(f"exact ties (gap = 0): {ties:.1%}\nx-axis clipped at 99.5th percentile",
                    (1.0, 0.85), xycoords="axes fraction", ha="right", color=TEXT_2, fontsize=9)
    _style(ax, "Top-2 BP pressure gap in contested epochs (all replications)",
           "|P₁ − P₂|  (1/min)", "epochs")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def contested_by_node(summary: pd.DataFrame, nodes: list[str], path: Path):
    rows = [lookup(summary, "contested_share", node=n) for n in nodes]
    y = np.nan_to_num(np.array([r["mean"] for r in rows])) * 100
    err = np.nan_to_num(np.array([r["ci_halfwidth"] for r in rows])) * 100
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    bars = ax.bar(nodes, y, yerr=err, color=SERIES[0], width=0.6, capsize=3,
                  error_kw={"ecolor": TEXT_2, "elinewidth": 1})
    for b, v, e in zip(bars, y, err):
        ax.annotate(f"{v:.1f}%", (b.get_x() + b.get_width() / 2, v + e), xytext=(0, 4),
                    textcoords="offset points", ha="center", color=TEXT, fontsize=9)
    _style(ax, "Contested decision epochs by node", "", "contested share (%)  mean ± 95% CI")
    ax.grid(axis="x", visible=False)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
