# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator


RESULT_DIR = Path(__file__).resolve().parent
INPUT_CSV = RESULT_DIR / "coordinate_decimal_rounding_results.csv"

CASES = ["6MRR", "8HFE", "8VC8"]
DECIMALS = [3, 2, 1, 0]
COLORS = {
    "6MRR": "#0072B2",
    "8HFE": "#D55E00",
    "8VC8": "#009E73",
}


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "case_id": row["case_id"],
                    "decimals": int(row["decimals"]),
                    "bitacc": float(row["bitacc_percent"]),
                    "ber": float(row["ber_percent"]),
                    "exact": row["exact_recovery"].lower() == "true",
                    "ca_rmsd": float(row["clean_to_rounded_ca_rmsd_angstrom"]),
                    "max_ca_displacement": float(row["max_ca_rounding_displacement_angstrom"]),
                }
            )
    return rows


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="major", direction="out", length=4.2, width=1.0, labelsize=7.4, pad=2.5)
    ax.tick_params(axis="both", which="minor", direction="out", length=2.5, width=0.8)


def row_value(rows: list[dict[str, object]], case_id: str, decimals: int, key: str) -> float:
    row = next(row for row in rows if row["case_id"] == case_id and row["decimals"] == decimals)
    if key == "exact":
        return 100.0 if row["exact"] else 0.0
    return float(row[key])


def main() -> None:
    rows = load_rows()
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.linewidth": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )

    x_positions = np.arange(len(DECIMALS))
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.45), sharex=True)
    panels = [
        ("bitacc", "Bit accuracy (%)", (0, 105)),
        ("exact", "Exact recovery (%)", (0, 105)),
        ("ca_rmsd", "C-alpha RMSD (A)", (0, None)),
    ]

    for ax, (metric, ylabel, ylim) in zip(axes, panels):
        for case_id in CASES:
            values = [row_value(rows, case_id, decimals, metric) for decimals in DECIMALS]
            ax.plot(
                x_positions,
                values,
                color=COLORS[case_id],
                marker="o",
                markersize=4.4,
                markerfacecolor="white",
                markeredgewidth=1.0,
                linewidth=1.35,
                label=case_id,
                zorder=3,
            )
        ax.set_xticks(x_positions)
        ax.set_xticklabels([str(decimals) for decimals in DECIMALS])
        ax.set_xlabel("Coordinate precision (decimals)", fontsize=7.8, labelpad=3.5)
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=3.5)
        ax.set_ylim(*ylim) if ylim[1] is not None else ax.set_ylim(bottom=0)
        if metric in {"bitacc", "exact"}:
            ax.axhline(100, color="#A0A0A0", linewidth=0.75, linestyle=":", zorder=0)
            ax.yaxis.set_major_locator(MultipleLocator(20))
            ax.yaxis.set_minor_locator(MultipleLocator(10))
        style_axes(ax)

    axes[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=1.5)
    fig.tight_layout(w_pad=1.05, pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"coordinate_decimal_rounding_summary_curves{suffix}", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.35, 2.75))
    steps = np.array([10 ** (-decimals) for decimals in DECIMALS], dtype=float)
    xlabels = [f"{decimals} dec\n{step:g} A" for decimals, step in zip(DECIMALS, steps)]

    ax.axhline(100, color="#202020", linewidth=0.9, linestyle="--", dashes=(3.2, 2.0), zorder=1)

    for case_id in CASES:
        values = [row_value(rows, case_id, decimals, "bitacc") for decimals in DECIMALS]
        ax.plot(
            steps,
            values,
            color=COLORS[case_id],
            linewidth=1.55,
            marker="o",
            markersize=5.2,
            markerfacecolor="white",
            markeredgewidth=1.2,
            label=case_id,
            zorder=3,
        )

    ax.set_xscale("log")
    ax.set_xlim(steps[0] * 0.68, steps[-1] * 1.48)
    ax.set_ylim(44, 102.5)
    ax.set_xticks(steps)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Rounding step / coordinate precision", fontsize=7.8, labelpad=4.5)
    ax.set_ylabel("Bit accuracy (%)", fontsize=7.8, labelpad=3.5)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))

    case_handles = [
        Line2D([0], [0], color=COLORS[case_id], lw=1.55, marker="o", markersize=4.8, markerfacecolor="white")
        for case_id in CASES
    ]
    ax.legend(case_handles, CASES, frameon=False, fontsize=6.8, loc="lower left", handlelength=1.6)
    style_axes(ax)
    fig.tight_layout(pad=0.65)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"coordinate_decimal_rounding_bitacc_bar{suffix}", bbox_inches="tight")
        fig.savefig(RESULT_DIR / f"coordinate_decimal_rounding_bitacc_response{suffix}", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
