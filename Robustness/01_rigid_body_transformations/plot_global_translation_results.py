# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "Global_translation"
INPUT_CSV = RESULT_DIR / "global_translation_results.csv"

ORDER = ["6MRR", "8HFE", "8VC8"]
COLORS = {
    "6MRR": "#0072B2",
    "8HFE": "#D55E00",
    "8VC8": "#009E73",
}
MARKERS = {
    "6MRR": "o",
    "8HFE": "s",
    "8VC8": "^",
}
MARKER_SIZES = {
    "6MRR": 7.2,
    "8HFE": 5.5,
    "8VC8": 4.0,
}


def load_rows() -> list[dict[str, object]]:
    rows = []
    with INPUT_CSV.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "case_id": row["case_id"],
                    "translation": float(row["translation_x_angstrom"]),
                    "bitacc": float(row["bitacc_percent"]),
                    "exact": row["exact_recovery"].lower() == "true",
                }
            )
    return rows


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    ax.tick_params(axis="both", which="major", direction="out", length=5, width=1.2, labelsize=8, pad=3)
    ax.tick_params(axis="both", which="minor", direction="out", length=3, width=1.0)


def main() -> None:
    rows = load_rows()
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )

    fig, ax = plt.subplots(figsize=(3.55, 2.55))
    for case_id in ORDER:
        group = sorted([row for row in rows if row["case_id"] == case_id], key=lambda row: row["translation"])
        x = [row["translation"] for row in group]
        y = [row["bitacc"] for row in group]
        ax.plot(
            x,
            y,
            color=COLORS[case_id],
            marker=MARKERS[case_id],
            markersize=MARKER_SIZES[case_id],
            markerfacecolor="white",
            markeredgewidth=1.1,
            linewidth=1.2,
            label=case_id,
            zorder=3,
        )

    ax.axhline(100, color="#A0A0A0", linewidth=0.9, linestyle=":", zorder=1)
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlim(0.08, 130)
    ax.set_ylim(94, 101.2)
    ax.set_xlabel("Global translation magnitude (Å)", fontsize=9, labelpad=4)
    ax.set_ylabel("Bit accuracy (%)", fontsize=9, labelpad=4)
    ax.set_xticks([0.1, 0.5, 1, 5, 10, 50, 100])
    ax.set_xticklabels(["0.1", "0.5", "1", "5", "10", "50", "100"])
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    style_axes(ax)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[case_id],
            marker=MARKERS[case_id],
            markersize=4.8,
            markerfacecolor="white",
            markeredgewidth=1.1,
            linewidth=1.2,
            label=case_id,
        )
        for case_id in ORDER
    ]
    legend = ax.legend(
        handles=legend_handles,
        loc="lower right",
        frameon=False,
        fontsize=8,
        handlelength=1.3,
        handletextpad=0.35,
        borderaxespad=0.2,
    )
    for text in legend.get_texts():
        text.set_fontweight("normal")

    fig.tight_layout(pad=0.6)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"global_translation_bitacc{suffix}", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
