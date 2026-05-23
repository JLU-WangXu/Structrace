# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "Global_rotation_translation"
INPUT_CSV = RESULT_DIR / "global_rotation_translation_results.csv"


def load_summary_points() -> list[dict[str, float]]:
    grouped: dict[tuple[float, float, float], list[float]] = defaultdict(list)
    with INPUT_CSV.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                float(row["angle_x_deg"]),
                float(row["angle_z_deg"]),
                float(row["translation_x_angstrom"]),
            )
            grouped[key].append(float(row["bitacc_percent"]))
    points = []
    for (angle_x, angle_z, translation), values in grouped.items():
        points.append(
            {
                "angle_x": angle_x,
                "angle_z": angle_z,
                "translation": translation,
                "bitacc": float(np.mean(values)),
                "spread": float(np.max(values) - np.min(values)),
            }
        )
    return points


def compressed_translation(distance: float) -> float:
    return float(np.log10(distance + 1.0))


def main() -> None:
    points = load_summary_points()
    values = np.array([point["bitacc"] for point in points], dtype=float)
    vmin = min(94.0, float(values.min()))
    vmax = 100.0
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")

    x = np.array([point["angle_x"] for point in points], dtype=float)
    y = np.array([point["angle_z"] for point in points], dtype=float)
    z = np.array([compressed_translation(point["translation"]) for point in points], dtype=float)
    c = np.array([point["bitacc"] for point in points], dtype=float)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )

    fig = plt.figure(figsize=(4.25, 3.35))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        x,
        y,
        z,
        c=c,
        cmap=cmap,
        norm=norm,
        s=3.2,
        marker="o",
        edgecolors="none",
        alpha=0.72,
        depthshade=False,
    )

    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    ax.set_zlim(compressed_translation(0), compressed_translation(100))
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 90, 180, 270, 360])
    translation_ticks = [0, 1, 10, 100]
    ax.set_zticks([compressed_translation(value) for value in translation_ticks])
    ax.set_zticklabels([str(value) for value in translation_ticks])
    ax.set_xlabel("Rotation about x axis (deg)", fontsize=6.8, labelpad=4)
    ax.set_ylabel("Rotation about z axis (deg)", fontsize=6.8, labelpad=5)
    ax.set_zlabel("Translation (Å)", fontsize=6.8, labelpad=4)
    ax.view_init(elev=22, azim=-50)
    ax.tick_params(axis="both", which="major", labelsize=5.8, pad=0)
    ax.set_box_aspect((1.1, 1.0, 0.62))

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.86, 0.86, 0.86, 1.0))
    ax.grid(True, linewidth=0.38, color="#D8D8D8")

    scalar = ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar = fig.colorbar(
        scalar,
        ax=ax,
        fraction=0.04,
        pad=0.09,
        shrink=0.66,
        aspect=22,
    )
    colorbar.set_label("Bit accuracy (%)", fontsize=6.8, labelpad=4)
    colorbar.ax.tick_params(labelsize=5.8, length=2.4, width=0.7)

    fig.subplots_adjust(left=0.01, right=0.82, bottom=0.03, top=0.99)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(
            RESULT_DIR / f"global_rotation_translation_heatmap_3d{suffix}",
            bbox_inches="tight",
            pad_inches=0.18,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
