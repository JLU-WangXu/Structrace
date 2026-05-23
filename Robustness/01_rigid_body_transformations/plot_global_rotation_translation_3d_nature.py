# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


ROOT = Path(__file__).resolve().parent
INPUT_CSV = ROOT / "global_rotation_translation_results.csv"
ORDER = ["6MRR", "8HFE", "8VC8"]
OUTPUT_PREFIX = "global_rotation_translation_3d_nature"


def load_rows() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with INPUT_CSV.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "case_id": row["case_id"],
                    "angle_x": float(row["angle_x_deg"]),
                    "angle_z": float(row["angle_z_deg"]),
                    "translation_x": float(row["translation_x_angstrom"]),
                    "original_rmsd": float(row["original_alignment_ca_rmsd"]),
                    "watermarked_rmsd": float(row["watermarked_alignment_ca_rmsd"]),
                    "score": np.log10(
                        float(row["watermarked_alignment_ca_rmsd"]) / float(row["original_alignment_ca_rmsd"])
                    ),
                    "bitacc": float(row["bitacc_percent"]),
                }
            )
    return rows


def grouped_grid(rows: list[dict[str, float | str]], case_id: str):
    case_rows = [row for row in rows if row["case_id"] == case_id]
    angle_x = sorted({float(row["angle_x"]) for row in case_rows})
    angle_z = sorted({float(row["angle_z"]) for row in case_rows})
    translation_x = sorted({float(row["translation_x"]) for row in case_rows})

    lookup = {
        (float(row["angle_x"]), float(row["angle_z"]), float(row["translation_x"])): float(row["score"])
        for row in case_rows
    }

    for tx in translation_x:
        grid = np.array(
            [[lookup[(ax, az, tx)] for az in angle_z] for ax in angle_x],
            dtype=float,
        )
        yield tx, np.array(angle_x, dtype=float), np.array(angle_z, dtype=float), grid


def z_transform(values):
    return np.log10(np.asarray(values, dtype=float) + 1.0)


def style_3d_axis(ax) -> None:
    ax.grid(True)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.set_facecolor((1, 1, 1, 0))
        pane.set_edgecolor((0.72, 0.72, 0.72, 1))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = (0.82, 0.82, 0.82, 1)
        axis._axinfo["grid"]["linewidth"] = 0.6
        axis._axinfo["tick"]["inward_factor"] = 0.0
        axis._axinfo["tick"]["outward_factor"] = 0.25
    ax.tick_params(axis="both", which="major", labelsize=8, pad=1)


def build_figure(rows: list[dict[str, float | str]], case_id: str, vmin: float, vmax: float) -> plt.Figure:
    case_rows = [row for row in rows if row["case_id"] == case_id]
    bitacc = int(round(float(case_rows[0]["bitacc"]))) if case_rows else 0
    fig = plt.figure(figsize=(5.9, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=vmin, vmax=vmax)

    for tx, angle_x, angle_z, grid in grouped_grid(rows, case_id):
        x, y = np.meshgrid(angle_z, angle_x)
        z = np.full_like(x, z_transform([tx])[0], dtype=float)
        face = grid[:-1, :-1]
        ax.plot_surface(
            x,
            y,
            z,
            facecolors=cmap(norm(face)),
            rstride=1,
            cstride=1,
            linewidth=0.12,
            antialiased=False,
            shade=False,
            alpha=0.92,
        )

    ax.view_init(elev=23, azim=-58)
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    translation_ticks = [0, 0.5, 1, 5, 10, 50, 100]
    ax.set_zlim(z_transform([0])[0], z_transform([100])[0])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 90, 180, 270, 360])
    ax.set_zticks(z_transform(translation_ticks))
    ax.set_zticklabels(["0", "0.5", "1", "5", "10", "50", "100"])
    ax.set_xlabel("Rotation about z axis (deg)", labelpad=5, fontsize=8)
    ax.set_ylabel("Rotation about x axis (deg)", labelpad=5, fontsize=8)
    ax.set_zlabel("Translation x (A)", labelpad=5, fontsize=8)
    ax.set_title(case_id, fontsize=11, pad=10, weight="bold")
    style_3d_axis(ax)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.03, fraction=0.04, shrink=0.83)
    cbar.set_label("log10(watermarked / original RMSD)", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5, width=0.8, length=4)

    ax.text2D(
        0.03,
        0.97,
        f"BitAcc = {bitacc:.0f}%  |  n = {len(case_rows)}",
        transform=ax.transAxes,
        fontsize=8,
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=2.2),
    )
    ax.text2D(
        0.03,
        0.02,
        "Static 3D response map",
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        ha="left",
        color="#444444",
    )

    fig.subplots_adjust(left=0.03, right=0.86, bottom=0.06, top=0.94)
    return fig


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(ROOT / f"{stem}.{ext}", bbox_inches="tight")


def main() -> None:
    rows = load_rows()
    global_values = [float(row["score"]) for row in rows]
    vmin = min(global_values)
    vmax = max(global_values)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 600,
        }
    )

    for case_id in ORDER:
        fig = build_figure(rows, case_id, vmin, vmax)
        save_figure(fig, f"{OUTPUT_PREFIX}_{case_id}")
        plt.close(fig)


if __name__ == "__main__":
    main()
