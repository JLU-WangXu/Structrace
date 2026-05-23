# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


RESULT_DIR = Path(__file__).resolve().parent
INPUT_CSV = RESULT_DIR / "local_distortion_sigma_results.csv"

CASES = ["6MRR", "8HFE", "8VC8"]
SIGMAS = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
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
                    "sigma": float(row["sigma_angstrom"]),
                    "bitacc": float(row["bitacc_percent"]),
                    "global_ca_rmsd": float(row["global_ca_rmsd_angstrom"]),
                    "local_ca_rmsd": float(row["local_ca_rmsd_angstrom"]),
                    "exact": row["exact_recovery"].lower() == "true",
                }
            )
    return rows


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="major", direction="out", length=4.2, width=1.0, labelsize=7.4, pad=2.5)
    ax.tick_params(axis="both", which="minor", direction="out", length=2.5, width=0.8)


def gaussian_kde_manual(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    if len(values) == 0:
        return np.zeros_like(grid)
    diff = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * diff**2).sum(axis=1)
    density /= len(values) * bandwidth * np.sqrt(2 * np.pi)
    return density


def draw_raincloud_panel(
    ax,
    rows: list[dict[str, object]],
    case_id: str,
    ylabel: str | None = None,
    spacious: bool = False,
) -> None:
    color = COLORS[case_id]
    spacing = 1.48 if spacious else 1.0
    positions = np.arange(1, len(SIGMAS) + 1) * spacing
    grid = np.linspace(0, 100, 320)
    rng = np.random.default_rng(2000 + CASES.index(case_id))
    scatter_left = -0.44 if spacious else -0.34
    scatter_right = -0.25 if spacious else -0.17
    box_offset = 0.02 if spacious else 0.03
    box_width = 0.16 if spacious else 0.18
    density_start = 0.28 if spacious else 0.22
    density_width = 0.26 if spacious else 0.32
    point_size = 8.0 if spacious else 11.0

    for position, sigma in zip(positions, SIGMAS):
        values = np.array(
            [row["bitacc"] for row in rows if row["case_id"] == case_id and np.isclose(row["sigma"], sigma)],
            dtype=float,
        )
        density = gaussian_kde_manual(values, grid, bandwidth=3.2)
        if density.max() > 0:
            density = density / density.max() * density_width
        ax.fill_betweenx(
            grid,
            position + density_start,
            position + density_start + density,
            color=color,
            alpha=0.20,
            linewidth=0,
            zorder=1,
        )

        jitter = rng.uniform(scatter_left, scatter_right, size=len(values))
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=point_size,
            color=color,
            alpha=0.46,
            linewidths=0,
            zorder=2,
        )

        ax.boxplot(
            [values],
            positions=[position + box_offset],
            widths=box_width,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.05, "zorder": 5},
            boxprops={"linewidth": 0.95, "color": "black", "facecolor": color, "alpha": 0.55, "zorder": 4},
            whiskerprops={"linewidth": 0.95, "color": "black", "zorder": 4},
            capprops={"linewidth": 0.95, "color": "black", "zorder": 4},
        )

        mean_value = float(np.mean(values))
        ax.scatter(
            [position + box_offset],
            [mean_value],
            marker="D",
            s=10 if spacious else 12,
            color="white",
            edgecolor="black",
            linewidth=0.7,
            zorder=6,
        )

    ax.set_title(case_id, fontsize=8.5, fontweight="bold", pad=3)
    ax.set_xlim(positions[0] - 0.72, positions[-1] + 0.82)
    ax.set_ylim(0, 102)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{sigma:g}" for sigma in SIGMAS])
    ax.set_xlabel("Gaussian noise sigma (A)", fontsize=7.8, labelpad=3.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=3.5)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.axhline(100, color="#A0A0A0", linewidth=0.75, linestyle=":", zorder=0)
    style_axes(ax)


def draw_box_panel(ax, rows: list[dict[str, object]], case_id: str, metric: str, ylabel: str | None = None) -> None:
    data = [
        [row[metric] for row in rows if row["case_id"] == case_id and np.isclose(row["sigma"], sigma)]
        for sigma in SIGMAS
    ]
    positions = np.arange(1, len(SIGMAS) + 1)
    color = COLORS[case_id]
    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.15, "zorder": 4},
        boxprops={"linewidth": 1.0, "color": color, "zorder": 3},
        whiskerprops={"linewidth": 1.0, "color": color, "zorder": 3},
        capprops={"linewidth": 1.0, "color": color, "zorder": 3},
    )
    for patch in box["boxes"]:
        patch.set_facecolor(color)
        patch.set_alpha(0.18)

    rng = np.random.default_rng(1000 + CASES.index(case_id))
    for index, values in enumerate(data, start=1):
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=6,
            color=color,
            alpha=0.28,
            linewidths=0,
            zorder=2,
        )

    ax.set_title(case_id, fontsize=8.5, fontweight="bold", pad=3)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{sigma:g}" for sigma in SIGMAS], rotation=0)
    ax.set_xlabel("Gaussian noise σ (Å)", fontsize=7.8, labelpad=3.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=3.5)
    style_axes(ax)


def draw_density_overlay_panel(ax, rows: list[dict[str, object]], case_id: str, ylabel: str | None = None) -> float:
    color = COLORS[case_id]
    grid = np.linspace(0, 100, 760)
    alphas = np.linspace(0.22, 0.78, len(SIGMAS))
    max_density = 0.0

    for sigma, alpha in zip(SIGMAS, alphas):
        values = np.array(
            [row["bitacc"] for row in rows if row["case_id"] == case_id and np.isclose(row["sigma"], sigma)],
            dtype=float,
        )
        density = gaussian_kde_manual(values, grid, bandwidth=2.0)
        max_density = max(max_density, float(density.max()))
        ax.fill_between(
            grid,
            0,
            density,
            color=color,
            alpha=alpha * 0.30,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            grid,
            density,
            color=color,
            alpha=alpha,
            linewidth=1.35,
            zorder=2,
        )

    ax.axvline(100, color="#202020", linewidth=1.45, linestyle="--", dashes=(3.2, 2.0), zorder=4)
    ax.set_title(case_id, fontsize=8.8, fontweight="bold", pad=3)
    ax.set_xlim(30, 100)
    ax.set_ylim(0, max_density * 1.14 if max_density > 0 else 1)
    ax.set_xlabel("Bit accuracy (%)", fontsize=7.8, labelpad=3.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=3.5)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    style_axes(ax)
    return max_density


def draw_sigma_alpha_bar(bar_ax) -> None:
    alphas = np.linspace(0.22, 0.78, len(SIGMAS))
    y_positions = np.arange(len(SIGMAS))
    bar_ax.barh(
        y_positions,
        np.ones(len(SIGMAS)),
        height=0.82,
        color="#4A4A4A",
        alpha=1.0,
        edgecolor="none",
    )
    for patch, alpha in zip(bar_ax.patches, alphas):
        patch.set_alpha(alpha)
    bar_ax.set_xlim(0, 1)
    bar_ax.set_ylim(-0.55, len(SIGMAS) - 0.45)
    bar_ax.set_xticks([])
    bar_ax.set_yticks(y_positions)
    bar_ax.set_yticklabels([f"{sigma:g}" for sigma in SIGMAS], fontsize=6.8)
    bar_ax.yaxis.tick_right()
    bar_ax.tick_params(axis="y", direction="out", length=2.4, width=0.75, pad=2.0)
    bar_ax.set_title("sigma (A)", fontsize=6.9, pad=3)
    for spine in bar_ax.spines.values():
        spine.set_visible(False)


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

    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.35), sharey=True)
    for ax, case_id in zip(axes, CASES):
        draw_box_panel(
            ax,
            rows,
            case_id,
            "bitacc",
            "Bit accuracy (%)" if case_id == CASES[0] else None,
        )
        ax.set_ylim(0, 102)
        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(10))
        ax.axhline(100, color="#A0A0A0", linewidth=0.75, linestyle=":", zorder=1)

    fig.tight_layout(w_pad=1.2, pad=0.5)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"local_distortion_sigma_bitacc_boxplot{suffix}", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.35), sharey=True)
    for ax, case_id in zip(axes, CASES):
        draw_box_panel(
            ax,
            rows,
            case_id,
            "global_ca_rmsd",
            "Global Cα-RMSD (Å)" if case_id == CASES[0] else None,
        )
        ax.set_ylim(bottom=0)
    fig.tight_layout(w_pad=1.2, pad=0.5)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"local_distortion_sigma_global_rmsd_boxplot{suffix}", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(6.9, 4.15),
        sharex="col",
        gridspec_kw={"height_ratios": [2.15, 1.0]},
    )
    positions = np.arange(1, len(SIGMAS) + 1)
    for col, case_id in enumerate(CASES):
        ax_top = axes[0, col]
        draw_box_panel(
            ax_top,
            rows,
            case_id,
            "bitacc",
            "Bit accuracy (%)" if col == 0 else None,
        )
        ax_top.set_xlabel("")
        ax_top.set_ylim(0, 102)
        ax_top.yaxis.set_major_locator(MultipleLocator(20))
        ax_top.yaxis.set_minor_locator(MultipleLocator(10))
        ax_top.axhline(100, color="#A0A0A0", linewidth=0.75, linestyle=":", zorder=1)

        ax_bottom = axes[1, col]
        rates = []
        for sigma in SIGMAS:
            group = [
                row
                for row in rows
                if row["case_id"] == case_id and np.isclose(row["sigma"], sigma)
            ]
            rates.append(sum(bool(row["exact"]) for row in group) / len(group) * 100.0)
        color = COLORS[case_id]
        ax_bottom.plot(
            positions,
            rates,
            color=color,
            marker="o",
            markersize=4.2,
            markerfacecolor="white",
            markeredgewidth=1.0,
            linewidth=1.35,
        )
        ax_bottom.set_ylim(0, 105)
        ax_bottom.yaxis.set_major_locator(MultipleLocator(50))
        ax_bottom.yaxis.set_minor_locator(MultipleLocator(25))
        ax_bottom.set_xticks(positions)
        ax_bottom.set_xticklabels([f"{sigma:g}" for sigma in SIGMAS])
        ax_bottom.set_xlabel("Gaussian noise σ (Å)", fontsize=7.8, labelpad=3.5)
        if col == 0:
            ax_bottom.set_ylabel("Exact recovery (%)", fontsize=7.8, labelpad=3.5)
        style_axes(ax_bottom)

    fig.tight_layout(w_pad=1.2, h_pad=0.85, pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"local_distortion_sigma_bitacc_exact_combined{suffix}", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(7.05, 2.55), sharey=True)
    for ax, case_id in zip(axes, CASES):
        draw_raincloud_panel(
            ax,
            rows,
            case_id,
            "Bit accuracy (%)" if case_id == CASES[0] else None,
        )
    fig.tight_layout(w_pad=1.25, pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"local_distortion_sigma_bitacc_raincloud{suffix}", bbox_inches="tight")
    plt.close(fig)

    for case_id in CASES:
        fig, ax = plt.subplots(figsize=(7.6, 3.15))
        draw_raincloud_panel(ax, rows, case_id, "Bit accuracy (%)", spacious=True)
        ax.set_title(case_id, fontsize=9.0, fontweight="bold", pad=4)
        fig.tight_layout(pad=0.65)
        for suffix in [".png", ".pdf", ".svg"]:
            fig.savefig(
                RESULT_DIR / f"local_distortion_sigma_bitacc_raincloud_{case_id}{suffix}",
                bbox_inches="tight",
            )
        plt.close(fig)

    fig = plt.figure(figsize=(7.65, 2.35))
    grid_spec = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.075], wspace=0.34)
    axes = [fig.add_subplot(grid_spec[0, index]) for index in range(3)]
    axes[1].sharey(axes[0])
    axes[2].sharey(axes[0])
    bar_ax = fig.add_subplot(grid_spec[0, 3])
    density_maxima = []
    for ax, case_id in zip(axes, CASES):
        density_maxima.append(
            draw_density_overlay_panel(
                ax,
                rows,
                case_id,
                "Density" if case_id == CASES[0] else None,
            )
        )
    shared_density_top = max(density_maxima) * 1.16
    for ax in axes:
        ax.set_ylim(0, shared_density_top)
    draw_sigma_alpha_bar(bar_ax)
    fig.tight_layout(pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"local_distortion_sigma_bitacc_density_overlay{suffix}", bbox_inches="tight")
    plt.close(fig)

    for case_id in CASES:
        fig = plt.figure(figsize=(3.75, 2.45))
        grid_spec = fig.add_gridspec(1, 2, width_ratios=[1, 0.085], wspace=0.22)
        ax = fig.add_subplot(grid_spec[0, 0])
        bar_ax = fig.add_subplot(grid_spec[0, 1])
        draw_density_overlay_panel(ax, rows, case_id, "Density")
        draw_sigma_alpha_bar(bar_ax)
        fig.tight_layout(pad=0.55)
        for suffix in [".png", ".pdf", ".svg"]:
            fig.savefig(
                RESULT_DIR / f"local_distortion_sigma_bitacc_density_overlay_{case_id}{suffix}",
                bbox_inches="tight",
            )
        plt.close(fig)


if __name__ == "__main__":
    main()
