# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


RESULT_DIR = Path(__file__).resolve().parent
INPUT_CSV = RESULT_DIR / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv"

CASES = ["6MRR", "8HFE", "8VC8"]
SIGMAS = [0.0010, 0.0012, 0.0014, 0.0016, 0.0018, 0.0020]
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
                    "ber": float(row["ber_percent"]),
                    "exact": row["exact_recovery"].lower() == "true",
                    "clean_to_noisy_ca_rmsd": float(row["clean_to_noisy_ca_rmsd_angstrom"]),
                    "original_to_noisy_ca_rmsd": float(row["original_to_noisy_ca_rmsd_angstrom"]),
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


def draw_sigma_alpha_bar(bar_ax) -> None:
    alphas = np.linspace(0.22, 0.78, len(SIGMAS))
    y_positions = np.arange(len(SIGMAS))
    bars = bar_ax.barh(
        y_positions,
        np.ones(len(SIGMAS)),
        height=0.82,
        color="#4A4A4A",
        edgecolor="none",
    )
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)
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
        ax.fill_between(grid, 0, density, color=color, alpha=alpha * 0.30, linewidth=0, zorder=1)
        ax.plot(grid, density, color=color, alpha=alpha, linewidth=1.35, zorder=2)

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


def wilson_interval(successes: int, total: int, z_value: float = 1.96) -> tuple[float, float]:
    if total == 0:
        return 0.0, 0.0
    proportion = successes / total
    denominator = 1.0 + z_value**2 / total
    centre = proportion + z_value**2 / (2.0 * total)
    margin = z_value * np.sqrt((proportion * (1.0 - proportion) + z_value**2 / (4.0 * total)) / total)
    lower = (centre - margin) / denominator
    upper = (centre + margin) / denominator
    return lower * 100.0, upper * 100.0


def draw_summary_curves(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.85, 2.55), sharex=True)
    metrics = [
        ("bitacc", "Bit accuracy (%)", (0, 105)),
        ("exact", "Exact recovery (%)", (0, 105)),
        ("clean_to_noisy_ca_rmsd", "C-alpha RMSD (A)", (0, None)),
    ]
    for ax, (metric, ylabel, ylim) in zip(axes, metrics):
        for case_id in CASES:
            centers = []
            lower_errors = []
            upper_errors = []
            for sigma in SIGMAS:
                group = [row for row in rows if row["case_id"] == case_id and np.isclose(row["sigma"], sigma)]
                if metric == "exact":
                    successes = sum(bool(row["exact"]) for row in group)
                    center = successes / len(group) * 100.0
                    lower, upper = wilson_interval(successes, len(group))
                elif metric == "bitacc":
                    values = np.array([row[metric] for row in group], dtype=float)
                    lower, center, upper = np.percentile(values, [25, 50, 75])
                else:
                    values = np.array([row[metric] for row in group], dtype=float)
                    center = float(values.mean())
                    spread = float(values.std(ddof=1)) if len(values) > 1 else 0.0
                    lower = center - spread
                    upper = center + spread
                centers.append(float(center))
                lower_errors.append(float(center - lower))
                upper_errors.append(float(upper - center))
            color = COLORS[case_id]
            ax.errorbar(
                SIGMAS,
                centers,
                yerr=np.vstack([lower_errors, upper_errors]),
                color=color,
                marker="o",
                markersize=3.6,
                markerfacecolor="white",
                markeredgewidth=0.9,
                linewidth=1.25,
                capsize=1.8,
                label=case_id,
            )
        ax.set_xlabel("Gaussian noise sigma (10^-3 A)", fontsize=7.8, labelpad=3.5)
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=3.5)
        ax.set_ylim(*ylim) if ylim[1] is not None else ax.set_ylim(bottom=0)
        ax.set_xlim(SIGMAS[0] - 0.00008, SIGMAS[-1] + 0.00008)
        ax.set_xticks(SIGMAS)
        ax.set_xticklabels([f"{sigma * 1000:.1f}" for sigma in SIGMAS], rotation=0)
        if metric in {"bitacc", "exact"}:
            ax.axhline(100, color="#A0A0A0", linewidth=0.75, linestyle=":", zorder=0)
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=1.4)
    fig.tight_layout(w_pad=1.1, pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"gaussian_noise_all_atoms_sigma_0.001_to_0.002_summary_curves{suffix}", bbox_inches="tight")
    plt.close(fig)


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

    fig = plt.figure(figsize=(7.65, 2.35))
    grid_spec = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.075], wspace=0.34)
    axes = [fig.add_subplot(grid_spec[0, index]) for index in range(3)]
    axes[1].sharey(axes[0])
    axes[2].sharey(axes[0])
    bar_ax = fig.add_subplot(grid_spec[0, 3])
    density_maxima = []
    for ax, case_id in zip(axes, CASES):
        density_maxima.append(
            draw_density_overlay_panel(ax, rows, case_id, "Density" if case_id == CASES[0] else None)
        )
    shared_density_top = max(density_maxima) * 1.16
    for ax in axes:
        ax.set_ylim(0, shared_density_top)
    draw_sigma_alpha_bar(bar_ax)
    fig.tight_layout(pad=0.55)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"gaussian_noise_all_atoms_sigma_0.001_to_0.002_bitacc_density_overlay{suffix}", bbox_inches="tight")
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
                RESULT_DIR / f"gaussian_noise_all_atoms_sigma_0.001_to_0.002_bitacc_density_overlay_{case_id}{suffix}",
                bbox_inches="tight",
            )
        plt.close(fig)

    draw_summary_curves(rows)


if __name__ == "__main__":
    main()
