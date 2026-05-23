# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


ATTACK_ROOT = Path(__file__).resolve().parents[1]
CASES_ROOT = ATTACK_ROOT / "3cases"
OUT_DIR = Path(__file__).resolve().parent

CASES = ["6MRR", "8HFE", "8VC8"]
COLORS = {
    "6MRR": "#0072B2",
    "8HFE": "#D55E00",
    "8VC8": "#009E73",
}
SIGMAS_LOCAL = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", direction="out", length=3.2, width=0.8, labelsize=6.0, pad=2.0)
    ax.tick_params(axis="both", which="minor", direction="out", length=1.9, width=0.7)


def panel_label(ax, label: str) -> None:
    ax.text(-0.15, 1.07, label, transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="top", ha="left")


def gaussian_kde_manual(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    if len(values) == 0:
        return np.zeros_like(grid)
    diff = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * diff**2).sum(axis=1)
    density /= len(values) * bandwidth * np.sqrt(2 * np.pi)
    return density


def draw_translation(ax) -> None:
    rows = read_csv(CASES_ROOT / "Global_translation" / "global_translation_results.csv")
    for case_id in CASES:
        group = [row for row in rows if row["case_id"] == case_id]
        x = np.asarray([float(row["translation_x_angstrom"]) for row in group], dtype=float)
        y = np.asarray([float(row["bitacc_percent"]) for row in group], dtype=float)
        order = np.argsort(x)
        ax.plot(
            x[order],
            y[order],
            color=COLORS[case_id],
            marker="o",
            markersize=3.5,
            markerfacecolor="white",
            markeredgewidth=0.9,
            linewidth=1.15,
            label=case_id,
        )
    ax.set_xscale("log")
    ax.set_ylim(96, 101.2)
    ax.set_xlabel("Translation (A)", fontsize=6.3, labelpad=2.4)
    ax.set_ylabel("Bit accuracy (%)", fontsize=6.3, labelpad=2.4)
    ax.axhline(100, color="#202020", linewidth=0.8, linestyle="--", dashes=(3.0, 2.0), zorder=0)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.text(0.06, 0.13, "100% recovery", transform=ax.transAxes, fontsize=5.9)
    ax.legend(frameon=False, fontsize=5.8, loc="lower right", handlelength=1.15)
    ax.set_title("Rigid-body translation", fontsize=7.0, pad=3)
    style_axes(ax)
    panel_label(ax, "A")


def draw_rotation(ax) -> None:
    rows = read_csv(CASES_ROOT / "Global_rotation" / "global_rotation_2angle_results.csv")
    summary = []
    for case_id in CASES:
        group = [row for row in rows if row["case_id"] == case_id]
        exact_fraction = sum(row["exact_recovery"].lower() == "true" for row in group) / len(group) * 100.0
        bitacc_min = min(float(row["bitacc_percent"]) for row in group)
        summary.append((case_id, exact_fraction, bitacc_min, len(group)))
    x = np.arange(len(CASES))
    ax.bar(
        x,
        [item[1] for item in summary],
        color=[COLORS[item[0]] for item in summary],
        alpha=0.55,
        edgecolor=[COLORS[item[0]] for item in summary],
        linewidth=1.0,
        width=0.58,
    )
    ax.scatter(x, [item[2] for item in summary], s=18, color="white", edgecolor="#202020", linewidth=0.8, zorder=3)
    for xpos, (_, _, _, n) in zip(x, summary):
        ax.text(xpos, 97.0, f"n={n}", ha="center", va="bottom", fontsize=5.6)
    ax.set_ylim(94, 102)
    ax.set_xticks(x)
    ax.set_xticklabels(CASES)
    ax.set_ylabel("Exact recovery (%)", fontsize=6.3, labelpad=2.4)
    ax.axhline(100, color="#202020", linewidth=0.8, linestyle="--", dashes=(3.0, 2.0), zorder=0)
    ax.set_title("Two-axis rigid-body rotation", fontsize=7.0, pad=3)
    style_axes(ax)
    panel_label(ax, "B")


def draw_local_density_panel(ax, rows: list[dict[str, str]], case_id: str, ylabel: str | None = None) -> float:
    color = COLORS[case_id]
    grid = np.linspace(30, 100, 760)
    alphas = np.linspace(0.22, 0.78, len(SIGMAS_LOCAL))
    max_density = 0.0
    for sigma, alpha in zip(SIGMAS_LOCAL, alphas):
        values = np.asarray(
            [
                float(row["bitacc_percent"])
                for row in rows
                if row["case_id"] == case_id and np.isclose(float(row["sigma_angstrom"]), sigma)
            ],
            dtype=float,
        )
        density = gaussian_kde_manual(values, grid, bandwidth=2.0)
        max_density = max(max_density, float(density.max()))
        ax.fill_between(grid, 0, density, color=color, alpha=alpha * 0.30, linewidth=0, zorder=1)
        ax.plot(grid, density, color=color, alpha=alpha, linewidth=1.05, zorder=2)
    ax.axvline(100, color="#202020", linewidth=1.1, linestyle="--", dashes=(3.0, 2.0), zorder=4)
    ax.set_title(case_id, fontsize=6.6, fontweight="bold", pad=2)
    ax.set_xlim(30, 100)
    ax.set_xlabel("Bit accuracy (%)", fontsize=6.2, labelpad=2.2)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=6.2, labelpad=2.4)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    style_axes(ax)
    return max_density


def draw_selected_ridge_panel(ax, rows: list[dict[str, str]], case_id: str, show_labels: bool = False) -> None:
    color = COLORS[case_id]
    grid = np.linspace(30, 100, 760)
    layer_y = {"selected_CA": 0.0, "nonselected_CA": 1.02}
    ridge_height = 0.62
    means = {}
    densities = {}
    max_density = 0.0
    for group in ["selected_CA", "nonselected_CA"]:
        values = np.asarray(
            [float(row["bitacc_percent"]) for row in rows if row["case_id"] == case_id and row["target_group"] == group],
            dtype=float,
        )
        means[group] = float(np.mean(values))
        density = gaussian_kde_manual(values, grid, bandwidth=2.0)
        densities[group] = density
        max_density = max(max_density, float(density.max()))

    bracket_y = layer_y["nonselected_CA"] + ridge_height * 1.04
    bracket_h = ridge_height * 0.08
    for group in ["selected_CA", "nonselected_CA"]:
        baseline = layer_y[group]
        density = densities[group] / max_density * ridge_height if max_density > 0 else densities[group]
        line_top = bracket_y if group == "selected_CA" else bracket_y + bracket_h
        ax.fill_between(grid, baseline, baseline + density, color=color, alpha=0.30, linewidth=0)
        ax.plot(grid, baseline + density, color=color, linewidth=1.15)
        ax.hlines(baseline, 30, 100, color="#202020", linewidth=0.65)
        ax.vlines(means[group], baseline, line_top, color="#202020", linewidth=1.15, linestyle="--", zorder=4)
    x1, x2 = means["selected_CA"], means["nonselected_CA"]
    ax.plot([x1, x1], [bracket_y, bracket_y + bracket_h], color="black", linewidth=0.8, linestyle="--", dashes=(3, 2))
    ax.plot([x2, x2], [bracket_y + bracket_h, bracket_y], color="black", linewidth=0.8, linestyle="--", dashes=(3, 2))
    ax.plot([x1, x2], [bracket_y + bracket_h, bracket_y + bracket_h], color="black", linewidth=0.8)
    ax.text((x1 + x2) / 2, bracket_y + bracket_h + 0.005, "***", ha="center", va="bottom", fontsize=6.9)
    ax.set_xlim(30, 100)
    ax.set_ylim(-0.05, layer_y["nonselected_CA"] + ridge_height * 1.34)
    ax.set_xlabel("Bit accuracy (%)", fontsize=6.2, labelpad=2.2)
    ax.set_title(case_id, fontsize=6.6, fontweight="bold", pad=2)
    ax.set_yticks([layer_y["selected_CA"] + ridge_height * 0.42, layer_y["nonselected_CA"] + ridge_height * 0.42])
    if show_labels:
        ax.set_yticklabels(["Selected", "Non-selected"], fontsize=5.9, rotation=90, va="center")
        ax.set_ylabel("Density", fontsize=6.2, labelpad=15)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    style_axes(ax)


def draw_low_sigma(ax) -> None:
    rows = read_csv(
        CASES_ROOT
        / "Gaussian_coordinate_noise_all_atoms_low_sigma"
        / "sigma_0.001_to_0.002"
        / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv"
    )
    sigmas = sorted({float(row["sigma_angstrom"]) for row in rows})
    for case_id in CASES:
        exact_rates = []
        med_bitacc = []
        for sigma in sigmas:
            group = [row for row in rows if row["case_id"] == case_id and np.isclose(float(row["sigma_angstrom"]), sigma)]
            exact_rates.append(sum(row["exact_recovery"].lower() == "true" for row in group) / len(group) * 100.0)
            med_bitacc.append(float(np.median([float(row["bitacc_percent"]) for row in group])))
        x = np.asarray(sigmas) * 1000.0
        ax.plot(
            x,
            exact_rates,
            color=COLORS[case_id],
            marker="o",
            markersize=3.4,
            markerfacecolor="white",
            markeredgewidth=0.9,
            linewidth=1.1,
            label=case_id,
        )
        ax.plot(x, med_bitacc, color=COLORS[case_id], linewidth=0.75, linestyle=":", alpha=0.85)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Noise sigma (10$^{-3}$ A)", fontsize=6.3, labelpad=2.4)
    ax.set_ylabel("Exact recovery (%)", fontsize=6.3, labelpad=2.4)
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.axhline(100, color="#202020", linewidth=0.8, linestyle="--", dashes=(3.0, 2.0), zorder=0)
    ax.set_title("Low-amplitude all-atom noise", fontsize=7.0, pad=3)
    ax.legend(frameon=False, fontsize=5.7, loc="lower left", handlelength=1.15)
    style_axes(ax)
    panel_label(ax, "E")


def write_caption() -> None:
    caption = """Figure X | Robustness and tamper-boundary analysis of StrucTrace decoding.
a, Rigid-body translation of the watermarked structures along the x axis over 0.1-100 A did not affect payload recovery in the three representative proteins (6MRR, 8HFE and 8VC8), with 100% bit accuracy and exact recovery in all cases. b, Two-axis rigid-body rotations around the x and z axes followed by Kabsch alignment also preserved exact recovery across the full rotation grid (n = 169 angle combinations per protein). c, Local coordinate distortion was tested by randomly selecting a C-alpha-centred spherical region (radius, 2 A) and applying Gaussian coordinate noise with sigma values from 0.05 to 2 A over 50 repeats per condition. Density overlays show the resulting bit-accuracy distributions; line opacity increases with perturbation strength. d, Random perturbation of 12 watermark-bearing selected C-alpha atoms (sigma = 0.05 A, 100 repeats) substantially reduced bit accuracy, whereas perturbing the same number of non-selected C-alpha atoms preserved complete recovery. Dashed vertical lines denote group means; *** denotes P < 0.001 by a two-sided permutation test. e, Low-amplitude all-atom Gaussian coordinate noise defined a fine robustness threshold: exact recovery remained complete at 0.001 A but declined as sigma approached 0.002 A, despite high median bit accuracy (dotted lines). All panels report decoding of the 56-bit UTF-8 payload "npj SB" plus a null terminator."""
    (OUT_DIR / "main_robustness_figure_caption.txt").write_text(caption + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 6.2,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )

    fig = plt.figure(figsize=(7.25, 5.85))
    outer = fig.add_gridspec(3, 3, height_ratios=[0.92, 1.0, 1.0], hspace=0.70, wspace=0.46)

    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    ax_e = fig.add_subplot(outer[0, 2])
    draw_translation(ax_a)
    draw_rotation(ax_b)
    draw_low_sigma(ax_e)

    local_rows = read_csv(CASES_ROOT / "Local_coordinate_distortion" / "local_distortion_sigma_results.csv")
    selected_rows = read_csv(
        ATTACK_ROOT
        / "Selected_vs_nonselected_CA_perturbation"
        / "selected_vs_nonselected_ca_perturbation_results.csv"
    )

    local_grid = outer[1, :].subgridspec(1, 3, wspace=0.28)
    local_axes = [fig.add_subplot(local_grid[0, i]) for i in range(3)]
    maxima = []
    for ax, case_id in zip(local_axes, CASES):
        maxima.append(draw_local_density_panel(ax, local_rows, case_id, "Density" if case_id == CASES[0] else None))
    shared_top = max(maxima) * 1.16
    for ax in local_axes:
        ax.set_ylim(0, shared_top)
    panel_label(local_axes[0], "C")
    local_axes[1].text(
        0.5,
        1.12,
        "Local coordinate distortion",
        transform=local_axes[1].transAxes,
        ha="center",
        va="bottom",
        fontsize=7.0,
    )

    selected_grid = outer[2, :].subgridspec(1, 3, wspace=0.28)
    selected_axes = [fig.add_subplot(selected_grid[0, i]) for i in range(3)]
    for i, (ax, case_id) in enumerate(zip(selected_axes, CASES)):
        draw_selected_ridge_panel(ax, selected_rows, case_id, show_labels=i == 0)
    panel_label(selected_axes[0], "D")
    selected_axes[1].text(
        0.5,
        1.12,
        "Selected versus non-selected C-alpha perturbation",
        transform=selected_axes[1].transAxes,
        ha="center",
        va="bottom",
        fontsize=7.0,
    )

    fig.savefig(OUT_DIR / "main_robustness_figure.png", bbox_inches="tight")
    fig.savefig(OUT_DIR / "main_robustness_figure.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "main_robustness_figure.svg", bbox_inches="tight")
    plt.close(fig)
    write_caption()


if __name__ == "__main__":
    main()
