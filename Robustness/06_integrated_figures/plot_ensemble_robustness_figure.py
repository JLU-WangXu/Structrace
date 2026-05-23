# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


ATTACK_ROOT = Path(__file__).resolve().parents[1]
CASES_ROOT = ATTACK_ROOT / "3cases"
SELECTED_ROOT = ATTACK_ROOT / "Selected_vs_nonselected_CA_perturbation"
OUT_DIR = Path(__file__).resolve().parent

CASES = ["6MRR", "8HFE", "8VC8"]
COLORS = {
    "6MRR": "#0072B2",
    "8HFE": "#D55E00",
    "8VC8": "#009E73",
}
MARKERS = {"6MRR": "o", "8HFE": "s", "8VC8": "^"}
LOCAL_SIGMAS = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
ALL_ATOM_SIGMAS = [0.0010, 0.0012, 0.0014, 0.0016, 0.0018, 0.0020]

LABEL_SIZE = 6.4
TICK_SIZE = 5.8
LEGEND_SIZE = 5.8
CASE_SIZE = 6.2
PANEL_SIZE = 9.0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def set_bold_ticks(ax) -> None:
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")
        tick.set_fontfamily("Arial")


def style_axes(ax, hide_top_right: bool = False) -> None:
    if hide_top_right:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.75)
        spine.set_alpha(0.58)
    ax.tick_params(axis="both", which="major", direction="out", length=3.0, width=0.75, labelsize=TICK_SIZE, pad=1.8)
    ax.tick_params(axis="both", which="minor", direction="out", length=1.8, width=0.65)
    set_bold_ticks(ax)


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.07,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        fontfamily="Arial",
        ha="left",
        va="bottom",
    )


def case_label(ax, case_id: str) -> None:
    ax.text(
        0.04,
        0.92,
        case_id,
        transform=ax.transAxes,
        fontsize=CASE_SIZE,
        fontweight="bold",
        fontfamily="Arial",
        ha="left",
        va="top",
    )


def gaussian_kde_manual(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    if len(values) == 0:
        return np.zeros_like(grid)
    diff = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * diff**2).sum(axis=1)
    density /= len(values) * bandwidth * np.sqrt(2.0 * np.pi)
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
            marker=MARKERS[case_id],
            markersize=3.2,
            markerfacecolor="white",
            markeredgewidth=0.8,
            linewidth=1.0,
            label=case_id,
        )
    ax.set_xscale("log")
    ax.set_xlim(0.08, 130)
    ax.set_ylim(94, 101.2)
    ax.set_xlabel("Translation (A)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.0)
    ax.set_ylabel("Bit accuracy (%)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.2)
    ax.axhline(100, color="#2B2B2B", linewidth=0.75, linestyle="--", dashes=(3, 2), zorder=0)
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    leg = ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower right", handlelength=1.25, borderpad=0.1)
    for text in leg.get_texts():
        text.set_fontweight("bold")
    style_axes(ax)
    panel_label(ax, "a")


def draw_rotation(ax_x, ax_z) -> None:
    rows = read_csv(CASES_ROOT / "Global_rotation" / "global_rotation_results.csv")
    for ax, axis_name in [(ax_x, "x"), (ax_z, "z")]:
        for case_id in CASES:
            group = [row for row in rows if row["case_id"] == case_id and row["rotation_axis"] == axis_name]
            x = np.asarray([float(row["angle_deg"]) for row in group], dtype=float)
            y = np.asarray([float(row["bitacc_percent"]) for row in group], dtype=float)
            order = np.argsort(x)
            ax.plot(
                x[order],
                y[order],
                color=COLORS[case_id],
                marker=MARKERS[case_id],
                markersize=3.0,
                markerfacecolor="white",
                markeredgewidth=0.8,
                linewidth=0.95,
                label=case_id,
            )
        ax.set_xlim(-12, 372)
        ax.set_ylim(94, 101.2)
        ax.axhline(100, color="#2B2B2B", linewidth=0.75, linestyle="--", dashes=(3, 2), zorder=0)
        ax.xaxis.set_major_locator(MultipleLocator(90))
        ax.xaxis.set_minor_locator(MultipleLocator(30))
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.text(
            0.04,
            0.12,
            f"{axis_name}-axis",
            transform=ax.transAxes,
            fontsize=LEGEND_SIZE,
            fontweight="bold",
            fontfamily="Arial",
            ha="left",
            va="bottom",
        )
        style_axes(ax)
    ax_x.tick_params(axis="x", labelbottom=False)
    ax_x.set_xlabel("")
    ax_z.set_xlabel("Rotation angle (deg)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.0)
    ax_x.set_ylabel("")
    ax_z.set_ylabel("")
    ax_x.text(
        -0.17,
        -0.12,
        "Bit accuracy (%)",
        transform=ax_x.transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        fontfamily="Arial",
    )
    leg = ax_z.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower right", handlelength=1.2, borderpad=0.1)
    for text in leg.get_texts():
        text.set_fontweight("bold")
    panel_label(ax_x, "b")


def draw_alpha_bar(ax, values: list[float], label: str) -> None:
    alphas = np.linspace(0.22, 0.78, len(values))
    y_positions = np.arange(len(values))
    bars = ax.barh(y_positions, np.ones(len(values)), height=0.82, color="#4A4A4A", edgecolor="none")
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.55, len(values) - 0.45)
    ax.set_xticks([])
    ax.set_yticks(y_positions)
    labels = [f"{value:g}" for value in values]
    if max(values) < 0.01:
        labels = [f"{value * 1000:.1f}" for value in values]
    ax.set_yticklabels(labels, fontsize=TICK_SIZE, fontweight="bold")
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", direction="out", length=2.0, width=0.65, pad=1.8)
    ax.text(
        0.5,
        1.02,
        label,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=LEGEND_SIZE,
        fontweight="bold",
        fontfamily="Arial",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    set_bold_ticks(ax)


def draw_alpha_bar_horizontal(ax, values: list[float], label: str) -> None:
    alphas = np.linspace(0.22, 0.78, len(values))
    x_positions = np.arange(len(values))
    bars = ax.bar(x_positions, np.ones(len(values)), width=0.82, color="#4A4A4A", edgecolor="none")
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)
    labels = [f"{value:g}" for value in values]
    if max(values) < 0.01:
        labels = [f"{value * 1000:.1f}" for value in values]
    ax.set_xlim(-0.55, len(values) - 0.45)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for xpos, text, alpha in zip(x_positions, labels, alphas):
        text_color = "black" if alpha < 0.58 else "white"
        ax.text(
            xpos,
            0.50,
            text,
            ha="center",
            va="center",
            fontsize=TICK_SIZE,
            fontweight="bold",
            fontfamily="Arial",
            color=text_color,
        )
    ax.text(
        0.5,
        1.25,
        label,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=LEGEND_SIZE,
        fontweight="bold",
        fontfamily="Arial",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    set_bold_ticks(ax)


def draw_colored_alpha_inset(ax, values: list[float], label: str, color: str) -> None:
    labels = [f"{value:g}" for value in values]
    if max(values) < 0.01:
        labels = [f"{value * 1000:.1f}" for value in values]

    x0, y0 = 0.51, 0.845
    total_w, h = 0.45, 0.092
    rgb = np.asarray(mpl.colors.to_rgb(color), dtype=float)
    gradient = np.ones((1, 192, 4), dtype=float)
    gradient[:, :, :3] = rgb
    gradient[:, :, 3] = np.linspace(0.22, 0.78, gradient.shape[1]) * 0.62
    ax.imshow(
        gradient,
        extent=(x0, x0 + total_w, y0, y0 + h),
        transform=ax.transAxes,
        aspect="auto",
        interpolation="bicubic",
        zorder=9,
    )
    ax.text(
        x0,
        y0 - 0.014,
        labels[0],
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=max(4.2, TICK_SIZE - 1.2),
        fontweight="bold",
        fontfamily="Arial",
        zorder=10,
    )
    ax.text(
        x0 + total_w,
        y0 - 0.014,
        labels[-1],
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=max(4.2, TICK_SIZE - 1.2),
        fontweight="bold",
        fontfamily="Arial",
        zorder=10,
    )
    ax.text(
        x0 + total_w / 2,
        y0 - 0.064,
        label,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=LEGEND_SIZE - 0.4,
        fontweight="bold",
        fontfamily="Arial",
        zorder=10,
    )


def draw_density_overlay(
    ax,
    rows: list[dict[str, str]],
    case_id: str,
    sigma_key: str,
    sigmas: list[float],
    ylabel: bool,
    alpha_label: str | None = None,
) -> float:
    color = COLORS[case_id]
    grid = np.linspace(30, 100, 760)
    alphas = np.linspace(0.22, 0.78, len(sigmas))
    max_density = 0.0
    for sigma, alpha in zip(sigmas, alphas):
        values = np.asarray(
            [
                float(row["bitacc_percent"])
                for row in rows
                if row["case_id"] == case_id and np.isclose(float(row[sigma_key]), sigma)
            ],
            dtype=float,
        )
        density = gaussian_kde_manual(values, grid, bandwidth=2.0)
        max_density = max(max_density, float(density.max()))
        ax.fill_between(grid, 0, density, color=color, alpha=alpha * 0.30, linewidth=0, zorder=1)
        ax.plot(grid, density, color=color, alpha=alpha, linewidth=1.0, zorder=2)
    ax.axvline(100, color="#2B2B2B", linewidth=0.9, linestyle="--", dashes=(3, 2), zorder=4)
    ax.set_xlim(30, 100)
    ax.set_xlabel("Bit accuracy (%)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.0)
    if ylabel:
        ax.set_ylabel("Density", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.2)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    case_label(ax, case_id)
    if alpha_label:
        draw_colored_alpha_inset(ax, sigmas, alpha_label, color)
    style_axes(ax)
    return max_density


def draw_selected_ridge(ax, rows: list[dict[str, str]], case_id: str, show_labels: bool) -> None:
    color = COLORS[case_id]
    grid = np.linspace(30, 100, 760)
    layer_y = {"selected_CA": 0.0, "nonselected_CA": 1.05}
    ridge_height = 0.68
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

    bracket_y = layer_y["nonselected_CA"] + ridge_height * 1.055
    bracket_h = ridge_height * 0.070
    for group in ["selected_CA", "nonselected_CA"]:
        baseline = layer_y[group]
        density = densities[group] / max_density * ridge_height if max_density > 0 else densities[group]
        mean_line_top = bracket_y if group == "selected_CA" else bracket_y + bracket_h
        ax.fill_between(grid, baseline, baseline + density, color=color, alpha=0.30, linewidth=0, zorder=1)
        ax.plot(grid, baseline + density, color=color, linewidth=1.05, zorder=2)
        ax.hlines(baseline, 30, 100, color="#2B2B2B", linewidth=0.65, zorder=3)
        ax.vlines(means[group], baseline, mean_line_top, color="#2B2B2B", linewidth=1.0, linestyle="--", zorder=4)
    x1, x2 = means["selected_CA"], means["nonselected_CA"]
    ax.plot([x1, x1], [bracket_y, bracket_y + bracket_h], color="black", linewidth=0.75, linestyle="--", dashes=(3, 2))
    ax.plot([x2, x2], [bracket_y + bracket_h, bracket_y], color="black", linewidth=0.75, linestyle="--", dashes=(3, 2))
    bracket_top = bracket_y + bracket_h
    ylim_top = layer_y["nonselected_CA"] + ridge_height * 1.42
    ax.plot([x1, x2], [bracket_top, bracket_top], color="black", linewidth=0.75)
    ax.text(
        (x1 + x2) / 2,
        bracket_top + (ylim_top - bracket_top) * 0.36,
        "***",
        ha="center",
        va="center",
        fontsize=LEGEND_SIZE,
        fontweight="bold",
        fontfamily="Arial",
    )
    ax.set_xlim(30, 100)
    ax.set_ylim(-0.06, ylim_top)
    ax.set_xlabel("Bit accuracy (%)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.0)
    ax.set_yticks([layer_y["selected_CA"] + ridge_height * 0.50, layer_y["nonselected_CA"] + ridge_height * 0.50])
    if show_labels:
        ax.set_yticklabels(["Selected", "Non-selected"], fontsize=TICK_SIZE, fontweight="bold", rotation=90, va="center")
        ax.set_ylabel("Density", fontsize=LABEL_SIZE, fontweight="bold", labelpad=14)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    case_label(ax, case_id)
    style_axes(ax)


def draw_rounding(ax) -> None:
    rows = read_csv(CASES_ROOT / "Coordinate_decimal_rounding" / "coordinate_decimal_rounding_results.csv")
    precision_order = [0.001, 0.01, 0.1, 1.0]
    x = np.arange(len(precision_order))
    for case_id in CASES:
        y = []
        for precision in precision_order:
            match = [
                row
                for row in rows
                if row["case_id"] == case_id and np.isclose(float(row["coordinate_step_angstrom"]), precision)
            ][0]
            y.append(float(match["bitacc_percent"]))
        ax.plot(
            x,
            y,
            color=COLORS[case_id],
            marker=MARKERS[case_id],
            markersize=3.4,
            markerfacecolor="white",
            markeredgewidth=0.85,
            linewidth=1.05,
            label=case_id,
        )
    ax.axhline(100, color="#2B2B2B", linewidth=0.75, linestyle="--", dashes=(3, 2), zorder=0)
    ax.set_xlim(-0.15, len(precision_order) - 0.85)
    ax.set_ylim(45, 102.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["0.001 A", "0.01 A", "0.1 A", "1 A"])
    ax.set_xlabel("Coordinate precision", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.0)
    ax.set_ylabel("Bit accuracy (%)", fontsize=LABEL_SIZE, fontweight="bold", labelpad=2.2)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    leg = ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left", handlelength=1.25, borderpad=0.1)
    for text in leg.get_texts():
        text.set_fontweight("bold")
    style_axes(ax)
    panel_label(ax, "c")


def write_caption() -> None:
    caption = """Figure X | Multi-scale robustness and precision limits of StrucTrace watermark decoding.
a, Rigid-body translation of watermarked structures along the x axis over 0.1-100 A preserved complete payload recovery in all three representative proteins. b, Two-angle rigid-body rotation was assessed by independently rotating structures around the x and z axes over 0-360 degrees, followed by Kabsch alignment; bit accuracy and exact recovery remained 100% throughout. c, Coordinate-rounding analysis shows that exact recovery is retained at the standard PDB coordinate precision of 0.001 A, whereas coarser, non-standard precisions progressively compromise bit recovery. d, Low-amplitude all-atom Gaussian coordinate noise was evaluated over sigma = 0.0010-0.0020 A with 50 repeats per condition; density overlays show the resulting bit-accuracy distributions, with darker curves denoting larger sigma. e, Local coordinate distortion was tested by applying Gaussian noise to a randomly selected C-alpha-centred spherical region (radius, 2 A; sigma = 0.05-2 A; 50 repeats per condition). f, White-box C-alpha perturbation defined the tamper boundary: perturbing 12 non-selected C-alpha atoms preserved recovery, whereas perturbing 12 watermark-bearing selected C-alpha atoms substantially reduced bit accuracy (sigma = 0.05 A; 100 repeats; dashed vertical lines denote group means; ***P < 0.001 by two-sided permutation test). All panels decode the 56-bit UTF-8 payload "npj SB" plus a null terminator."""
    (OUT_DIR / "ensemble_robustness_figure_caption.txt").write_text(caption + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 600,
        }
    )

    fig = plt.figure(figsize=(7.25, 7.05))
    outer = fig.add_gridspec(4, 1, height_ratios=[1.28, 1.0, 1.0, 1.08], hspace=0.43)

    top_grid = outer[0].subgridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.28)
    ax_a = fig.add_subplot(top_grid[0, 0])
    rot_grid = top_grid[0, 1].subgridspec(2, 1, hspace=0.10)
    ax_bx = fig.add_subplot(rot_grid[0, 0])
    ax_bz = fig.add_subplot(rot_grid[1, 0], sharex=ax_bx, sharey=ax_bx)
    ax_c = fig.add_subplot(top_grid[0, 2])
    draw_translation(ax_a)
    draw_rotation(ax_bx, ax_bz)
    draw_rounding(ax_c)

    all_atom_rows = read_csv(
        CASES_ROOT
        / "Gaussian_coordinate_noise_all_atoms_low_sigma"
        / "sigma_0.001_to_0.002"
        / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv"
    )
    d_grid = outer[1].subgridspec(1, 3, wspace=0.28)
    d_axes = [fig.add_subplot(d_grid[0, i]) for i in range(3)]
    d_max = [
        draw_density_overlay(
            ax,
            all_atom_rows,
            case_id,
            "sigma_angstrom",
            ALL_ATOM_SIGMAS,
            i == 0,
            r"$\sigma$ ($10^{-3}$ A)",
        )
        for i, (ax, case_id) in enumerate(zip(d_axes, CASES))
    ]
    for ax in d_axes:
        ax.set_ylim(0, max(d_max) * 1.34)
    panel_label(d_axes[0], "d")

    local_rows = read_csv(CASES_ROOT / "Local_coordinate_distortion" / "local_distortion_sigma_results.csv")
    e_grid = outer[2].subgridspec(1, 3, wspace=0.28)
    e_axes = [fig.add_subplot(e_grid[0, i]) for i in range(3)]
    e_max = [
        draw_density_overlay(
            ax,
            local_rows,
            case_id,
            "sigma_angstrom",
            LOCAL_SIGMAS,
            i == 0,
            r"$\sigma$ (A)",
        )
        for i, (ax, case_id) in enumerate(zip(e_axes, CASES))
    ]
    for ax in e_axes:
        ax.set_ylim(0, max(e_max) * 1.34)
    panel_label(e_axes[0], "e")

    selected_rows = read_csv(SELECTED_ROOT / "selected_vs_nonselected_ca_perturbation_results.csv")
    f_grid = outer[3].subgridspec(1, 3, wspace=0.28)
    f_axes = [fig.add_subplot(f_grid[0, i]) for i in range(3)]
    for i, (ax, case_id) in enumerate(zip(f_axes, CASES)):
        draw_selected_ridge(ax, selected_rows, case_id, show_labels=i == 0)
    panel_label(f_axes[0], "f")

    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(OUT_DIR / f"ensemble_robustness_figure{suffix}", bbox_inches="tight")
    plt.close(fig)
    write_caption()


if __name__ == "__main__":
    main()
