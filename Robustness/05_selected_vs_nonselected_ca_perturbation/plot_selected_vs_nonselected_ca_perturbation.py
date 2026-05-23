# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


RESULT_DIR = Path(__file__).resolve().parent
INPUT_CSV = RESULT_DIR / "selected_vs_nonselected_ca_perturbation_results.csv"

CASES = ["6MRR", "8HFE", "8VC8"]
GROUPS = ["nonselected_CA", "selected_CA"]
GROUP_LABELS = {
    "nonselected_CA": "Non-selected",
    "selected_CA": "Selected",
}
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
                    "target_group": row["target_group"],
                    "bitacc": float(row["bitacc_percent"]),
                    "exact": row["exact_recovery"].lower() == "true",
                    "sigma": float(row["sigma_angstrom"]),
                    "target_count": int(row["target_ca_count"]),
                }
            )
    return rows


def style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="major", direction="out", length=4.0, width=1.0, labelsize=7.4, pad=2.5)
    ax.tick_params(axis="both", which="minor", direction="out", length=2.4, width=0.8)


def gaussian_kde_manual(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    if len(values) == 0:
        return np.zeros_like(grid)
    diff = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * diff**2).sum(axis=1)
    density /= len(values) * bandwidth * np.sqrt(2 * np.pi)
    return density


def permutation_p_value(values_a: np.ndarray, values_b: np.ndarray, n_perm: int = 20000, seed: int = 9183) -> float:
    observed = abs(float(np.mean(values_a) - np.mean(values_b)))
    combined = np.concatenate([values_a, values_b])
    n_a = len(values_a)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        permuted = rng.permutation(combined)
        diff = abs(float(np.mean(permuted[:n_a]) - np.mean(permuted[n_a:])))
        if diff >= observed - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)


def significance_label(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "n.s."


def add_significance(ax, x1: float, x2: float, y: float, height: float, label: str) -> None:
    y_top = y + height
    ax.plot([x1, x1], [y, y_top], color="black", linewidth=0.9, linestyle="--", dashes=(3.2, 2.0), clip_on=False)
    ax.plot([x2, x2], [y_top, y], color="black", linewidth=0.9, linestyle="--", dashes=(3.2, 2.0), clip_on=False)
    ax.plot([x1, x2], [y_top, y_top], color="black", linewidth=0.9, linestyle="-", clip_on=False)
    if label:
        ax.text((x1 + x2) / 2, y + height * 1.08, label, ha="center", va="bottom", fontsize=8.4)


def draw_density_overlay_panel(
    ax,
    rows: list[dict[str, object]],
    case_id: str,
    p_value: float,
    ylabel: str | None = None,
    show_layer_labels: bool = False,
) -> float:
    color = COLORS[case_id]
    grid = np.linspace(30, 100, 760)
    layer_y = {
        "selected_CA": 0.0,
        "nonselected_CA": 1.05,
    }
    layer_labels = {
        "selected_CA": "Selected",
        "nonselected_CA": "Non-selected",
    }
    ridge_height = 0.68
    means = {}
    densities = {}
    max_density = 0.0

    for group in ["selected_CA", "nonselected_CA"]:
        values = np.asarray(
            [row["bitacc"] for row in rows if row["case_id"] == case_id and row["target_group"] == group],
            dtype=float,
        )
        means[group] = float(np.mean(values))
        bandwidth = 2.0
        density = gaussian_kde_manual(values, grid, bandwidth=bandwidth)
        densities[group] = density
        max_density = max(max_density, float(density.max()))

    top = layer_y["nonselected_CA"] + ridge_height * 1.42
    bracket_y = layer_y["nonselected_CA"] + ridge_height * 1.055
    bracket_height = ridge_height * 0.070

    for group in ["selected_CA", "nonselected_CA"]:
        baseline = layer_y[group]
        density = densities[group] / max_density * ridge_height if max_density > 0 else densities[group]
        mean_line_top = bracket_y if group == "selected_CA" else bracket_y + bracket_height
        ax.fill_between(
            grid,
            baseline,
            baseline + density,
            color=color,
            alpha=0.30,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            grid,
            baseline + density,
            color=color,
            alpha=1.0,
            linewidth=1.35,
            zorder=2,
        )
        ax.hlines(baseline, 30, 100, color="#202020", linewidth=0.75, zorder=3)
        ax.vlines(
            means[group],
            baseline,
            mean_line_top,
            color="#202020",
            linewidth=1.45,
            linestyle="--",
            zorder=4,
        )

    add_significance(
        ax,
        means["selected_CA"],
        means["nonselected_CA"],
        bracket_y,
        bracket_height,
        "",
    )
    ax.text(
        (means["selected_CA"] + means["nonselected_CA"]) / 2,
        bracket_y + bracket_height + 0.006,
        significance_label(p_value),
        ha="center",
        va="bottom",
        fontsize=8.8,
        color="black",
    )

    ax.set_title(case_id, fontsize=8.8, fontweight="bold", pad=3)
    ax.set_xlim(30, 100)
    ax.set_ylim(-0.06, top)
    ax.set_xlabel("Bit accuracy (%)", fontsize=7.8, labelpad=3.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7.8, labelpad=18)
    ax.set_yticks([layer_y["selected_CA"] + ridge_height * 0.42, layer_y["nonselected_CA"] + ridge_height * 0.42])
    if show_layer_labels:
        ax.set_yticklabels(["Selected", "Non-selected"], fontsize=7.4, rotation=90, va="center")
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    style_axes(ax)
    return top


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

    p_values = {}
    for case_id in CASES:
        nonselected = np.asarray(
            [row["bitacc"] for row in rows if row["case_id"] == case_id and row["target_group"] == "nonselected_CA"],
            dtype=float,
        )
        selected = np.asarray(
            [row["bitacc"] for row in rows if row["case_id"] == case_id and row["target_group"] == "selected_CA"],
            dtype=float,
        )
        p_values[case_id] = permutation_p_value(nonselected, selected, seed=9183 + CASES.index(case_id))

    with (RESULT_DIR / "selected_vs_nonselected_ca_perturbation_statistics.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "test",
                "p_value",
                "significance",
                "nonselected_mean_bitacc",
                "selected_mean_bitacc",
                "nonselected_median_bitacc",
                "selected_median_bitacc",
            ],
        )
        writer.writeheader()
        for case_id in CASES:
            nonselected = np.asarray(
                [
                    row["bitacc"]
                    for row in rows
                    if row["case_id"] == case_id and row["target_group"] == "nonselected_CA"
                ],
                dtype=float,
            )
            selected = np.asarray(
                [row["bitacc"] for row in rows if row["case_id"] == case_id and row["target_group"] == "selected_CA"],
                dtype=float,
            )
            writer.writerow(
                {
                    "case_id": case_id,
                    "test": "two-sided permutation test on mean BitAcc difference",
                    "p_value": p_values[case_id],
                    "significance": significance_label(p_values[case_id]),
                    "nonselected_mean_bitacc": float(np.mean(nonselected)),
                    "selected_mean_bitacc": float(np.mean(selected)),
                    "nonselected_median_bitacc": float(np.median(nonselected)),
                    "selected_median_bitacc": float(np.median(selected)),
                }
            )

    fig, axes = plt.subplots(1, 3, figsize=(7.65, 2.35), sharey=False)
    density_tops = []
    for ax, case_id in zip(axes, CASES):
        density_tops.append(
            draw_density_overlay_panel(
                ax,
                rows,
                case_id,
                p_values[case_id],
                "Density" if case_id == CASES[0] else None,
                show_layer_labels=case_id == CASES[0],
            )
        )
    fig.tight_layout(w_pad=1.25, pad=0.55)

    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(RESULT_DIR / f"selected_vs_nonselected_ca_bitacc_density_overlay{suffix}", bbox_inches="tight")
        fig.savefig(RESULT_DIR / f"selected_vs_nonselected_ca_bitacc_boxplot{suffix}", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
