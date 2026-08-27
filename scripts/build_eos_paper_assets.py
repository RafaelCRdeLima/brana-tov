#!/usr/bin/env python3
"""Build multi-EoS figures and a LaTeX table from a campaign summary."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


DISPLAY_NAMES = {
    "BSK24": "BSk24",
    "DD2": "DD2",
    "NL3": "NL3",
    "QHC19A": "QHC19-A",
    "QHC21D": "QHC21-D",
    "SLy": "SLy",
    "SPGM4": "SPG(M4)",
}

ORDER = ["SLy", "BSK24", "DD2", "NL3", "QHC19A", "QHC21D", "SPGM4"]


def load_summary(path: Path) -> list[dict[str, float | str]]:
    with path.open(encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    converted: list[dict[str, float | str]] = []
    for row in rows:
        converted.append({key: value if key == "eos" else float(value) for key, value in row.items()})
    lookup = {str(row["eos"]): row for row in converted}
    return [lookup[name] for name in ORDER if name in lookup]


def load_sequence(path: Path) -> tuple[list[float], list[float]]:
    masses: list[float] = []
    radii: list[float] = []
    with path.open(encoding="utf-8") as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            radius = float(fields[3])
            mass = float(fields[5])
            if 7.0 <= radius <= 20.0 and mass >= 0.5:
                radii.append(radius)
                masses.append(mass)
    return radii, masses


def plot_gr_sequences(rows: list[dict[str, float | str]], sequence_root: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    cmap = plt.get_cmap("tab10")
    for index, row in enumerate(rows):
        eos = str(row["eos"])
        radii, masses = load_sequence(sequence_root / eos / "sequence_gr.dat")
        ax.plot(radii, masses, lw=2.0, color=cmap(index), label=DISPLAY_NAMES[eos])
        ax.scatter([float(row["R14_GR"])], [1.4], s=18, color=cmap(index), zorder=3)
    ax.axhline(2.0, color="0.55", lw=1.0, ls=":", label=r"$2\,M_\odot$")
    ax.set_xlim(9.0, 16.2)
    ax.set_ylim(0.7, 3.0)
    ax.set_xlabel("Radius [km]")
    ax.set_ylabel(r"Mass [$M_\odot$]")
    ax.grid(alpha=0.2)
    ax.legend(ncol=2, fontsize=9, frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_posterior_summary(rows: list[dict[str, float | str]], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.8), sharey=True)
    y = list(range(len(rows)))
    labels = [DISPLAY_NAMES[str(row["eos"])] for row in rows]
    quantities = [
        ("alpha_U", r"$\alpha_{\mathcal{U}}$", 0.0),
        ("Mmax", r"$M_{\max}$ [$M_\odot$]", None),
        ("R14", r"$R_{1.4}$ [km]", None),
    ]
    for ax, (key, xlabel, reference) in zip(axes, quantities):
        medians = [float(row[f"{key}_q50"]) for row in rows]
        lower = [median - float(row[f"{key}_q16"]) for median, row in zip(medians, rows)]
        upper = [float(row[f"{key}_q84"]) - median for median, row in zip(medians, rows)]
        ax.errorbar(medians, y, xerr=[lower, upper], fmt="o", capsize=3, color="#1f5a94")
        if reference is not None:
            ax.axvline(reference, color="0.35", ls="--", lw=1.2)
        if key == "Mmax":
            ax.scatter([float(row["Mmax_GR"]) for row in rows], y, marker="x", color="#b23a48", label="GR")
            ax.legend(frameon=False, fontsize=9, loc="upper left")
        if key == "R14":
            ax.scatter([float(row["R14_GR"]) for row in rows], y, marker="x", color="#b23a48")
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.2)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def interval(row: dict[str, float | str], key: str, digits: int = 2) -> str:
    median = float(row[f"{key}_q50"])
    lower = median - float(row[f"{key}_q16"])
    upper = float(row[f"{key}_q84"]) - median
    return f"${median:.{digits}f}^{{+{upper:.{digits}f}}}_{{-{lower:.{digits}f}}}$"


def write_latex_table(rows: list[dict[str, float | str]], output: Path, campaign_label: str) -> None:
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\caption{Multi-EoS posterior summary, ordered by the general-relativistic canonical radius for the " + campaign_label + r" campaign. Intervals are 68\% credible intervals; $M_{\max}^{\rm GR}$ and $R_{1.4}^{\rm GR}$ are the corresponding GR baselines.}",
        r"\label{tab:multi_eos}",
        r"\begin{tabular}{lccccccc}",
        r"\hline",
        r"EoS & $M_{\max}^{\rm GR}$ [$M_\odot$] & $R_{1.4}^{\rm GR}$ [km] & $\alpha_{\mathcal U}$ & $w_{\mathcal U}$ & $P(\alpha_{\mathcal U}<0)$ & $M_{\max}$ [$M_\odot$] & $R_{1.4}$ [km] \\",
        r"\hline",
    ]
    for row in rows:
        eos = DISPLAY_NAMES[str(row["eos"])]
        lines.append(
            f"{eos} & {float(row['Mmax_GR']):.3f} & {float(row['R14_GR']):.2f} & "
            f"{interval(row, 'alpha_U', 2)} & {interval(row, 'w_U', 2)} & "
            f"{100.0 * float(row['P_alpha_lt_0']):.1f}\\% & "
            f"{interval(row, 'Mmax', 2)} & {interval(row, 'R14', 2)} \\\\"
        )
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table*}", ""])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")


def write_diagnostics_table(path: Path, output: Path) -> None:
    with path.open(encoding="utf-8") as source:
        lookup = {row["eos"]: row for row in csv.DictReader(source)}
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\caption{Production-chain diagnostics after the 25\% burn-in.}",
        r"\label{tab:multi_eos_diagnostics}",
        r"\begin{tabular}{lccccc}",
        r"\hline",
        r"EoS & steps & acceptance & $\hat R_{\log_{10}\lambda}$ & $\hat R_{\alpha_{\mathcal U}}$ & $\hat R_{w_{\mathcal U}}$ \\",
        r"\hline",
    ]
    for eos in ORDER:
        if eos not in lookup:
            continue
        row = lookup[eos]
        lines.append(
            f"{DISPLAY_NAMES[eos]} & {int(float(row['n_steps']))} & {float(row['acceptance_rate_global']):.3f} & "
            f"{float(row['rhat_log10_lambda']):.4f} & {float(row['rhat_alpha_U']):.4f} & "
            f"{float(row['rhat_w_U']):.4f} \\\\"
        )
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table}", ""])
    output.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("output/campaign/smoke/summary.csv"))
    parser.add_argument("--sequences", type=Path, default=Path("output/campaign/sequences"))
    parser.add_argument("--output-dir", type=Path, default=Path("paper2"))
    parser.add_argument("--campaign-label", default="exploratory")
    parser.add_argument("--diagnostics", type=Path)
    args = parser.parse_args()
    rows = load_summary(args.summary)
    # Order by the GR canonical radius. The inferred Weyl coupling is monotonic in
    # this quantity across the whole set, so the ordering makes the paper's central
    # result visible in both the table and the summary figure instead of stating it.
    rows = sorted(rows, key=lambda r: float(r["R14_GR"]))
    plot_gr_sequences(rows, args.sequences, args.output_dir / "multi_eos_gr_sequences.png")
    plot_posterior_summary(rows, args.output_dir / "multi_eos_posterior_summary.png")
    write_latex_table(rows, args.output_dir / "multi_eos_results.tex", args.campaign_label)
    if args.diagnostics:
        write_diagnostics_table(args.diagnostics, args.output_dir / "multi_eos_diagnostics.tex")
    print(f"wrote multi-EoS paper assets to {args.output_dir}")


if __name__ == "__main__":
    main()
