#!/usr/bin/env python3
"""
Generate publication-quality and presentation-ready DAG / Bayesian Hierarchical Model diagrams.
Outputs:
  - results/figures/figure_dag_model.png (Presentation Infographic, 300 DPI)
  - results/figures/figure_dag_model.pdf
  - results/figures/figure_dag_model.svg
  - results/figures/figure_dag_plate.png (Formal Publication Plate DAG, 300 DPI)
  - results/figures/figure_dag_plate.pdf
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def create_presentation_dag():
    """Generates a modern, elegant, card-based DAG optimized for slide presentations and executive reading."""
    fig, ax = plt.subplots(figsize=(15.2, 9.8), dpi=300)
    ax.set_xlim(0, 15.2)
    ax.set_ylim(0, 9.8)
    ax.axis("off")

    fig.patch.set_facecolor("#FFFFFF")
    ax.set_facecolor("#FFFFFF")

    # Clean Color Palette
    c_banner_bg = "#0F172A"      # Slate 900
    c_hyper_bg = "#F8FAFC"       # Slate 50
    c_hyper_border = "#94A3B8"   # Slate 400
    c_param_bg = "#F0F7FF"       # Soft Blue
    c_param_border = "#2563EB"   # Blue 600
    c_pred_bg = "#ECFDF5"        # Soft Emerald
    c_pred_border = "#059669"    # Emerald 600
    c_disp_bg = "#FFFBEB"        # Soft Amber
    c_disp_border = "#D97706"    # Amber 600
    c_obs_bg = "#FAF5FF"         # Soft Violet
    c_obs_border = "#7C3AED"     # Violet 600
    c_arrow = "#334155"          # Slate 700

    def draw_card(x, y, w, h, bg_color, border_color, category=None, cat_bg="#64748B", title="", subtitle=None):
        """Draws a card with a clean category badge at top-left and title underneath, preventing any collisions."""
        card = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.06,rounding_size=0.16",
            facecolor=bg_color,
            edgecolor=border_color,
            linewidth=1.8,
            zorder=2
        )
        ax.add_patch(card)

        curr_y = y + h - 0.28
        # Category Badge
        if category:
            badge_w = len(category) * 0.11 + 0.45
            b_patch = FancyBboxPatch(
                (x + 0.28, curr_y - 0.22), badge_w, 0.25,
                boxstyle="round,pad=0.02,rounding_size=0.07",
                facecolor=cat_bg,
                edgecolor="none",
                zorder=4
            )
            ax.add_patch(b_patch)
            ax.text(
                x + 0.28 + badge_w / 2.0, curr_y - 0.09,
                category,
                fontsize=7.2,
                fontweight="bold",
                color="#FFFFFF",
                ha="center",
                va="center",
                zorder=5
            )
            curr_y -= 0.38

        # Title
        if title:
            ax.text(
                x + 0.28, curr_y,
                title,
                fontsize=11.0,
                fontweight="bold",
                color="#0F172A",
                va="center",
                zorder=4
            )
            curr_y -= 0.32

        # Optional Subtitle
        if subtitle:
            ax.text(
                x + 0.28, curr_y,
                subtitle,
                fontsize=8.2,
                color="#64748B",
                va="center",
                zorder=4
            )

    # -------------------------------------------------------------
    # 0. HEADER BANNER
    # -------------------------------------------------------------
    banner = FancyBboxPatch(
        (0.6, 8.75), 14.0, 0.85,
        boxstyle="round,pad=0.06,rounding_size=0.14",
        facecolor=c_banner_bg,
        edgecolor="#020617",
        linewidth=1.2,
        zorder=2
    )
    ax.add_patch(banner)
    ax.text(
        7.6, 9.3,
        "GENERATIVE DIRECTED ACYCLIC GRAPH (DAG) OF BAYESIAN MODEL",
        fontsize=13.5,
        fontweight="bold",
        color="#F8FAFC",
        ha="center",
        va="center",
        zorder=3
    )
    ax.text(
        7.6, 8.97,
        "Hierarchical Formulation for Brain Morphometry under Clinical Slice Degradation & Gadolinium Confounding",
        fontsize=9.2,
        color="#94A3B8",
        ha="center",
        va="center",
        zorder=3
    )

    # -------------------------------------------------------------
    # 1. LAYER 1: HYPERPRIORS
    # -------------------------------------------------------------
    # Disease Hyperprior
    draw_card(1.0, 7.15, 5.8, 1.4, c_hyper_bg, c_hyper_border,
              category="HYPERPRIOR", cat_bg="#64748B",
              title="Disease Scale Hyperprior",
              subtitle="Controls magnitude of diagnostic phenotype offsets")
    ax.text(1.3, 7.35, r"$\tau_\beta \sim \mathrm{Half\text{-}Normal}(0, \, 0.05)$",
            fontsize=12.2, color="#0F172A", zorder=4)

    # Site Hyperprior
    draw_card(8.4, 7.15, 5.8, 1.4, c_hyper_bg, c_hyper_border,
              category="HYPERPRIOR", cat_bg="#64748B",
              title="Scanner Site Scale Hyperprior",
              subtitle="Controls inter-scanner calibration variance across 6 hospitals")
    ax.text(8.7, 7.35, r"$\tau_\gamma \sim \mathrm{Half\text{-}Normal}(0, \, 0.05)$",
            fontsize=12.2, color="#0F172A", zorder=4)

    # -------------------------------------------------------------
    # 2. LAYER 2: NON-CENTERED RANDOM EFFECTS
    # -------------------------------------------------------------
    # Disease Effects
    draw_card(1.0, 5.0, 5.8, 1.8, c_param_bg, c_param_border,
              category="GROUP EFFECT (NON-CENTERED)", cat_bg="#2563EB",
              title="Disease Diagnostic Offsets (5 Classes)",
              subtitle=r"Baseline: $\beta_{\mathrm{CONTROL}} \equiv 0$ | Partial pooling for small groups")
    ax.text(1.3, 5.5,
            r"$z_\beta \sim \mathcal{N}(0, 1) \quad \longrightarrow \quad \beta_d = \tau_\beta \cdot z_\beta$",
            fontsize=11.5, fontweight="bold", color="#1E3A8A", zorder=4)
    ax.text(1.3, 5.18,
            r"$\beta_d \in \{\mathrm{CONTROL}, \mathrm{DEMENTIA}, \mathrm{EPILEPSY}, \mathrm{HYDROCEPHALUS}, \mathrm{PARKINSON}\}$",
            fontsize=8.0, color="#1D4ED8", zorder=4)

    # Site Intercepts
    draw_card(8.4, 5.0, 5.8, 1.8, c_param_bg, c_param_border,
              category="RANDOM INTERCEPT (NON-CENTERED)", cat_bg="#2563EB",
              title="Scanner Site Intercepts (6 Centers)",
              subtitle=r"Absorbs regional $B_0$ field and RF coil variations")
    ax.text(8.7, 5.5,
            r"$z_\gamma \sim \mathcal{N}(0, 1) \quad \longrightarrow \quad \gamma_s = \tau_\gamma \cdot z_\gamma$",
            fontsize=11.5, fontweight="bold", color="#1E3A8A", zorder=4)
    ax.text(8.7, 5.18,
            r"$\gamma_s \in \{\mathrm{RSUTH}, \mathrm{UPTH}, \mathrm{IDC}, \mathrm{AKTH}, \mathrm{BMH}, \mathrm{LifeBridge}\}$",
            fontsize=8.0, color="#1D4ED8", zorder=4)

    # Flanking Fixed Intercept & Contrast Effects
    draw_card(0.6, 3.1, 2.5, 1.5, "#F8FAFC", "#CBD5E1",
              category="FIXED PRIOR", cat_bg="#64748B",
              title="Global Baseline",
              subtitle="Expected 1mm Control")
    ax.text(0.88, 3.32, r"$\alpha \sim \mathcal{N}(0, 1)$", fontsize=11.5, fontweight="bold", color="#334155", zorder=4)

    draw_card(12.1, 3.1, 2.5, 1.5, "#F8FAFC", "#CBD5E1",
              category="FIXED PRIOR", cat_bg="#64748B",
              title="Gadolinium Shift",
              subtitle=r"$\text{Contrast}_i \in \{0, 1\}$")
    ax.text(12.38, 3.42, r"$\delta \cdot \mathrm{Contrast}_i$", fontsize=11, fontweight="bold", color="#334155", zorder=4)
    ax.text(12.38, 3.2, r"$\delta \sim \mathcal{N}(0, 0.1)$", fontsize=8.2, color="#64748B", zorder=4)

    # -------------------------------------------------------------
    # 3. LAYER 3: LINEAR PREDICTOR
    # -------------------------------------------------------------
    draw_card(3.5, 3.05, 8.2, 1.65, c_pred_bg, c_pred_border,
              category="DETERMINISTIC LINK", cat_bg="#059669",
              title="Linear Predictor (Expected Metric Value per Subject)",
              subtitle="Integrates baseline geometry + disease deviation + scanner shift + contrast enhancement")
    ax.text(7.6, 3.45,
            r"$\mu_i = \alpha + \beta_{\mathrm{diagnosis}[i]} + \gamma_{\mathrm{site}[i]} + \delta \cdot \mathrm{Contrast}_i$",
            fontsize=13.0, fontweight="bold", color="#065F46", ha="center", zorder=4)

    # -------------------------------------------------------------
    # 4. LAYER 4: HETEROSKEDASTIC DISPERSION
    # -------------------------------------------------------------
    draw_card(8.4, 0.95, 5.8, 1.75, c_disp_bg, c_disp_border,
              category="SCALE MODEL (MEASUREMENT NOISE)", cat_bg="#D97706",
              title="Heteroskedastic Dispersion",
              subtitle="Exponential slice thickness blur elasticity")
    ax.text(8.7, 1.48,
            r"$\sigma_i = \sigma_0 \cdot \exp\left( \lambda \cdot [\mathrm{SliceThickness}_i - 1.0] \right)$",
            fontsize=11.0, fontweight="bold", color="#92400E", zorder=4)
    ax.text(8.7, 1.12,
            r"$\sigma_0 \sim \mathrm{Half\text{-}Normal}(0, 0.05) \quad | \quad \lambda \sim \mathcal{N}(0, 0.15)$",
            fontsize=8.8, color="#B45309", zorder=4)

    # -------------------------------------------------------------
    # 5. LAYER 5: LIKELIHOOD / OBSERVED DATA (With Subject Plate)
    # -------------------------------------------------------------
    plate_i = FancyBboxPatch(
        (0.8, 0.25), 6.8, 2.55,
        boxstyle="round,pad=0.08,rounding_size=0.16",
        facecolor="#FDF4FF",
        edgecolor="#C084FC",
        linewidth=1.6,
        linestyle="--",
        zorder=1
    )
    ax.add_patch(plate_i)
    ax.text(1.08, 2.62, r"Subject Plate:  $i \in \{1, \dots, N\} \quad (N = 218)$",
            fontsize=9.2, fontweight="bold", color="#9333EA", zorder=3)

    # Likelihood Card
    draw_card(1.1, 0.45, 6.2, 1.95, c_obs_bg, c_obs_border,
              category="OBSERVED DATA / LIKELIHOOD", cat_bg="#7C3AED",
              title="Gaussian Likelihood",
              subtitle="Evans' Index (EI), Parenchymal Envelope (PEF), or VBR")
    ax.text(4.2, 1.15,
            r"$y_i \sim \mathcal{N}\left( \, \mu_i, \,\, \sigma_i^2 \, \right)$",
            fontsize=16, fontweight="bold", color="#4C1D95", ha="center", zorder=4)
    ax.text(4.2, 0.72,
            r"Heteroskedastic variance $\sigma_i^2$ penalizes low-confidence thick slices",
            fontsize=8.2, color="#6D28D9", ha="center", zorder=4)

    # -------------------------------------------------------------
    # ARROWS (Directed Generative Flow)
    # -------------------------------------------------------------
    def draw_arrow(start, end, rad=0.0):
        arrow = FancyArrowPatch(
            start, end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.6,
            color=c_arrow,
            connectionstyle=f"arc3,rad={rad}",
            zorder=6
        )
        ax.add_patch(arrow)

    # Hyperpriors down to random effects
    draw_arrow((3.9, 7.15), (3.9, 6.8))
    draw_arrow((11.3, 7.15), (11.3, 6.8))

    # Random effects down to linear predictor
    draw_arrow((3.9, 5.0), (5.5, 4.7))
    draw_arrow((11.3, 5.0), (9.7, 4.7))

    # Baseline & Contrast into linear predictor
    draw_arrow((3.1, 3.85), (3.5, 3.85))
    draw_arrow((12.1, 3.85), (11.7, 3.85))

    # Linear Predictor down to Likelihood
    draw_arrow((6.0, 3.05), (4.5, 2.4))

    # Heteroskedastic Dispersion into Likelihood
    draw_arrow((8.4, 1.8), (7.3, 1.45))

    out_png = FIGURES_DIR / "figure_dag_model.png"
    out_pdf = FIGURES_DIR / "figure_dag_model.pdf"
    out_svg = FIGURES_DIR / "figure_dag_model.svg"

    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="#FFFFFF", edgecolor="none")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ Generated Presentation DAG: {out_png}")


def create_publication_plate_dag():
    """Generates a formal, standard Bayesian Plate Notation DAG for journal publications."""
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8.5)
    ax.axis("off")

    fig.patch.set_facecolor("#FFFFFF")
    ax.set_facecolor("#FFFFFF")

    # Helper for node creation
    def draw_node(x, y, r, label, sublabel=None, observed=False, deterministic=False):
        if deterministic:
            box = FancyBboxPatch(
                (x - r, y - r * 0.75), r * 2, r * 1.5,
                boxstyle="round,pad=0.04,rounding_size=0.1",
                facecolor="#F0FDF4",
                edgecolor="#16A34A",
                linewidth=1.8,
                zorder=3
            )
            ax.add_patch(box)
        else:
            circle = Circle(
                (x, y), r,
                facecolor="#E2E8F0" if observed else "#FFFFFF",
                edgecolor="#0F172A",
                linewidth=1.8,
                zorder=3
            )
            ax.add_patch(circle)

        ax.text(x, y + (0.05 if sublabel else 0), label, fontsize=12, fontweight="bold",
                ha="center", va="center", color="#0F172A", zorder=4)
        if sublabel:
            ax.text(x, y - 0.28, sublabel, fontsize=7.5, fontstyle="italic",
                    ha="center", va="center", color="#475569", zorder=4)

    # 1. Subject Plate [N = 218]
    plate_subj = FancyBboxPatch(
        (1.0, 0.6), 10.0, 4.6,
        boxstyle="round,pad=0.08,rounding_size=0.12",
        facecolor="#F8FAFC",
        edgecolor="#64748B",
        linewidth=1.4,
        linestyle="-",
        zorder=1
    )
    ax.add_patch(plate_subj)
    ax.text(10.7, 0.85, r"$i \in \{1, \dots, N\}$", fontsize=10, fontweight="bold",
            ha="right", color="#334155", zorder=2)

    # 2. Site Plate [S = 6]
    plate_site = FancyBboxPatch(
        (6.8, 5.4), 3.2, 2.5,
        boxstyle="round,pad=0.08,rounding_size=0.12",
        facecolor="#F1F5F9",
        edgecolor="#94A3B8",
        linewidth=1.2,
        linestyle="--",
        zorder=1
    )
    ax.add_patch(plate_site)
    ax.text(9.8, 5.65, r"$s \in \{1, \dots, 6\}$", fontsize=9, fontweight="bold",
            ha="right", color="#475569", zorder=2)

    # 3. Diagnosis Plate [D = 5]
    plate_diag = FancyBboxPatch(
        (2.0, 5.4), 3.2, 2.5,
        boxstyle="round,pad=0.08,rounding_size=0.12",
        facecolor="#F1F5F9",
        edgecolor="#94A3B8",
        linewidth=1.2,
        linestyle="--",
        zorder=1
    )
    ax.add_patch(plate_diag)
    ax.text(5.0, 5.65, r"$d \in \{1, \dots, 5\}$", fontsize=9, fontweight="bold",
            ha="right", color="#475569", zorder=2)

    # Nodes
    # Hyperpriors
    draw_node(3.6, 7.3, 0.45, r"$\tau_\beta$", "disease scale")
    draw_node(8.4, 7.3, 0.45, r"$\tau_\gamma$", "site scale")

    # Random effects
    draw_node(3.6, 6.1, 0.45, r"$\beta_d$", "disease effect")
    draw_node(8.4, 6.1, 0.45, r"$\gamma_s$", "site intercept")

    # Baseline & Fixed effects inside / entering plate
    draw_node(1.8, 4.2, 0.45, r"$\alpha$", "intercept")
    draw_node(3.8, 4.2, 0.45, r"$\delta$", "contrast coef")
    draw_node(5.2, 4.2, 0.40, r"$C_i$", "contrast", observed=True)

    # Deterministic linear predictor
    draw_node(3.6, 2.6, 0.65, r"$\mu_i$", r"$\alpha + \beta_{d[i]} + \gamma_{s[i]} + \delta C_i$", deterministic=True)

    # Heteroskedastic Dispersion nodes
    draw_node(7.8, 3.8, 0.45, r"$\sigma_0$", "base noise")
    draw_node(9.6, 3.8, 0.45, r"$\lambda$", "elasticity")
    draw_node(10.0, 2.6, 0.40, r"$T_i$", "slice thick.", observed=True)
    draw_node(8.4, 2.6, 0.65, r"$\sigma_i$", r"$\sigma_0 \exp(\lambda [T_i - 1])$", deterministic=True)

    # Observed Likelihood Node
    draw_node(6.0, 1.4, 0.60, r"$y_i$", "observed metric", observed=True)

    # Directed Arrows
    def draw_arrow(start, end):
        arrow = FancyArrowPatch(
            start, end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.4,
            color="#334155",
            zorder=5
        )
        ax.add_patch(arrow)

    # Flows
    draw_arrow((3.6, 6.85), (3.6, 6.55))
    draw_arrow((8.4, 6.85), (8.4, 6.55))

    draw_arrow((3.6, 5.65), (3.6, 3.1))
    draw_arrow((8.4, 5.65), (4.2, 3.1))

    draw_arrow((1.8, 3.75), (3.0, 2.8))
    draw_arrow((3.8, 3.75), (3.7, 3.1))
    draw_arrow((5.2, 3.8), (4.1, 3.0))

    draw_arrow((7.8, 3.35), (8.1, 3.0))
    draw_arrow((9.6, 3.35), (8.7, 3.0))
    draw_arrow((9.6, 2.6), (9.05, 2.6))

    draw_arrow((3.6, 2.1), (5.5, 1.6))
    draw_arrow((8.4, 2.1), (6.5, 1.6))

    # Title
    ax.text(6.0, 8.2, "Formal Bayesian Plate Notation: Brain Morphometry Generative DAG",
            fontsize=13, fontweight="bold", ha="center", color="#0F172A")

    out_png = FIGURES_DIR / "figure_dag_plate.png"
    out_pdf = FIGURES_DIR / "figure_dag_plate.pdf"

    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="#FFFFFF")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ Generated Publication Plate DAG: {out_png}")


if __name__ == "__main__":
    create_presentation_dag()
    create_publication_plate_dag()
