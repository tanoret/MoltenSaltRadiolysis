"""Reactor-scale parametric molten-salt radiolysis studies."""
from __future__ import annotations
from msr_radiolysis_field_core import *

def run_msr_offgas_study(figures: Path, results: Path) -> dict[str, Any]:
    q_W_cc = np.logspace(-2, 2, 121)
    half_life_h = np.logspace(-2, 4, 121)
    Q, H = np.meshgrid(q_W_cc, half_life_h)
    V_liq = 50.0
    V_gas = 5.0
    T = 873.15
    pressure_maps: dict[str, np.ndarray] = {}
    selected_rows: list[dict[str, Any]] = []

    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.05), sharex=True, sharey=True)
    for ax, comp in zip(axes, (CHLORIDE, FLUORIDE), strict=True):
        chi = V_gas / (R_GAS * T) + comp.henry * V_liq
        total_source_mol_h = g_source(comp.stable_g_value, Q * 1.0e6) * V_liq * 3600.0
        removal_h = math.log(2.0) / H
        pressure = total_source_mol_h / (chi * removal_h)
        pressure_maps[comp.key] = pressure
        mesh = ax.pcolormesh(
            Q,
            H,
            pressure,
            shading="auto",
            norm=LogNorm(vmin=1.0, vmax=1.0e8),
            cmap="viridis",
        )
        contours = ax.contour(
            Q,
            H,
            pressure,
            levels=[10.0, 100.0, 1.0e3, 1.0e4, 1.0e5, 1.0e6],
            colors="black",
            linewidths=0.55,
        )
        ax.clabel(contours, fmt={10: "10 Pa", 100: "100 Pa", 1e3: "1 kPa", 1e4: "10 kPa", 1e5: "0.1 MPa", 1e6: "1 MPa"}, fontsize=6.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(comp.label)
        ax.set_xlabel(r"Deposited power density [W cm$^{-3}$]")
        ax.grid(which="both", alpha=0.12, linewidth=0.4)
        cbar = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.055)
        cbar.set_label("Conditional steady halogen pressure [Pa]")
        for q_sel in (0.1, 1.0, 10.0, 30.0, 100.0):
            for h_sel in (0.1, 1.0, 10.0, 100.0, 1000.0):
                source = float(g_source(comp.stable_g_value, q_sel * 1.0e6) * V_liq * 3600.0)
                p = source / (chi * (math.log(2.0) / h_sel))
                selected_rows.append(
                    {
                        "salt": comp.key,
                        "power_W_cc": q_sel,
                        "cleanup_half_life_h": h_sel,
                        "stable_G_molecules_100eV": comp.stable_g_value,
                        "source_mol_h": source,
                        "inventory_capacitance_mol_Pa": chi,
                        "steady_pressure_Pa": p,
                    }
                )
    axes[0].set_ylabel("Effective removal half-life [h]")
    fig.suptitle("Off-gas source-removal envelopes for a 50 m$^3$ salt inventory", y=1.02)
    fig.tight_layout()
    _save_figure(fig, figures / "fig_msr_offgas_operating_map")
    _write_csv(results / "msr_offgas_selected_cases.csv", selected_rows)
    return {"selected_rows": selected_rows, "pressure_maps": pressure_maps}


def run_msr_redox_study(figures: Path, results: Path) -> dict[str, Any]:
    q_W_cc = np.logspace(-2, 2, 121)
    beta = np.logspace(-10, -3, 141)
    Q, B = np.meshgrid(q_W_cc, beta)
    T = 873.15
    time_s = SECONDS_PER_YEAR
    rows: list[dict[str, Any]] = []

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.45), sharex="col", sharey="row")
    for col, comp in enumerate((CHLORIDE, FLUORIDE)):
        state = comp.composition_state()
        total_u = float(state["total_u_mol_m3"])
        f3 = float(state["nominal_u3_fraction"])
        primary = g_source(comp.radical_g_value, Q * 1.0e6)
        transfer = B * primary * time_s
        u3_0 = total_u * f3
        u4_0 = total_u - u3_0
        depleted = transfer >= 0.999 * u3_0
        transfer_clip = np.minimum(transfer, 0.999 * u3_0)
        ratio0 = u4_0 / u3_0
        ratio = (u4_0 + transfer_clip) / (u3_0 - transfer_clip)
        delta_mV = 1000.0 * R_GAS * T / FARADAY * np.log(ratio / ratio0)
        delta_plot = np.ma.masked_where(depleted, delta_mV)
        top = axes[0, col]
        im = top.pcolormesh(
            Q,
            B,
            delta_plot,
            shading="auto",
            cmap="magma",
            norm=LogNorm(vmin=0.1, vmax=500.0),
        )
        top.contourf(Q, B, depleted.astype(float), levels=[0.5, 1.5], colors=["0.82"], alpha=0.8)
        c = top.contour(Q, B, delta_mV, levels=[1, 10, 50, 100, 200], colors="black", linewidths=0.5)
        top.clabel(c, fmt=lambda x: f"{x:g} mV", fontsize=5.9)
        top.set_xscale("log")
        top.set_yscale("log")
        top.set_title(comp.label)
        cb = fig.colorbar(im, ax=top, pad=0.02, fraction=0.05)
        cb.set_label(r"$\Delta E$ after one full-power year [mV]")

        time_10_s = 0.1 * u3_0 / np.maximum(B * primary, 1.0e-300)
        time_10_y = time_10_s / SECONDS_PER_YEAR
        bottom = axes[1, col]
        im2 = bottom.pcolormesh(
            Q,
            B,
            time_10_y,
            shading="auto",
            cmap="cividis_r",
            norm=LogNorm(vmin=1.0e-3, vmax=1.0e4),
        )
        cc = bottom.contour(
            Q,
            B,
            time_10_y,
            levels=[1.0 / 365.25, 30.0 / 365.25, 1.0, 10.0, 60.0],
            colors="black",
            linewidths=0.5,
        )
        bottom.clabel(
            cc,
            fmt={1.0 / 365.25: "1 d", 30.0 / 365.25: "30 d", 1.0: "1 y", 10.0: "10 y", 60.0: "60 y"},
            fontsize=5.9,
        )
        bottom.set_xscale("log")
        bottom.set_yscale("log")
        bottom.set_xlabel(r"Deposited power density [W cm$^{-3}$]")
        cb2 = fig.colorbar(im2, ax=bottom, pad=0.02, fraction=0.05)
        cb2.set_label("Time to oxidize 10% of initial U(III) [y]")
        for q_sel in (0.1, 1.0, 10.0, 30.0):
            for beta_sel in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
                rate = float(beta_sel * g_source(comp.radical_g_value, q_sel * 1e6))
                redox = uranium_redox_state(
                    np.array([total_u]),
                    np.array([f3]),
                    np.array([T]),
                    np.array([rate]),
                    time_s,
                )
                rows.append(
                    {
                        "salt": comp.key,
                        "power_W_cc": q_sel,
                        "net_redox_branch_beta": beta_sel,
                        "initial_U3_fraction": f3,
                        "net_oxidation_rate_mol_m3_s": rate,
                        "delta_E_after_1y_mV": 1000.0 * float(redox["delta_E_V"][0]),
                        "fraction_U3_consumed_after_1y": float(redox["fraction_u3_consumed"][0]),
                        "time_to_10pct_U3_year": 0.1 * total_u * f3 / max(rate, 1e-300) / SECONDS_PER_YEAR,
                    }
                )
    axes[0, 0].set_ylabel(r"Net radiolytic redox branch, $\beta_{\rm net}$")
    axes[1, 0].set_ylabel(r"Net radiolytic redox branch, $\beta_{\rm net}$")
    for ax in axes.flat:
        ax.grid(which="both", alpha=0.10, linewidth=0.4)
    fig.suptitle("Uranium redox sensitivity to a small imbalance in primary-pair fate", y=0.998)
    fig.tight_layout()
    _save_figure(fig, figures / "fig_msr_redox_potential_map")
    _write_csv(results / "msr_redox_selected_cases.csv", rows)
    return {"selected_rows": rows}


def run_msr_intermediate_study(figures: Path, results: Path) -> dict[str, Any]:
    scavenging_lifetime_s = np.logspace(-9, -2, 121)
    q_W_cc = np.logspace(-1, 2, 101)
    TAU, Q = np.meshgrid(scavenging_lifetime_s, q_W_cc)
    temperature_K = 1000.0
    rows: list[dict[str, Any]] = []
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.35), sharex=True, sharey=True)
    for col, comp in enumerate((CHLORIDE, FLUORIDE)):
        state = comp.composition_state()
        halide = float(state["halide_mol_m3"])
        field = steady_intermediates(
            comp, temperature_K, Q * 1e6, halide, tau_e_s=TAU, tau_oxidant_s=TAU
        )
        electron_uM = field["electron_mol_m3"] * 1000.0
        oxidant_uM = field["oxidizing_intermediate_mol_m3"] * 1000.0
        labels = (
            r"$e_s^-$ [$\mu$M]",
            (r"Cl$_2^{\bullet-}$ [$\mu$M]" if comp.kernel == "chloride" else r"F$^{\bullet}$ [$\mu$M]"),
        )
        for row_idx, values in enumerate((electron_uM, oxidant_uM)):
            ax = axes[row_idx, col]
            positive = values[values > 0.0]
            vmin = max(float(np.quantile(positive, 0.005)), 1e-12)
            vmax = max(float(np.quantile(positive, 0.995)), vmin * 10.0)
            im = ax.pcolormesh(
                TAU, Q, values, shading="auto", cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax)
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(which="both", alpha=0.10, linewidth=0.4)
            cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
            cb.set_label(labels[row_idx])
            if row_idx == 0:
                ax.set_title(comp.label)
            if row_idx == 1:
                ax.set_xlabel("Unresolved scavenging lifetime [s]")
        for tau_sel in (1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2):
            for q_sel in (0.1, 1.0, 10.0, 40.0, 100.0):
                f = steady_intermediates(
                    comp, temperature_K, q_sel * 1e6, halide, tau_sel, tau_sel
                )
                rows.append(
                    {
                        "salt": comp.key,
                        "temperature_K": temperature_K,
                        "power_W_cc": q_sel,
                        "tau_e_s": tau_sel,
                        "tau_oxidant_s": tau_sel,
                        "electron_mol_m3": float(f["electron_mol_m3"]),
                        "primary_radical_mol_m3": float(f["primary_radical_mol_m3"]),
                        "oxidizing_intermediate_mol_m3": float(f["oxidizing_intermediate_mol_m3"]),
                        "polyhalide_mol_m3": float(f["polyhalide_mol_m3"]),
                        "network_stable_halogen_mol_m3_s": float(f["network_stable_halogen_mol_m3_s"]),
                        "electron_capture_fraction": float(f["electron_capture_fraction"]),
                        "oxidant_capture_fraction": float(f["oxidant_capture_fraction"]),
                        "tau_e_effective_s": float(f["tau_e_s"]),
                        "tau_oxidant_effective_s": float(f["tau_oxidant_s"]),
                    }
                )
    axes[0, 0].set_ylabel(r"Deposited power density [W cm$^{-3}$]")
    axes[1, 0].set_ylabel(r"Deposited power density [W cm$^{-3}$]")
    fig.suptitle(
        "Intermediate-species sensitivity to unresolved scavenging at 1000 K",
        y=0.998,
    )
    fig.tight_layout()
    _save_figure(fig, figures / "fig_msr_intermediate_species_map")
    _write_csv(results / "msr_intermediate_selected_cases.csv", rows)
    return {"selected_rows": rows}

