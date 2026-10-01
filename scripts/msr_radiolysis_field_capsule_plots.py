"""Publication plots for the field-resolved capsule studies."""
from __future__ import annotations
from msr_radiolysis_field_core import *
from msr_radiolysis_field_capsule_extract import _field_collection, _cell_polygons

def _case_title(data: dict[str, Any]) -> str:
    case: CapsuleCase = data["case"]
    salt_name = "chloride" if case.composition.kernel == "chloride" else "fluoride"
    return f"{case.capsule} {salt_name}"


def plot_capsule_temperature(cases: list[dict[str, Any]], figures: Path) -> None:
    all_T = np.concatenate([d["temperature_K"] for d in cases])
    norm = Normalize(vmin=float(np.min(all_T)), vmax=float(np.max(all_T)))
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 3.75), sharex=False, sharey=False)
    for panel, (ax, data) in enumerate(zip(axes.flat, cases, strict=True)):
        coll = _field_collection(ax, data["salt"], data["temperature_K"], cmap="inferno", norm=norm)
        ax.set_title(_case_title(data))
        if panel % 2 == 0:
            ax.set_ylabel("Radius [mm]")
        else:
            ax.set_ylabel("")
        if panel < 2:
            ax.set_xlabel("")
        summary = data["salt"].summary()
        ax.text(
            0.02,
            0.04,
            f"$T_V$={summary['temperature_volume_average_K']:.0f} K\n"
            f"range {summary['temperature_min_K']:.0f}--{summary['temperature_max_K']:.0f} K",
            transform=ax.transAxes,
            fontsize=6.4,
            color="white",
            ha="left",
            va="bottom",
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 2.0},
        )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.22, top=0.86, wspace=0.15, hspace=0.45)
    cax = fig.add_axes([0.20, 0.075, 0.60, 0.028])
    cbar = fig.colorbar(coll, cax=cax, orientation="horizontal")
    cbar.set_label("STAR-CCM+ salt temperature [K]", labelpad=3)
    fig.suptitle("Reconstructed saved salt-temperature fields", y=0.97)
    _save_figure(fig, figures / "fig_capsule_temperature_fields")


def plot_capsule_intermediates(cases: list[dict[str, Any]], figures: Path) -> None:
    values_all = np.concatenate(
        [d["intermediate"]["oxidizing_intermediate_mol_m3"] * 1000.0 for d in cases]
    )
    positive = values_all[values_all > 0.0]
    norm = LogNorm(vmin=float(np.quantile(positive, 0.01)), vmax=float(np.quantile(positive, 0.99)))
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 3.75))
    for panel, (ax, data) in enumerate(zip(axes.flat, cases, strict=True)):
        values = data["intermediate"]["oxidizing_intermediate_mol_m3"] * 1000.0
        coll = _field_collection(ax, data["salt"], values, cmap="viridis", norm=norm)
        ax.set_title(_case_title(data))
        if panel % 2 == 0:
            ax.set_ylabel("Radius [mm]")
        else:
            ax.set_ylabel("")
        if panel < 2:
            ax.set_xlabel("")
        label = r"Cl$_2^{\bullet-}$" if data["case"].composition.kernel == "chloride" else r"F$^{\bullet}$"
        ax.text(
            0.02,
            0.04,
            f"{label}: {volume_quantile(data['salt'], values, 0.05):.3g}--"
            f"{volume_quantile(data['salt'], values, 0.95):.3g} $\\mu$M",
            transform=ax.transAxes,
            fontsize=6.4,
            color="white",
            ha="left",
            va="bottom",
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 2.0},
        )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.22, top=0.86, wspace=0.15, hspace=0.45)
    cax = fig.add_axes([0.20, 0.075, 0.60, 0.028])
    cbar = fig.colorbar(coll, cax=cax, orientation="horizontal")
    cbar.set_label(r"Quasi-steady oxidizing intermediate [$\mu$M]", labelpad=3)
    fig.suptitle(
        r"Field-resolved intermediate concentrations ($\tau_e=\tau_{ox}=10^{-4}$ s)",
        y=0.97,
    )
    _save_figure(fig, figures / "fig_capsule_intermediate_fields")


def plot_capsule_stable_sources(cases: list[dict[str, Any]], figures: Path) -> None:
    """Plot the gross empirical stable-halogen source before downstream losses."""
    values_all = np.concatenate(
        [d["empirical_stable_source_mol_m3_s"] for d in cases]
    )
    positive = values_all[values_all > 0.0]
    norm = LogNorm(
        vmin=float(np.quantile(positive, 0.005)),
        vmax=float(np.quantile(positive, 0.995)),
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 3.75))
    for panel, (ax, data) in enumerate(zip(axes.flat, cases, strict=True)):
        values = data["empirical_stable_source_mol_m3_s"]
        coll = _field_collection(ax, data["salt"], values, cmap="plasma", norm=norm)
        ax.set_title(_case_title(data))
        if panel % 2 == 0:
            ax.set_ylabel("Radius [mm]")
        else:
            ax.set_ylabel("")
        if panel < 2:
            ax.set_xlabel("")
        avg = volume_average(data["salt"], values)
        p05 = volume_quantile(data["salt"], values, 0.05)
        p95 = volume_quantile(data["salt"], values, 0.95)
        ax.text(
            0.02,
            0.04,
            rf"$\langle S_{{X_2}}\rangle_V$={avg:.2e} mol m$^{{-3}}$ s$^{{-1}}$" "\n"
            rf"P05--P95: {p05:.2e}--{p95:.2e}",
            transform=ax.transAxes,
            fontsize=6.1,
            color="white",
            ha="left",
            va="bottom",
            bbox={"facecolor": "black", "alpha": 0.48, "edgecolor": "none", "pad": 2.0},
        )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.22, top=0.86, wspace=0.15, hspace=0.45)
    cax = fig.add_axes([0.20, 0.075, 0.60, 0.028])
    cbar = fig.colorbar(coll, cax=cax, orientation="horizontal")
    cbar.set_label(r"Gross empirical stable-halogen source, $S_{X_2}$ [mol m$^{-3}$ s$^{-1}$]", labelpad=3)
    fig.suptitle(
        "Stable-halogen source field before redox buffering, wall loss, and venting",
        y=0.97,
    )
    _save_figure(fig, figures / "fig_capsule_stable_source_fields")


def plot_capsule_redox(cases: list[dict[str, Any]], figures: Path) -> None:
    values_all = np.concatenate([1000.0 * d["redox4"]["delta_E_V"] for d in cases])
    max_abs = max(float(np.quantile(np.abs(values_all), 0.99)), 1.0)
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 3.75))
    for panel, (ax, data) in enumerate(zip(axes.flat, cases, strict=True)):
        values = 1000.0 * data["redox4"]["delta_E_V"]
        coll = _field_collection(ax, data["salt"], values, cmap="coolwarm", norm=norm)
        ax.set_title(_case_title(data))
        if panel % 2 == 0:
            ax.set_ylabel("Radius [mm]")
        else:
            ax.set_ylabel("")
        if panel < 2:
            ax.set_xlabel("")
        avg = volume_average(data["salt"], values)
        p05 = volume_quantile(data["salt"], values, 0.05)
        p95 = volume_quantile(data["salt"], values, 0.95)
        ax.text(
            0.02,
            0.04,
            f"$\\langle\\Delta E\\rangle_V$={avg:.1f} mV\nP05--P95: {p05:.1f}--{p95:.1f} mV",
            transform=ax.transAxes,
            fontsize=6.4,
            color="black",
            ha="left",
            va="bottom",
            bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 2.0},
        )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.22, top=0.86, wspace=0.15, hspace=0.45)
    cax = fig.add_axes([0.20, 0.075, 0.60, 0.028])
    cbar = fig.colorbar(coll, cax=cax, orientation="horizontal")
    cbar.set_label(r"Four-cycle relative U(IV)/U(III) shift, $\Delta E$ [mV]", labelpad=3)
    fig.suptitle(
        r"Conditional redox field for $\langle\beta_{net}\rangle_V=10^{-5}$",
        y=0.97,
    )
    _save_figure(fig, figures / "fig_capsule_redox_fields")


def plot_capsule_fast_evolution(cases: list[dict[str, Any]], figures: Path) -> None:
    times = np.logspace(-13, 0, 450)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05))
    for data in cases:
        mesh = data["salt"]
        inter = data["intermediate"]
        e_history = approach_to_quasi_steady(inter["electron_mol_m3"], inter["tau_e_s"], times)
        x_history = approach_to_quasi_steady(
            inter["oxidizing_intermediate_mol_m3"], inter["tau_oxidant_s"], times
        )
        weights = mesh.cell_volume_m3 / mesh.total_volume_m3
        e_avg = np.sum(e_history * weights[:, None], axis=0) * 1000.0
        x_avg = np.sum(x_history * weights[:, None], axis=0) * 1000.0
        label = _case_title(data)
        axes[0].plot(times, e_avg, label=label)
        axes[1].plot(times, x_avg, label=label)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Time after irradiation begins [s]")
        ax.grid(which="both", alpha=0.16, linewidth=0.45)
    axes[0].set_ylabel(r"Volume-averaged $e_s^-$ [$\mu$M]")
    axes[1].set_ylabel(r"Volume-averaged oxidizing intermediate [$\mu$M]")
    axes[1].legend(loc="lower right", ncol=1)
    fig.suptitle("Fast establishment of the frozen-field radiolysis state", y=1.01)
    fig.tight_layout()
    _save_figure(fig, figures / "fig_capsule_fast_species_evolution")


def plot_capsule_long_evolution(cases: list[dict[str, Any]], figures: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
    markers = {"3": "o", "4": "s"}
    for data in cases:
        case: CapsuleCase = data["case"]
        t_end = max(data["time4_s"], data["time3_s"])
        time = np.linspace(0.0, t_end, 320)
        days = time / SECONDS_PER_DAY
        gross = data["stable_source_total_mol_s"] * time
        axes[0].plot(days, gross, label=_case_title(data))
        axes[0].scatter(
            [data["time3_s"] / SECONDS_PER_DAY, data["time4_s"] / SECONDS_PER_DAY],
            [case.gas_3_mol, case.gas_4_mol],
            marker=markers["3"],
            s=18,
            facecolors="none",
            edgecolors=axes[0].lines[-1].get_color(),
            linewidths=0.9,
            zorder=4,
        )
        # Redox evolution is evaluated on a common cell field at each time.
        avg_delta = []
        for t in time:
            redox = uranium_redox_state(
                data["composition_state"]["total_u_mol_m3"],
                data["composition_state"]["nominal_u3_fraction"],
                data["temperature_K"],
                data["net_redox_rate_mol_m3_s"],
                float(t),
            )
            avg_delta.append(1000.0 * volume_average(data["salt"], redox["delta_E_V"]))
        axes[1].plot(days, avg_delta, label=_case_title(data))
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Equivalent full-power irradiation time [d]")
    axes[0].set_ylabel("Inventory [mol]")
    axes[0].grid(which="both", alpha=0.16, linewidth=0.45)
    axes[0].text(
        0.02,
        0.03,
        "Lines: gross empirical Cl$_2$/F$_2$ source\nOpen points: reported total gas (composition unspecified)",
        transform=axes[0].transAxes,
        fontsize=6.4,
        va="bottom",
    )
    axes[1].set_xlabel("Equivalent full-power irradiation time [d]")
    axes[1].set_ylabel(r"Volume-averaged relative redox shift, $\Delta E$ [mV]")
    axes[1].grid(which="both", alpha=0.16, linewidth=0.45)
    axes[1].legend(loc="best")
    fig.suptitle("Slow inventory and redox evolution in the sealed capsules", y=1.01)
    fig.tight_layout()
    _save_figure(fig, figures / "fig_capsule_long_term_evolution")

