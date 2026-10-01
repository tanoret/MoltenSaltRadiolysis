"""Kinetic and inventory closures for field-resolved radiolysis studies."""
from __future__ import annotations
from msr_radiolysis_field_base import *

def steady_intermediates(
    composition: SaltComposition,
    temperature_K: float | np.ndarray,
    power_density_W_m3: float | np.ndarray,
    halide_mol_m3: float | np.ndarray,
    tau_e_s: float = TAU_E_REF_S,
    tau_oxidant_s: float = TAU_OX_REF_S,
) -> dict[str, np.ndarray]:
    """Return a local quasi-steady solution of the repository reaction kernel.

    The base kernels lack sufficient composition-specific sinks to bound all
    primary species in uranium-bearing salts.  Two explicit pseudo-first-order
    closures are therefore added: ``tau_e_s`` for reducing-electron capture and
    ``tau_oxidant_s`` for capture of the dominant oxidizing intermediate.  The
    returned capture fractions quantify how strongly those closures affect the
    result and are varied in the parametric study.
    """
    tau_e_input = np.asarray(tau_e_s, dtype=float)
    tau_ox_input = np.asarray(tau_oxidant_s, dtype=float)
    if np.any(tau_e_input <= 0.0) or np.any(tau_ox_input <= 0.0):
        raise ValueError("scavenging lifetimes must be positive")
    T, q, halide, tau_e_input, tau_ox_input = np.broadcast_arrays(
        np.asarray(temperature_K, float),
        np.asarray(power_density_W_m3, float),
        np.asarray(halide_mol_m3, float),
        tau_e_input,
        tau_ox_input,
    )
    k_e_scav = 1.0 / tau_e_input
    k_ox_scav = 1.0 / tau_ox_input

    if composition.kernel == "chloride":
        S_e = g_source(G_E_CHLORIDE, q)
        S_ox = g_source(G_CL_RADICAL, q)
        k_capture = arrhenius_from_reference(1.0e7, 673.15, 2.0e4, T)
        k_recomb = arrhenius_from_reference(1.0e7, 673.15, 2.5e4, T)
        k_radical_pair = arrhenius_from_reference(5.0e6, 673.15, 2.0e4, T)
        k_disprop = arrhenius_from_reference(2.2e6, 673.15, 2.6e4, T)
        k_breakup = arrhenius_from_reference(1.0, 673.15, 1.0e4, T)

        electron = S_e / k_e_scav
        chlorine_radical = np.zeros_like(S_ox)
        for _ in range(80):
            linear_cl = k_capture * halide + k_recomb * electron
            new_cl = _positive_quadratic_root(S_ox, linear_cl, 2.0 * k_radical_pair)
            new_e = np.divide(
                S_e,
                k_e_scav + k_recomb * new_cl,
                out=np.zeros_like(S_e),
                where=(k_e_scav + k_recomb * new_cl) > 0.0,
            )
            if np.max(np.abs(new_e - electron) / np.maximum(new_e, 1.0e-30)) < 1.0e-11:
                electron = new_e
                chlorine_radical = new_cl
                break
            electron = 0.5 * electron + 0.5 * new_e
            chlorine_radical = new_cl

        cl2_radical_source = k_capture * halide * chlorine_radical
        cl2_radical = _positive_quadratic_root(
            cl2_radical_source,
            np.full_like(T, k_ox_scav),
            2.0 * k_disprop,
        )
        cl3 = np.divide(
            k_disprop * cl2_radical**2,
            k_breakup,
            out=np.zeros_like(cl2_radical),
            where=k_breakup > 0.0,
        )
        stable_rate = k_radical_pair * chlorine_radical**2 + k_breakup * cl3
        electron_capture_fraction = np.divide(
            k_e_scav * electron,
            S_e,
            out=np.zeros_like(S_e),
            where=S_e > 0.0,
        )
        oxidant_capture_fraction = np.divide(
            k_ox_scav * cl2_radical,
            S_ox,
            out=np.zeros_like(S_ox),
            where=S_ox > 0.0,
        )
        tau_e = 1.0 / np.maximum(k_e_scav + k_recomb * chlorine_radical, 1.0e-300)
        tau_primary = 1.0 / np.maximum(
            k_capture * halide + k_recomb * electron + 4.0 * k_radical_pair * chlorine_radical,
            1.0e-300,
        )
        tau_oxidant = 1.0 / np.maximum(k_ox_scav + 4.0 * k_disprop * cl2_radical, 1.0e-300)
        tau_polyhalide = 1.0 / np.maximum(k_breakup, 1.0e-300)
        return {
            "source_pair_mol_m3_s": S_ox,
            "electron_mol_m3": electron,
            "primary_radical_mol_m3": chlorine_radical,
            "oxidizing_intermediate_mol_m3": cl2_radical,
            "polyhalide_mol_m3": cl3,
            "network_stable_halogen_mol_m3_s": stable_rate,
            "electron_capture_fraction": np.clip(electron_capture_fraction, 0.0, 1.0),
            "oxidant_capture_fraction": np.clip(oxidant_capture_fraction, 0.0, 1.0),
            "tau_e_s": tau_e,
            "tau_primary_s": tau_primary,
            "tau_oxidant_s": tau_oxidant,
            "tau_polyhalide_s": tau_polyhalide,
        }

    S_e = g_source(G_E_FLUORIDE, q)
    S_ox = g_source(G_F_RADICAL, q)
    k_recomb = arrhenius_from_reference(1.0e7, 673.15, 2.5e4, T)
    k_pair = arrhenius_from_reference(5.0e6, 673.15, 2.0e4, T)
    electron = S_e / k_e_scav
    fluorine_radical = np.sqrt(np.maximum(S_ox / (2.0 * k_pair), 0.0))
    for _ in range(80):
        linear_f = k_ox_scav + k_recomb * electron
        new_f = _positive_quadratic_root(S_ox, linear_f, 2.0 * k_pair)
        new_e = np.divide(
            S_e,
            k_e_scav + k_recomb * new_f,
            out=np.zeros_like(S_e),
            where=(k_e_scav + k_recomb * new_f) > 0.0,
        )
        if np.max(np.abs(new_e - electron) / np.maximum(new_e, 1.0e-30)) < 1.0e-11:
            electron = new_e
            fluorine_radical = new_f
            break
        electron = 0.5 * electron + 0.5 * new_e
        fluorine_radical = new_f
    stable_rate = k_pair * fluorine_radical**2
    electron_capture_fraction = np.divide(
        k_e_scav * electron,
        S_e,
        out=np.zeros_like(S_e),
        where=S_e > 0.0,
    )
    oxidant_capture_fraction = np.divide(
        k_ox_scav * fluorine_radical,
        S_ox,
        out=np.zeros_like(S_ox),
        where=S_ox > 0.0,
    )
    tau_e = 1.0 / np.maximum(k_e_scav + k_recomb * fluorine_radical, 1.0e-300)
    tau_primary = 1.0 / np.maximum(
        k_ox_scav + k_recomb * electron + 4.0 * k_pair * fluorine_radical,
        1.0e-300,
    )
    return {
        "source_pair_mol_m3_s": S_ox,
        "electron_mol_m3": electron,
        "primary_radical_mol_m3": fluorine_radical,
        "oxidizing_intermediate_mol_m3": fluorine_radical,
        "polyhalide_mol_m3": np.zeros_like(fluorine_radical),
        "network_stable_halogen_mol_m3_s": stable_rate,
        "electron_capture_fraction": np.clip(electron_capture_fraction, 0.0, 1.0),
        "oxidant_capture_fraction": np.clip(oxidant_capture_fraction, 0.0, 1.0),
        "tau_e_s": tau_e,
        "tau_primary_s": tau_primary,
        "tau_oxidant_s": tau_primary,
        "tau_polyhalide_s": np.zeros_like(tau_primary),
    }


def intermediate_balance_residuals(
    composition: SaltComposition,
    temperature_K: float,
    power_density_W_m3: float,
    halide_mol_m3: float,
    tau_e_s: float = TAU_E_REF_S,
    tau_oxidant_s: float = TAU_OX_REF_S,
) -> dict[str, float]:
    """Return normalized steady-state residuals for regression testing."""
    result = steady_intermediates(
        composition,
        temperature_K,
        power_density_W_m3,
        halide_mol_m3,
        tau_e_s,
        tau_oxidant_s,
    )
    T = float(temperature_K)
    q = float(power_density_W_m3)
    e = float(result["electron_mol_m3"])
    r = float(result["primary_radical_mol_m3"])
    x = float(result["oxidizing_intermediate_mol_m3"])
    ke = 1.0 / tau_e_s
    kx = 1.0 / tau_oxidant_s
    if composition.kernel == "chloride":
        S = float(g_source(G_CL_RADICAL, q))
        k1 = float(arrhenius_from_reference(1.0e7, 673.15, 2.0e4, T))
        k2 = float(arrhenius_from_reference(1.0e7, 673.15, 2.5e4, T))
        k3 = float(arrhenius_from_reference(5.0e6, 673.15, 2.0e4, T))
        k4 = float(arrhenius_from_reference(2.2e6, 673.15, 2.6e4, T))
        source_x = k1 * halide_mol_m3 * r
        return {
            "electron": (S - ke * e - k2 * e * r) / max(S, 1e-30),
            "primary": (
                S - k1 * halide_mol_m3 * r - k2 * e * r - 2.0 * k3 * r * r
            )
            / max(S, 1e-30),
            "oxidant": (source_x - kx * x - 2.0 * k4 * x * x) / max(source_x, 1e-30),
        }
    S = float(g_source(G_F_RADICAL, q))
    k1 = float(arrhenius_from_reference(1.0e7, 673.15, 2.5e4, T))
    k2 = float(arrhenius_from_reference(5.0e6, 673.15, 2.0e4, T))
    return {
        "electron": (S - ke * e - k1 * e * r) / max(S, 1e-30),
        "primary": (S - kx * r - k1 * e * r - 2.0 * k2 * r * r) / max(S, 1e-30),
        "oxidant": (S - kx * x - k1 * e * x - 2.0 * k2 * x * x) / max(S, 1e-30),
    }


def approach_to_quasi_steady(steady: np.ndarray, tau_s: np.ndarray, time_s: float | np.ndarray) -> np.ndarray:
    t = np.asarray(time_s, dtype=float)
    if np.any(t < 0.0):
        raise ValueError("time must be nonnegative")
    return np.asarray(steady)[..., None] * (
        1.0 - np.exp(-t.reshape((1,) * np.asarray(steady).ndim + t.shape) / np.maximum(np.asarray(tau_s)[..., None], 1e-300))
    )


def uranium_redox_state(
    total_u_mol_m3: np.ndarray,
    initial_u3_fraction: np.ndarray,
    temperature_K: np.ndarray,
    net_oxidizing_equivalent_rate_mol_m3_s: np.ndarray,
    time_s: float,
) -> dict[str, np.ndarray]:
    """Map an equivalent balance into a relative U(IV)/U(III) Nernst shift."""
    total_u, f3, T, rate = np.broadcast_arrays(
        np.asarray(total_u_mol_m3, float),
        np.asarray(initial_u3_fraction, float),
        np.asarray(temperature_K, float),
        np.asarray(net_oxidizing_equivalent_rate_mol_m3_s, float),
    )
    if np.any(total_u <= 0.0):
        raise ValueError("total uranium concentration must be positive")
    if np.any((f3 <= 0.0) | (f3 >= 1.0)):
        raise ValueError("initial U(III) fraction must lie strictly between zero and one")
    u3_0 = total_u * f3
    u4_0 = total_u - u3_0
    transfer = rate * float(time_s)
    transfer = np.minimum(np.maximum(transfer, -0.999999 * u4_0), 0.999999 * u3_0)
    u3 = u3_0 - transfer
    u4 = u4_0 + transfer
    ratio0 = u4_0 / u3_0
    ratio = u4 / u3
    delta_E = R_GAS * T / FARADAY * np.log(ratio / ratio0)
    return {
        "u3_initial_mol_m3": u3_0,
        "u4_initial_mol_m3": u4_0,
        "u3_final_mol_m3": u3,
        "u4_final_mol_m3": u4,
        "u3_fraction_final": u3 / total_u,
        "delta_E_V": delta_E,
        "transfer_mol_m3": transfer,
        "fraction_u3_consumed": np.maximum(transfer, 0.0) / u3_0,
        "fraction_u4_reduced": np.maximum(-transfer, 0.0) / u4_0,
    }


def irradiation_time_s(
    burnup_GWd_tU: float,
    uranium_mass_kg: float,
    thermal_power_W: float,
) -> float:
    """Convert GWd/tU and capsule uranium mass to an equivalent full-power time."""
    if burnup_GWd_tU < 0.0 or uranium_mass_kg <= 0.0 or thermal_power_W <= 0.0:
        raise ValueError("burnup must be nonnegative and mass/power positive")
    uranium_mass_t = uranium_mass_kg / 1000.0
    return burnup_GWd_tU * GW_DAY_J * uranium_mass_t / thermal_power_W


def gas_capacitance_mol_Pa(
    salt_mesh: SaltMesh,
    plenum_mesh: SaltMesh,
    henry_mol_m3_Pa: float,
) -> float:
    """Inventory per pressure for nonuniform gas temperature and liquid partition."""
    gas_temperature = plenum_mesh.fields["Temperature"]
    gas_capacity = float(np.sum(plenum_mesh.cell_volume_m3 / (R_GAS * gas_temperature)))
    liquid_capacity = henry_mol_m3_Pa * salt_mesh.total_volume_m3
    return gas_capacity + liquid_capacity


def volume_average(mesh: SaltMesh, values: np.ndarray) -> float:
    return mesh.volume_average(np.asarray(values, dtype=float))


def volume_quantile(mesh: SaltMesh, values: np.ndarray, quantile: float) -> float:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    cumulative = np.cumsum(mesh.cell_volume_m3[order])
    target = quantile * cumulative[-1]
    return float(values[order[min(np.searchsorted(cumulative, target), len(values) - 1)]])


def volume_correlation(mesh: SaltMesh, a: np.ndarray, b: np.ndarray) -> float:
    """Return a cell-volume-weighted Pearson correlation coefficient."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.shape != mesh.cell_volume_m3.shape:
        raise ValueError("fields must match the cell-volume array")
    weights = mesh.cell_volume_m3 / mesh.total_volume_m3
    a_mean = float(np.sum(weights * a))
    b_mean = float(np.sum(weights * b))
    covariance = float(np.sum(weights * (a - a_mean) * (b - b_mean)))
    a_var = float(np.sum(weights * (a - a_mean) ** 2))
    b_var = float(np.sum(weights * (b - b_mean) ** 2))
    if a_var <= 0.0 or b_var <= 0.0:
        return math.nan
    return covariance / math.sqrt(a_var * b_var)


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _save_figure(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_base.with_suffix(".png")
    pdf_path = output_base.with_suffix(".pdf")
    fig.savefig(png_path, dpi=360, bbox_inches="tight", pad_inches=0.04)
    # Polygon-rich capsule maps are deliberately rasterized for compact, fast,
    # and reproducible PDF inclusion.  Line/contour figures remain vector PDF.
    if output_base.name.startswith("fig_capsule_"):
        from PIL import Image

        with Image.open(png_path) as image:
            image.convert("RGB").save(pdf_path, "PDF", resolution=360.0)
    else:
        fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)



__all__ = [name for name in globals() if not name.startswith("__")]
