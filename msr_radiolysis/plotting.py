
import matplotlib.pyplot as plt
import numpy as np

def quick_plot(t, C, system, species_to_plot=None, title="Radiolysis species (liq)"):
    if species_to_plot is None:
        # default: plot first 6 liquid species
        species_to_plot = [s.name for s in system.species if s.phase=="liq"][:6]
    plt.figure()
    for s in species_to_plot:
        if s not in system.species_index:
            continue
        idx = system.species_index[s]
        plt.plot(t, C[:, idx], label=s)
    plt.xlabel("Time (s)")
    plt.ylabel("Concentration (mol/m^3)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    return plt.gcf()

def plot_gas(t, C, system, gas_name="Cl2"):
    # plots gas moles and pressure for target gas
    liq = f"{gas_name}_diss"
    gas = f"{gas_name}_g"
    if gas not in system.species_index:
        return None
    idx_g = system.species_index[gas]
    n = C[:, idx_g]
    p = (n * system.R_gas * system.T) / max(system.V_gas_m3, 1e-30) if system.V_gas_m3>0 else np.zeros_like(n)
    plt.figure()
    plt.plot(t, n, label=f"{gas} mol (gas)")
    plt.xlabel("Time (s)")
    plt.ylabel("Moles (mol)")
    plt.title(f"{gas} gas accumulation")
    plt.legend()
    plt.tight_layout()
    return plt.gcf()
