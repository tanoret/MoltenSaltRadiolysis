
import json, argparse, pathlib
import numpy as np
from .builder import build_system
from .integrator import integrate_system
from .plotting import quick_plot, plot_gas

def main():
    parser = argparse.ArgumentParser(description="Run molten-salt radiolysis ODE model from a JSON config")
    parser.add_argument("config", type=str, help="Path to JSON config file")
    parser.add_argument("--tfinal", type=float, default=10.0, help="Final time (s)")
    parser.add_argument("--nsteps", type=int, default=200, help="Number of output steps")
    args = parser.parse_args()

    cfg_path = pathlib.Path(args.config)
    with open(cfg_path, "r") as f:
        config = json.load(f)

    system = build_system(config)
    t, C, extra = integrate_system(system, t_final=args.tfinal, n_steps=args.nsteps)

    # write CSV
    out_csv = cfg_path.with_suffix(".results.csv")
    header = ["t(s)"] + [s.name for s in system.species]
    arr = np.column_stack([t, C])
    np.savetxt(out_csv, arr, delimiter=",", header=",".join(header), comments="")
    print(f"Wrote {out_csv}")

    # plots
    fig1 = quick_plot(t, C, system)
    fig1_path = cfg_path.with_suffix(".species.png")
    fig1.savefig(fig1_path)
    print(f"Wrote {fig1_path}")

    # gas plot if present
    for g in ["Cl2","F2"]:
        fig = plot_gas(t, C, system, gas_name=g)
        if fig is not None:
            p = cfg_path.with_suffix(f".{g}_gas.png")
            fig.savefig(p)
            print(f"Wrote {p}")

if __name__ == "__main__":
    main()
