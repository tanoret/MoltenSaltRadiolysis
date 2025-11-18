
from typing import Dict, List, Optional, Tuple
import numpy as np

def integrate_system(system, t_final: float, n_steps: int = 200, method: str = "auto"):
    """
    Integrate the ODEs for the system over [0, t_final] with n_steps outputs.
    Prefers SciPy's stiff solvers, but falls back to a simple implicit Euler if SciPy is unavailable.
    Returns (t, C, extras_time_series) where:
      - t is array [n_steps]
      - C is array [n_steps, n_species]
      - extras_time_series is list of dicts (e.g., gas pressures)
    """
    ts = np.linspace(0.0, t_final, n_steps)
    y0 = getattr(system, "initial_concentrations", None)
    if y0 is None:
        y0 = np.zeros(len(system.species))

    # try SciPy
    if method == "auto":
        try:
            from scipy.integrate import solve_ivp
            def fun(t, y):
                dydt, extra = system.rhs(t, y)
                return dydt
            sol = solve_ivp(fun, (ts[0], ts[-1]), y0, method="BDF", t_eval=ts, vectorized=False, atol=1e-12, rtol=1e-7)
            # collect extras by evaluating RHS on solution points
            extras = []
            for i in range(len(ts)):
                _, ex = system.rhs(ts[i], sol.y[:,i])
                extras.append(ex)
            return sol.t, sol.y.T, extras
        except Exception as e:
            # fallback
            pass

    # Fallback: implicit Euler with fixed substeps (not fancy but stable for demo)
    y = y0.copy()
    Ys = [y.copy()]
    extras = []
    for k in range(1, len(ts)):
        dt = ts[k] - ts[k-1]
        y = _implicit_euler(system, ts[k-1], y, dt, newton_iters=12, tol=1e-10)
        Ys.append(y.copy())
        _, ex = system.rhs(ts[k], y)
        extras.append(ex)
    extras.insert(0, extras[0] if extras else {})
    return ts, np.vstack(Ys), extras

def _implicit_euler(system, t, y, dt, newton_iters=8, tol=1e-9):
    import numpy as np
    yk = y.copy()
    for it in range(newton_iters):
        f, _ = system.rhs(t+dt, yk)
        R = yk - y - dt * f  # residual
        # numerical Jacobian
        n = len(y)
        J = np.zeros((n,n))
        eps = 1e-8
        for j in range(n):
            ytmp = yk.copy()
            dy = eps*(1.0 + abs(yk[j]))
            ytmp[j] += dy
            f1, _ = system.rhs(t+dt, ytmp)
            J[:, j] = ( (ytmp - y) - dt*f1 - R ) / dy  # derivative of residual wrt y_j
        # Solve J * delta = -R
        try:
            delta = np.linalg.solve(J, -R)
        except np.linalg.LinAlgError:
            delta = np.linalg.lstsq(J, -R, rcond=None)[0]
        yk = yk + delta
        if np.linalg.norm(delta, ord=np.inf) < tol:
            break
        # enforce non-negativity for liquid species; gas mol can be >=0
        # (weak projection helps stability)
        # (We assume indices of gas species will be identified by phase if needed; omitted here for simplicity)
    # clamp tiny negatives
    yk[yk < 0] = 0.0
    return yk
