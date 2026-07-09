"""spring_mass_sim.py -- Core numerical simulation for the two-mass spring
system (train suspension model). Pure NumPy port of Final_matlab_code.m,
solving the same coupled ODEs with three integrators (Euler, Heun, RK4) so
their accuracy/error can be compared side by side.

System: two masses connected in series by two springs.
    m1 sits at rest length L1 from a fixed wall (spring k1).
    m2 sits beyond m1, connected to it by spring k2 (rest length L2).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SimParams:
    k1: float = 5.0    # Spring 1 constant (N/m)
    k2: float = 5.0    # Spring 2 constant (N/m)
    m1: float = 2.0    # Mass 1 (kg)
    m2: float = 2.0    # Mass 2 (kg)
    L1: float = 2.0    # Unstretched length of spring 1 (m)
    L2: float = 2.0    # Unstretched length of spring 2 (m)
    w1: float = 5.0    # Width of mass 1 (m)
    w2: float = 5.0    # Width of mass 2 (m) -- kept for parity with the MATLAB inputs
    h: float = 0.01    # Step size (s)
    t_final: float = 20.0  # Simulation horizon (s)


def _derivatives(params: SimParams):
    """Return the four derivative functions dx1/dt, dx2/dt, dv1/dt, dv2/dt."""
    k1, k2, m1, m2, L1, L2, w1 = (
        params.k1, params.k2, params.m1, params.m2,
        params.L1, params.L2, params.w1,
    )

    def dx1dt(v1):
        return v1

    def dx2dt(v2):
        return v2

    def dv1dt(x1, x2):
        return (-k1 / m1) * (x1 - L1) + (k2 / m1) * (x2 - x1 - w1 - L2)

    def dv2dt(x1, x2):
        return (-k2 / m2) * (x2 - x1 - w1 - L2)

    return dx1dt, dx2dt, dv1dt, dv2dt


def _initial_state(params: SimParams):
    x1_0 = params.L1
    x2_0 = params.L1 + params.w1 + params.L2 + 6
    return x1_0, x2_0, 0.0, 0.0


def _time_grid(params: SimParams) -> np.ndarray:
    n_steps = int(round(params.t_final / params.h)) + 1
    return np.linspace(0.0, params.t_final, n_steps)


def _approx_error(prev: np.ndarray, curr: np.ndarray) -> np.ndarray:
    """Percent approximate relative error between consecutive steps."""
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.abs((curr - prev) / curr) * 100
    return np.nan_to_num(err, nan=0.0, posinf=0.0, neginf=0.0)


def simulate_euler(params: SimParams) -> dict:
    """Explicit (forward) Euler integration."""
    dx1dt, dx2dt, dv1dt, dv2dt = _derivatives(params)
    t = _time_grid(params)
    n = len(t)
    x1, x2, v1, v2 = (np.zeros(n) for _ in range(4))
    x1[0], x2[0], v1[0], v2[0] = _initial_state(params)

    h = params.h
    for i in range(n - 1):
        x1[i + 1] = x1[i] + h * dx1dt(v1[i])
        x2[i + 1] = x2[i] + h * dx2dt(v2[i])
        v1[i + 1] = v1[i] + h * dv1dt(x1[i], x2[i])
        v2[i + 1] = v2[i] + h * dv2dt(x1[i], x2[i])

    return {"t": t, "x1": x1, "x2": x2, "v1": v1, "v2": v2}


def simulate_heun(params: SimParams) -> dict:
    """Heun's method (predictor-corrector, 2nd order)."""
    dx1dt, dx2dt, dv1dt, dv2dt = _derivatives(params)
    t = _time_grid(params)
    n = len(t)
    x1, x2, v1, v2 = (np.zeros(n) for _ in range(4))
    x1[0], x2[0], v1[0], v2[0] = _initial_state(params)

    h = params.h
    for i in range(n - 1):
        x1_pred = x1[i] + h * dx1dt(v1[i])
        x2_pred = x2[i] + h * dx2dt(v2[i])
        v1_pred = v1[i] + h * dv1dt(x1[i], x2[i])
        v2_pred = v2[i] + h * dv2dt(x1[i], x2[i])

        x1[i + 1] = x1[i] + (h / 2) * (dx1dt(v1[i]) + dx1dt(v1_pred))
        x2[i + 1] = x2[i] + (h / 2) * (dx2dt(v2[i]) + dx2dt(v2_pred))
        v1[i + 1] = v1[i] + (h / 2) * (dv1dt(x1[i], x2[i]) + dv1dt(x1_pred, x2_pred))
        v2[i + 1] = v2[i] + (h / 2) * (dv2dt(x1[i], x2[i]) + dv2dt(x1_pred, x2_pred))

    return {"t": t, "x1": x1, "x2": x2, "v1": v1, "v2": v2}


def simulate_rk4(params: SimParams) -> dict:
    """Classic 4th-order Runge-Kutta."""
    dx1dt, dx2dt, dv1dt, dv2dt = _derivatives(params)
    t = _time_grid(params)
    n = len(t)
    x1, x2, v1, v2 = (np.zeros(n) for _ in range(4))
    x1[0], x2[0], v1[0], v2[0] = _initial_state(params)

    h = params.h
    for i in range(n - 1):
        k1x1 = h * dx1dt(v1[i])
        k1x2 = h * dx2dt(v2[i])
        k1v1 = h * dv1dt(x1[i], x2[i])
        k1v2 = h * dv2dt(x1[i], x2[i])

        k2x1 = h * dx1dt(v1[i] + k1v1 / 2)
        k2x2 = h * dx2dt(v2[i] + k1v2 / 2)
        k2v1 = h * dv1dt(x1[i] + k1x1 / 2, x2[i] + k1x2 / 2)
        k2v2 = h * dv2dt(x1[i] + k1x1 / 2, x2[i] + k1x2 / 2)

        k3x1 = h * dx1dt(v1[i] + k2v1 / 2)
        k3x2 = h * dx2dt(v2[i] + k2v2 / 2)
        k3v1 = h * dv1dt(x1[i] + k2x1 / 2, x2[i] + k2x2 / 2)
        k3v2 = h * dv2dt(x1[i] + k2x1 / 2, x2[i] + k2x2 / 2)

        k4x1 = h * dx1dt(v1[i] + k3v1)
        k4x2 = h * dx2dt(v2[i] + k3v2)
        k4v1 = h * dv1dt(x1[i] + k3x1, x2[i] + k3x2)
        k4v2 = h * dv2dt(x1[i] + k3x1, x2[i] + k3x2)

        x1[i + 1] = x1[i] + (k1x1 + 2 * k2x1 + 2 * k3x1 + k4x1) / 6
        x2[i + 1] = x2[i] + (k1x2 + 2 * k2x2 + 2 * k3x2 + k4x2) / 6
        v1[i + 1] = v1[i] + (k1v1 + 2 * k2v1 + 2 * k3v1 + k4v1) / 6
        v2[i + 1] = v2[i] + (k1v2 + 2 * k2v2 + 2 * k3v2 + k4v2) / 6

    return {"t": t, "x1": x1, "x2": x2, "v1": v1, "v2": v2}


def run_all(params: SimParams) -> dict:
    """Run Euler, Heun, and RK4 and attach final-step approximate errors."""
    results = {
        "euler": simulate_euler(params),
        "heun": simulate_heun(params),
        "rk4": simulate_rk4(params),
    }
    for name, res in results.items():
        for var in ("x1", "x2", "v1", "v2"):
            err = _approx_error(res[var][:-1], res[var][1:])
            res[f"err_{var}"] = err
    return results


if __name__ == "__main__":
    params = SimParams()
    out = run_all(params)
    print("Final displacements/velocities (RK4 reference):")
    print(f"  x1={out['rk4']['x1'][-1]:.4f} m  x2={out['rk4']['x2'][-1]:.4f} m")
    print(f"  v1={out['rk4']['v1'][-1]:.4f} m/s  v2={out['rk4']['v2'][-1]:.4f} m/s")
