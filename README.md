# Train Suspension Spring-Mass System — Numerical Simulation

Numerical simulation of a two-mass, two-spring system modeling train suspension dynamics. Compares three ODE integrators — **Euler**, **Heun**, and **4th-order Runge-Kutta** — on accuracy and approximate error. Originally written in MATLAB; ported to Python with an interactive Streamlit/Plotly demo.

<p>
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white">
  <img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-original-orange?logo=mathworks&logoColor=white">
  <img alt="Streamlit" src="https://img.shields.io/badge/Demo-Streamlit-FF4B4B?logo=streamlit&logoColor=white">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-yellow.svg">
</p>

**[Live demo →](#deployment)**

## Overview

Two masses (`m1`, `m2`) connected in series by two springs (`k1`, `k2`) model a simplified train suspension. The system is governed by:

```
dx1/dt = v1
dx2/dt = v2
dv1/dt = -(k1/m1)(x1 - L1) + (k2/m1)(x2 - x1 - w1 - L2)
dv2/dt = -(k2/m2)(x2 - x1 - w1 - L2)
```

Each method integrates this coupled ODE system over a configurable time horizon and step size, and the approximate percent error `Ea = |x(i+1) - x(i)| / x(i+1) * 100` is computed at every step to compare convergence behavior.

## Repository Contents

| File | Description |
|---|---|
| `Final_matlab_code.m` | Original MATLAB implementation (interactive `input()` prompts, Euler/Heun/RK4, plots) |
| `Applied Numerical Methods project (Spring Mass System).pdf` | Written project report |
| `spring_mass_sim.py` | Python/NumPy port of the simulation core — `SimParams` dataclass + `simulate_euler` / `simulate_heun` / `simulate_rk4` / `run_all` |
| `spring_app.py` | Interactive Streamlit demo with live parameter sliders and Plotly charts |
| `requirements.txt` | Python dependencies for the demo |

## Features

- **Three integrators side by side** — Euler (1st order), Heun (2nd order predictor-corrector), and classic RK4, run on identical parameters for direct comparison.
- **Interactive controls** — tune spring constants, masses, unstretched lengths, step size, and simulation time from the sidebar; results update on demand.
- **Displacement, velocity, and error views** — tabbed Plotly charts for `x1`/`x2`, `v1`/`v2`, plus a table of final-step approximate errors per method.
- **Faithful port** — the Python core mirrors the MATLAB derivative functions and update equations line-for-line, so results match the original within floating-point precision.

## How to Run

```bash
# Clone
git clone https://github.com/Abdelrahmanemad0/Spring-Mass-System.git
cd Spring-Mass-System

# Install dependencies
pip install -r requirements.txt

# Launch the interactive demo
streamlit run spring_app.py

# Or run the core simulation from the command line
python spring_mass_sim.py
```

The original MATLAB script (`Final_matlab_code.m`) still runs as-is in MATLAB/Octave if you'd rather use the reference implementation — it prompts for each parameter interactively and produces the same displacement/velocity/error plots via `figure`/`plot`.

## Deployment

Deployable in minutes on **Streamlit Community Cloud**: fork this repo, connect it at [share.streamlit.io](https://share.streamlit.io), and set `spring_app.py` as the entry point. No secrets required.

## Tech Stack

Python, NumPy, Pandas, Plotly, Streamlit · MATLAB (original reference implementation)

## License

MIT — see [LICENSE](LICENSE).

#numericalmethods #matlab #python #simulation #eulermethod #rungekutta #engineering
