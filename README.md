# Train Suspension Spring-Mass System — Numerical Simulation

This repository contains a MATLAB-based numerical simulation of a two-mass, two-spring system representing train suspension dynamics. It implements **Euler's Method**, **Heun's Method**, and the **4th-order Runge-Kutta (RK4) Method** to solve the system of ordinary differential equations (ODEs) governing the motion, and compares their accuracy and stability.

## Features

- **Mathematical Modeling** — formulates the two-mass spring system as a system of ODEs
- **Numerical Methods** — implements Euler, Heun, and RK4 for solving the equations
- **Error Analysis** — evaluates the approximate percent error and stability of each method
- **Visualization** — plots displacement and velocity of both masses over time

## Technologies Used

- MATLAB
- Numerical methods for ODEs (Euler, Heun, RK4)
- Data analysis and error estimation

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/Abdelrahmanemad0/Spring-Mass-System.git
   cd Spring-Mass-System
   ```
2. Open MATLAB and run `spring_mass_simulation.m`.
3. Enter the requested system parameters when prompted: spring constants (`k1`, `k2`), masses (`m1`, `m2`), unstretched spring lengths (`L1`, `L2`), mass widths (`w1`, `w2`), and step size (`h`). Reasonable reference values are listed in a comment near the top of the script.
4. Compare the printed final displacements/velocities and approximate errors, and inspect the generated plots, to see how the three methods compare.

## Project Structure

- `spring_mass_simulation.m` — main MATLAB script (renamed from `Final_matlab_code.m` to match these run instructions)
- `Applied Numerical Methods project (Spring Mass System).docx` — written report covering the modeling approach, methods, and results
- `README.md` — this file
- `LICENSE` — MIT license
- `.gitignore` — MATLAB build/editor artifacts to keep out of version control

## Fixes in this revision

- **Renamed `Final_matlab_code.m` to `spring_mass_simulation.m`** — the README always instructed users to run `spring_mass_simulation.m`, but the actual file in the repo was named `Final_matlab_code.m`, so following the instructions verbatim would fail.
- Added `LICENSE` and `.gitignore`.

## Future Enhancements

- Add automated unit tests comparing each method's output against a known analytical/reference solution
- Parameterize the script so it can also be run non-interactively (e.g., as a function with default arguments)
- Extend the model to n-mass chains

This project highlights the importance of accurate numerical methods in real-world engineering applications, particularly in train suspension systems.

## License

MIT — see [LICENSE](LICENSE).
