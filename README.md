# README

[![Matlab Version](https://img.shields.io/badge/MATLAB-R2024%2B-blue)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

A repository on the modelling of predator-prey relationships. It aims to showcase the use of numerical ODE solvers to simulate predator-prey relationships in an ecological context.

## Quick start / Setup

In order to use this programs, you have to

1. Clone or download the repository
2. Add it to your MATLAB path
3. Verify installation (Run `startup` and open any of the notebooks)
4. Head over to the [documentation](Docs/Documentation.md) to see the writeup, or to the [index](Docs/INDEX.md)

The folder contains a `startup.m` file, so running `startup` in the terminal will add all subfolders to your working directory until the program is closed. This way, programs from any folder can be run from the MATLAB folder inside the repository root directory.

### Prerequisites

- MatLab R2024b or later

The code has been developed in MatLab R2024b, and requires no additional toolboxes or dependencies.

### Example usage

On top of tweaking parameters on the [notebooks](MATLAB/Notebooks/), running the models can be as simple as

```matlab
a = 0.1; b = 0.25; c = 0.2; d = 0.2;
prey = 1; predator = 0.2; n = 1000; t_end = 400;

f = LotkaVolterra(a,b,c,d);

[t,x] = Euler(f, 0, t_end, [prey; predator], n);

plot(t,x)
legend("Prey","Predator")
```

## Contents

The repository aims to cover the following topics

- Classic deterministic model (Lotka-Volterra)
- Numerical comparison of ODE solvers
- Logistic growth for prey
- Functional responses for predation
- Seasonality
- Spatial models (reaction-diffusion)

## License

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

See the [LICENSE](LICENSE.md) file for details.

If you use this repository, please cite following the [citation](citation.cff).
