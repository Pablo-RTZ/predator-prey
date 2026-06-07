# Documentation

## Introduction

This repository is meant as an introduction into dynamical systems modelling, as well as its computational simulation. It will focus on ecological modelling, and mainly on numerical ODE solvers for the simulations. The aim is to provide both explanations and comparisons of basic ODE solvers (Euler, RK4,...), as well as some basic models for predator-prey relationships.

## The Lotka-Volterra model

The simplest model for predator-prey relationships is the Lotka-Volterra model. It assumes two species with populations $x_1,x_2$ respectively, where $x_1$ is the prey and $x_2$ the predator. The evolution of the populations is given by the following ODE.

$$
\begin{matrix}
\frac{dx_1}{dt}=a x_1-bx_1x_2 \\
\frac{dx_2}{dt}= c x_1x_2-dx_2
\end{matrix}
$$

where

- $a$: Prey reproduction rate
- $b$: Depredation rate
- $c$: Predator reproduction rate
- $d$: Predator death rate

For any set of initial conditions (real and positive), the model will generally show out of sync oscillations, with predator populations peaking after the prey populations. Despite not having a closed form solution, the behavior is generally stable, and is non-stiff for moderate values for the parameters. Therefore, it will be used as a benchmark to compare the various numerical solvers used in further sections.

The `LotkaVolterra.mlx` notebook inside the Notebooks folder shows some simulations using various methods, and parameters can be tweaked to see how the model reacts. By taking the following assumptions,

- Prey grows by 50% monthly in absence of predators
- 0.02 Prey encounters per month per predator
- 1% Predator efficiency (in converting prey into offspring)
- Predator population declines by 40% monthly in absence of prey
- Initial population of 40 prey and 10 predators

the resulting simulation (using RK4, 1000 timesteps) is

![Lotka-Volterra model](../Assets/LotkaVolterra.png)

### Lotka-Volterra equilibrium points

By setting the time derivatives to 0, equilibrium points for the model can be found. In this case,

$$
\begin{rcases}
0 =ax_1-bx_1x_2\\0=cx_1x_2-dx_2
\end{rcases}
\implies \begin{cases}
x_1=x_2=0\\ x_1=\frac dc,\quad x_2=\frac ab
\end{cases}
$$

It is trivial to verify that for those pairs of initial conditions, the system will remain at equilibrium.

## ODE solvers comparison over Lotka-Volterra

Due to this model not having a closed form solution, self-convergence tests will be used. In this case, Richardson convergence analysis will be used, where succesively precise iterations are compared, in order to obtain an aproximation of the theoretical convergence order. The aproximate convergence order is given by

$$
p= \frac{\log \left(\frac{\|x_h - x_{h/2}\|}{\|x_{h/2} - x_{h/4}\|}\right)}{\log 2}
$$

A function is defined (in the `NumericalComparison.mlx`notebook) to, given $N$ and $iter$, calculate succesive iterations using $N,2N,4N,...$ steps, and get $iter$ different aproximations for the convergence orders.

Taking 10 initial timesteps, and doubling each iteration, the convergence orders are the following.

| Method         | Value 1 | Value 2 | Value 3 | Value 4 | Value 5 |
|----------------|---------|---------|---------|---------|---------|
| Euler          | 1.5680  | 1.2603  | 1.1420  | 1.0729  | 1.0378  |
| AdamsBashfort  | 1.9568  | 2.0980  | 1.9917  | 1.9872  | 1.9907  |
| RK4            | 3.9980  | 3.9791  | 3.9678  | 3.9832  | 3.9915  |
| DormandPrince  | 6.5154  | 5.8117  | 5.5196  | 5.3225  | 5.1860  |

The main limitation of this testing method is that, if variable precission arithmetic is not used, the higher order methods' calculations will eventually collapse, due to rounding error dominating over the method's error.

Compute time can also be compared for these methods. Despite the higher order methods having higher compute times (~5x at worst case), the massive precision increase is almost always worth it, due to the error being various orders of magnitude smaller. The results, for the implementations in this repository are the following.

![Computation time comparison](../Assets/ComputeTime.png)

## Logistic prey growth

As it has been found on the previous section, the compute cost of higher order methods is justifiable, so from now on, most simulation will onle be run on Runge-Kutta (RK4) or Dormand-Prince (DOPRIS5). The next simplest modification that can be done to the Lotka-Volterra model is to model prey evolution (in absence of predators), using logistic growth. This way, a carrying capacity $K$ can for the system can be added, so as prey population approaches it, its growth rate slows down, to simulate the decreasing resource availability. Once this aspect is considered, the resulting ODE is as follows

$$
\begin{matrix}
\frac{dx_1}{dt}=a x_1(1-\frac{x_1}{K})-bx_1x_2 \\
\frac{dx_2}{dt}= c x_1x_2-dx_2
\end{matrix}
$$

where

- $a$: Prey reproduction rate
- $b$: Depredation rate
- $c$: Predator reproduction rate
- $d$: Predator death rate
- $K$: Carrying capacity

### Logistic model equilibrium points

By setting the time derivatives to 0, the equilibrium points can be found.

$$
\begin{rcases}
0 =ax_1(1-\frac{x_1}K)-bx_1x_2\\0=cx_1x_2-dx_2
\end{rcases}
\implies \begin{cases}
x_1=x_2=0\\ x_1=\frac dc,\quad x_2=\frac ab(1-\frac{d}{cK})
\end{cases}
$$

### Numerical simulation

By using the same parameters as before, but adding a carrying capacity, the model behavior changes, and now shows damped oscillations around the equilibrium points for the system. This shows how reducing growth rate when population approaches the carrying capacity adds a damping factor that stabilizes behavior. Simulating with Runge-Kutta, the following results are found.

![Logistic prey growth](../Assets/Logistic.png)

## Predator functional response

Once prey growth rate has been limited, the next step to improve the current model, is to limit predator growth. Right now, sufficiently big prey populations would increase predator growth rate regardless of predator population. The equivalent to logistic growth, but applied to predator growth rate would be to add a functional response, so predator growth rate is capped by predator population. This makes ecological sense, because a point would be reached where increasing prey population would not further increase predator growth rate, as predators would be saturated (limited by reproduction or digestion rate).
