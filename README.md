# 2DIsingMC.jl
This small Julia script runs a Monte Carlo simulation of the two dimensional Ising model in absence of magnetic field. It uses the Metropolis algorithm. It provides a file with the values of the Energy and the Magnetization for the considered lattice size at each of the Monte Carlo steps per site of the simulation. The original purpose of this script is to provide raw data for the calculation of the critical properties of the Ising model, as a test for the implementation of a finite size scaling library written in Python.

It is also a first try on using Julia for something more or less useful.

# Usage

1. Make sure that Julia>=1.8.3 is installed on your computer

2. Clone the repository

```bash
git clone https://github.com/cobos99/2DIsingMC.jl.git
```

3. Change directory into the repository's folder

```bash
cd 2DIsingMC.jl
```

4. Provide the required options to run the simulation as arguments in the program call

```bash
julia run2disingmc.jl [-i INIT] [-s SEACH] L temp J exec_time
```

### Required arguments:
    
- L : Number of nodes on the side of the considered lattice

- temp : Temperature of the system (in units of $k_\mathrm{B} = 1$)

- J : Coupling constant. The energy of the system is $E = -J\sum_{\left\langle i,j \right\rangle} S_i S_j$ and its critical temperature $T_c = 2J/\ln (1 + \sqrt{2})$ [https://en.wikipedia.org/wiki/Ising_model#Onsager's_exact_solution]

- exec_time : Montecarlo steps per site $t_s$ for the simulation. This refers to a total of $L^2t_s$ iterations of the Metropolis algorithm

### Optional arguments:

- INIT [-i] : Initial state. Accepts either of the strings ['r', 'f', 'fn']. Each meaning that either a random or ferromagnetic initial state will be considered. It also accepts the path to a file containing a square table of $\pm 1$ values refering to the initial state of each of the spins.  Defaults to 'r'.

- SEACH [-s] : Change the frequency at which the energy and magnetization of the system are saved to the output file. Defaults to saving each step per site.
