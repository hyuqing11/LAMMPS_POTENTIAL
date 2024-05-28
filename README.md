# Lennard-Jones Potential for Small Displacements

The Lennard-Jones (LJ) potential describes the interaction between particles in a system. For small displacements about the equilibrium position $r_0 = 2^\frac{1}{6} \sigma$, the potential energy $\Delta u(r)$ is given by:

$$\Delta u(r) = \frac{1}{2} C (r - r_0)^2 + \frac{1}{3} \alpha (r - r_0)^3 + \frac{1}{4} \beta (r - r_0)^4$$

Where:
- $C = 36 \sqrt[3]{4} \epsilon \sigma^{-3}$
- $\alpha = -378 \sqrt{2} \sigma \epsilon^{-3}$
- $\beta = 2226 \sqrt[3]{2} \sigma^{-4}$

To utilize this potential in LAMMPS, follow these steps:

1. Add the following two files to the lammps/src directory:
   - pair_quart_cut.cpp
   - pair_quart_cut.h

2. Recompile LAMMPS.

In the input file for LAMMPS, use the following commands:

- `pair_style quart/cut 1.5`: This command sets the pair style to quartic potential with a cutoff distance of 1.5 units. It models the interaction between particles in the system.

- `pair_coeff * * 1.0 1.0 -1.0 0.0 1.5`: This command sets the coefficients for the pair interactions. Here's what each parameter represents:
  - The first two asterisks (*) indicate that these coefficients apply to all pairs of atoms in the system.
  - The first two numbers (1.0 1.0) represent the energy and length scale parameters, $\epsilon$ and $\sigma$, respectively.
  - The third number (-1.0) modifies the potential to eliminate discontinuities in both $u(r)$ and $-u'(r)$.
  - The fourth and fifth numbers (0.0 1.0) control whether to open or close the second and third terms in the LJ potential equation, respectively. In this case, In this case, the second term is excluded (0.0) and the third term is included (1.0).


