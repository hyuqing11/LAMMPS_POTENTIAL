**The is the potential for Lennard-Johns potential for small displacements about the equalibrium position $r_0 = 2^\frac{1}{6} \sigma$:**
$$ \Delta U_{LJ} = \frac{1}{2} C (r-r_0)^2 + \frac{1}{3}\alpha (r-r_0)^3 + \frac{1}{4} \beta (r-r_0)^4 $$
Where $`C = 36\sqrt[3]{4}\epsilon\sigma^{-3}`$. The $\alpha$ and $\beta$ satisfy the Taylor approximation of the LJ interaction and they take the constant values $` \alpha = -378 \sqrt{2}\sigma \epsilon^{-3} `$ and $\beta = 2226 \sqrt[3]{2}\sigma^{-4}$
In order to use this potential, you need to add these two file to the lammps/src and then recompile lammps
Then in the inputfile for lammps using the command:
pair_style      quart/cut 1.5
pair_coeff      * * 1.0 1.0 -1.0 0.0 1.0
in the pair_coeff, the first term is $epsilon$; The second term is $sigma$; "1" and "0" used to whether open open the second term and third term in the equation

**The Cauchy-Schwarz Inequality**
$$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$
