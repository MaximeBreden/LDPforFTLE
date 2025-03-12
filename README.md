# Large Deviation Principle for Finite Time Lyapunov Exponents.

This repository contains the Matlab codes used for the paper [*Detecting random bifurcations via rigorous enclosures of large
deviations rate functions*](https://arxiv.org/abs/2408.12556), by Alexandra Blessing (Neamtu), Alex Blumenthal, Maxime Breden and Maximilian Engel. The proofs make usage of
interval arithmetic via the [Intlab toolbox](http://www.ti3.tu-harburg.de/intlab/).

## Pitchfork bifurcation

The folder "Pitchfork" contains the codes related to the theorems discussed in Section 5.1 of the paper:
- The proof of Theorem 5.5 can be reproduced by running script_pitchfork_single.m (which requires Intlab);
- The proof illustrated in Figure 1 can be reproduced by running script_pitchfork_all.m (which requires Intlab);
- The proof of Theorem 5.7 can be reproduced by running script_pitchfork_minimum.m (which requires Intlab).
  
Two additional scripts containing non-rigorous computations are provided:
- script_pitchfork_numerics.m provides faster but not guaranteed approximations of the rate function at 0;
- script_pitchfork_asymptoticLE.m approximates the asymptotic Lyapunov exponent and can be used to reproduce Figure 2.

## Toy model of shear-induces chaos

The folder "Shear" contains the codes related to the theorems discussed in Section 5.2 of the paper:
- The proof of Theorem 5.14 can be reproduced by running script_shear_single.m (which will run whithout Intlab, but then the result is of course not guaranteed);
- The proof illustrated in Figure 4 can be reproduced by running script_shear_all.m (which will run whithout Intlab, but then the result is of course not guaranteed).
