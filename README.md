# Molecular-dynamics-simulation-NVE
Tutorial: Molecular dynamics simulation of Lennard-Jones particles in microcanonical ensemble.

Please find more detailed comments in jupyter notebook.

Question: how to choose the best time step $\Delta t$ in MD? A brief Numerical analysis of $\Delta t$.


## Overview

Define a system as a set of N particles
Follow the evolution of this system by integrating the equations of motion
For particle i, $F_i=m_ia_i$, note that $F_i$ depends on all other particles, so that for the entire system we have 3N coupled equations

The equation of motion of each particle is integrated numerically. We increase the time in small $\Delta t$ step by step, and determine the positions of all particles at the discrete times.

