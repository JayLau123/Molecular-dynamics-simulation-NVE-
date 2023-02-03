# Molecular-dynamics-simulation-NVE
Tutorial: Molecular dynamics simulation of Lennard-Jones particles in microcanonical ensemble.

Please find more detailed comments in jupyter notebook.

Question: how to choose the best time step $\Delta t$ in MD? A brief Numerical analysis of $\Delta t$.


## Overview

Define a system as a set of N particles
Follow the evolution of this system by integrating the equations of motion
For particle i, $F_i=m_ia_i$, note that $F_i$ depends on all other particles, so that for the entire system we have 3N coupled equations.

The equation of motion of each particle is integrated numerically. We increase the time in small $\Delta t$ step by step, and determine the positions of all particles at the discrete times.

This implies the following iteration:

Use positions at current time to calculate all the forces, $x(t) \Rightarrow F$

Forces yields the acceleration of each particle, $F \Rightarrow a$

The positions at the new time are computed from the current positions and the current velocities, $x(t),v(t) \Rightarrow x(t+\Delta t)$

The new velocities are computed from the current velocities and the current accelerations, $v(t),a(t) \Rightarrow v(t+\Delta t)$

### Verlet integrator

One of the most widely used "integrators" (schemes to compute the new positions and velocities from those at an earlier time) is the so-called Verlet algorithm. It follows simply from a Taylor expansion of the positions $x(t)$ as a function of time $t$


$x(t+\Delta t)=x(t)+\dfrac{dx(t)}{dt}\Delta t+\dfrac{1}{2}\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}+\dfrac{1}{6}\dfrac{d^{3}x(t)}{dt^{3}}(\Delta t)^{3}+...$

$x(t-\Delta t)=x(t)-\dfrac{dx(t)}{dt}\Delta t+\dfrac{1}{2}\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}-\dfrac{1}{6}\dfrac{d^{3}x(t)}{dt^{3}}(\Delta t)^{3}+...$

Adding these two equations yields:

$x(t+\Delta t)=2x(t)-x(t-\Delta t)+\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}+\mathcal{O}((\Delta t)^4)$

or

$x(t+\Delta t)=2x(t)-x(t-\Delta t)+a(t)(\Delta t)^{2}+\mathcal{O}((\Delta t)^4)$

The velocity $v(t)$ is abtained via:

$v(t)=\dfrac{x(t+\Delta t)-x(t-\Delta t)}{2\Delta t}-\frac{x^{(3)}(\xi)}{6}(\Delta t)^{3}$






















