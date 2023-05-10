# Molecular-dynamics-simulation-NVE
Tutorial: Molecular dynamics simulation of Lennard-Jones particles in microcanonical ensemble.

Please find more detailed comments in jupyter notebook.

Question: how to choose the best time step $\Delta t$ in MD? A brief Numerical analysis of $\Delta t$.


## Overview

Molecular dynamics simulations are an important tool for understanding dynamic processes and mechanisms on a microscopic level in various areas of chemistry, biology, and materials science. Besides the prediction of equili- brium properties of molecules and materials, they also offer the possibility to simulate excited states and non-equilibrium dynamics, as well as slow processes and rare events.

Define a system as a set of N particles
Follow the evolution of this system by integrating the equations of motion
For particle i, $F_i=m_ia_i$, note that $F_i$ depends on all other particles, so that for the entire system we have 3N coupled equations.
Here, we use the NVE ensemble, energy must be conserved.

The equation of motion of each particle is integrated numerically. We increase the time in small $\Delta t$ step by step, and determine the positions of all particles at the discrete times.

This implies the following iteration:

Use positions at current time to calculate all the forces, $x(t) \Rightarrow F(t)$

Forces yields the acceleration of each particle, $F(t) \Rightarrow a(t)$

The positions at the new time are computed from the current positions and the current velocities, $x(t),v(t) \Rightarrow x(t+\Delta t)$

The new velocities are computed from the current velocities and the current accelerations, $v(t),a(t) \Rightarrow v(t+\Delta t)$

## Verlet integrator

One of the most widely used "integrators" (schemes to compute the new positions and velocities from those at an earlier time) is the so-called Verlet algorithm. It follows simply from a Taylor expansion of the positions $x(t)$ as a function of time $t$


$$x(t+\Delta t)=x(t)+\dfrac{dx(t)}{dt}\Delta t+\dfrac{1}{2}\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}+\dfrac{1}{6}\dfrac{d^{3}x(t)}{dt^{3}}(\Delta t)^{3}+...$$

$$x(t-\Delta t)=x(t)-\dfrac{dx(t)}{dt}\Delta t+\dfrac{1}{2}\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}-\dfrac{1}{6}\dfrac{d^{3}x(t)}{dt^{3}}(\Delta t)^{3}+...$$

Adding these two equations yields:

$$x(t+\Delta t)=2x(t)-x(t-\Delta t)+\dfrac{d^{2}x(t)}{dt^{2}}(\Delta t)^{2}+\mathcal{O}((\Delta t)^4)$$

or

$$x(t+\Delta t)=2x(t)-x(t-\Delta t)+a(t)(\Delta t)^{2}+\mathcal{O}((\Delta t)^4)$$

The velocity $v(t)$ is abtained via:

$$v(t)=\dfrac{x(t+\Delta t)-x(t-\Delta t)}{2\Delta t}-\dfrac{x^{(3)}(\xi)}{6}(\Delta t)^{3}$$

or 

$$v(t)\approx \dfrac{x(t+\Delta t)-x(t-\Delta t)}{2\Delta t}$$

### Disadvantages: we have to store three sets of positions: $x(t-\Delta t), x(t), x(t+\Delta t)$

Here, we use the NVE ensemble, energy must be conserved!

Round-off error endanger this. In a "symmetric" algorithm, there will be no net drift of the energy.

This method is not "self-starting" ! Why not?

Use a backward Euler method method to get $x(-\Delta t)$

$$x(-\Delta t)=x(0)-v(0)\Delta t$$

Why $\Delta t$ matters? The calculation of the froces is very time consuming. The bigger steps we can make, the better. (There is an interesting numerical analysis of $\Delta t$, please see another file in the same directory)

Higher-order integrators, such as predictor-corrector methods, in principle permit bigger time steps.In practice, this is only true for really small $\Delta t$, for large steps, the simple methods win.

## Leap-frog algorithm

Rewrite the Taylor expansion for the position:

$$x(t+\Delta t)=x(t)+\Delta t\left[\dfrac{d x}{d t}+\dfrac{1}{2} \dfrac{d^2 x}{d t^2} \Delta t\right]+\ldots$$

The term in square brackets can be interpreted as an approximation for $v(t+\Delta t/2)$,so

$$x(t+\Delta t)=x(t)+v(t+\Delta t / 2) \Delta t$$

Meanwhile, the velocity at $t+\Delta t/2$ is obtained from:

$$v(t+\Delta t / 2)=v(t-\Delta t / 2)+a(t) \Delta t$$

### Advantages of Leap-frog algorithm

No worries because of inaccuracy of higher powers of \Delta t.

No need to store three sets of positions.

## Velocity Verlet

Write:

$$v\left(t+\dfrac{\Delta t}{2}\right)=v(t)+\dfrac{1}{2} a(t) \Delta t$$

$$v(t+\Delta t)=v\left(t+\dfrac{\Delta t}{2}\right)+\dfrac{1}{2} a(t+\Delta t) \Delta t$$

### note: we use $a(t+\Delta t)$ instead of $a(t+\Delta t/2)$ in above equation. When the $\Delta t$ is very small, $\Delta t/2$ is very small, $a(t+\Delta t) \approx a(t+\Delta t/2)$.

$$\Rightarrow v(t+\Delta t)=v(t)+\dfrac{1}{2} a(t) \Delta t+\dfrac{1}{2} a(t+\Delta t) \Delta t $$

$$\Rightarrow v(t+\Delta t)=v(t)+\frac{1}{2}[a(t+\Delta t)+a(t)] \Delta t$$


Hence,

$$v(t+\Delta t)=v(t)+\frac{1}{2}[a(t+\Delta t)+a(t)] \Delta t$$

$$x(t+\Delta t)=x(t)+v(t) \Delta t+\dfrac{1}{2} a(t)(\Delta t)^2$$

It's synchronized form, the time-step $\Delta t$ must be constant to maintain stability. If the initial velocity $v(t)=v(0)=0$, once we know the $x(t)$, we know the $a(t)$ from potential field, so we can calculate the next position $x(t+\Delta t)$, and $a(t+\Delta t)$ is also available. Accordingly, we can calculate the next velocity $v(t+\Delta t)$


### Advantages of Velocity Verlet

Positions $x(t)$, velocities $v(t)$, and accelerations $a(t)$ are available at the same instance.

There is no need to store more than one set of positions, velocities, and accelerations. We only need to know $v(0), x(0),a(0)$, and we can calculate the next step information: $v(\Delta t), x(\Delta t), a(\Delta t)$

## How does a system evolve

Integrate the equations of motion:

Calculate ${x(t+\Delta t),v(t+\Delta t), a(t+\Delta t)}$ from ${x(t),v(t), a(t)}$ 

Assign the initial velocities, such that:

The system has the desired temperature: $E_k=\frac{3}{2}NkT$, and the system has no net momentum: $\sum_{i=1}^{N}m_iv_i=0$

Let the system evolve, it will generally adjust the balance between kinetic and potential energy, so the temperature of the system changes.

Adjust the temperature by rescaling the velocities and let the system evolve again.

#### Note: violation of the following indecates a bug or a too large timestep $\Delta t$:

Global energy conservation: $E_k+E_p=constant$

Momentum conservation

Note: violation of the following indicates lack of equilibrium: 

Temperature fluctuates around a stable average

### How to measure properties

At regular intervals, record properties of interet, there are some usual static properties: $E_k, E_p, P$, and some unusual static properties: average force on a particle, and typical number and type of neighbors, and dynamic properties
















