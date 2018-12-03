Magnetic spring with oblique magnets
====================================

This directory contains code for producing the graphs presented in the
publication `Theoretical analysis of a non-contact spring with inclined permanent magnets for load-independent resonance frequency`.

The figures in the manuscript are generated with the following files:

* Figures 2–3: `oblique angles.m`
* Figure 4: `oblique gaps.m`
* Figure 5: `oblique h.m`
* Figure 6: `oblique gaps angles.m`
* Figure 7: `oblique vol.m`
* Figures  8–11: `oblique_torques.m`
* Figures 12–14: `oblique_stability.m`

Schematics drawn in `oblique_torques.m` require the [fletcher](https://github.com/zprime/fletcher) Matlab package to be installed.

Usage and explanatory comments
------------------------------

The examples above use three main functions to do their work:

* `oblique_forces.m`
* `oblique_forces3.m`
* `oblique_dynamics.m`

These will be described in more detail below.
For all three of these, `help` comments are provided; e.g., type `doc oblique_forces` for more information.
If things don't make sense, please ask me!
I write the comments largely for my own use and some things may not be clear.

### `oblique_forces.m` ###

This function calculates static forces for the oblique magnet spring assuming zero rotation; translation is possible in all three translational degrees of freedom.

The function uses the `magnetforces.m` function that is included in the `magcode/matlab/` directory.

### `oblique_forces3.m` ###

This possibly ill-named function calculates static forces and torques for the oblique magnet spring for planar motion and rotation; that is, allowed motions are horizontal and vertical translation and rotation around the out-of-plane direction.

This function also has a large amount of code written to produce nice graphics or schematics to help illustrate the forces and torque in a given position.
(As shown in the second example in `oblique_torques.m`, animations are also possible.)

Also, the forces and torques are hand-tuned for speed (essentially having been taken directly from `magnetforces.m` and optimised a little) to aid the execution speed of `oblique_dynamics.m`.

### `oblique_dynamics.m` ###

This function solves the dynamic equations of motion for the planar case, using Runge-Kutta (`ode45`) with the output of `oblique_forces3.m`.
Various constraints can be applied based on the initial perturbation of the system.
(Long story short: it's not an entirely stable system.)


Licence for modification and distribution
-----------------------------------------

This work is freely modifiable and distributable under the terms and conditions of the
[Apache License v2](http://www.apache.org/licenses/LICENSE-2.0).
In effect, you are free to do with this code as you wish for the development of free or proprietary software.
Distributions of modified works must retain the original copyright notices and contain a list of modifications made, but this is only an expository notice; please see the licence text for complete details.

----------------------------------
Copyright 2010-2011 Will Robertson  
wspr81 at gmail dot com