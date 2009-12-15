Calculating the force and torque between magnets
================================================

This repository contains code originally written for my PhD for calculating the forces (and stiffnesses) between permanent magnets.

It has been originally written in Matlab, and I plan on adding translations into other programming languages as time permits. Contributions happily accepted.

Research
--------

The underlying maths of this work was not done by me but primarily by a group of French researchers over the last twenty-five years. To the best of my knowledge, this is the first attempt at presenting the work in a unified way, wrapped into a single interface.

Philosophy
----------

As time goes by, it becomes increasingly hard to build on the work of our predecessors unless we build tools to abstract the ideas we invent in our research. This code marks my first attempt to freely share, in a useful way, the work of my PhD.

Installation and information
----------------------------

The `matlab/` subdirectory of this repository contains both the source file
`magnetforces.web` and the ready-to-use Matlab functions `magnetforces.m` and
`multipoleforces.m`. After cloning this Git repository, you can simply add the `matlab/` folder to your Matlab path in order to start using the code.

The documentation `magnetocode-matlab.pdf` contains both a user's guide and the documented source code in a literate programming style.

Examples and test suite
-----------------------

The `examples/` subdirectory contains a number of files to illustrate the use of the `magnetforces` code (and other related things). These files are:

- `allag_torques.m` is an implementation of the theory of [Allag and Yonnet (2009)][1] for calculating the torques between cuboid magnets. I believe their theory is incorrect and this example is supposed to illustrate why.
- `magnetforces_example.m` contains a reproduction of the results of [Akoun and Yonnet (1984)][2] and of [Janssen et al. (2009)][3].
- `multipole_compare.m` is an unfinished comparison between different configurations of multipole arrays.
- `multipole_example.m` is a reproduction of the multipole results shown by [Allag, Yonnet, and Latreche (2009)][4].

[1]: http://dx.doi.org/10.1109/TMAG.2009.2025047
[2]: http://dx.doi.org/10.1109/TMAG.1984.1063554
[3]: http://dx.doi.org/10.1166/sl.2009.1049
[4]: http://dx.doi.org/10.1109/ELECTROMOTION.2009.5259084

The source code contains a number of automated tests to ensure that future changes don't break existing functionality or start producing incorrect results. These are not included in the repository for clarity; extract them with `mtangle`. 

Licence
-------

This work is freely modifiable and distributable under the terms and conditions of the Apache License v2. In effect, you are free to do with this code as you wish.

---------------------
Will Robertson  
wspr81 at gmail dot com