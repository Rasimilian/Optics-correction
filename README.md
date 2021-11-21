# Optics-correction

A program for the magnetic structure correction (linear optics) of the VEPP-4M collider.

The program is intended for beam-based correction (errors finding) and optimization (fitting to a design model). Actually, orbit correction is possible as well.

It uses a response matrix approach with optimization algorithms: Gauss-Newton, Marquardt-Levenberg and supplemented with regularization (constraints).

Possible parameters to be varied are quadrupole gradients/shifts/rolls, dipole gradients, BPM shifts/rolls.
