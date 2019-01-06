# xStellTools, v0.1b

A collection of scripts, codes and documents related to stellarators.

See also: http://vmecwiki.pppl.wikispaces.net/

* Biot Savart: A set of matlab codes that can be adapted to calculate the magnetic field components (Bx, By, Bz, Br, Bphi) at any point in space due to currents in a set of coils. Coils are described by sets of connected straight-line elements, although analytic functions could be implemented.  Examples for the Helically Symmetric Experiments, the Lithium Tokamak Experiment are shown. The code can easily be ported to support other experiments.
UPDATE: Old examples are current being updated to incorporate recent changes to the Biot-Savart function support. Some examples may not run correctly.

    The method of [J. Hanson, S. Hirshman PoP 9, 4410 (2002)] is used to rapidly and efficiently calculate the magnetic field components due to currents in the coilset. Check are included to see if the distance between the observxeation point and a coil is less than some specified tolerance.

