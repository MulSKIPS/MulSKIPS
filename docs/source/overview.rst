The MulSKIPS code
=================

`MulSKIPS` is a Kinetic Monte Carlo super-Lattice code, designed to study with an atomic resolution
the growth kinetics of elements, alloys and compounds characterized by the sp3 bond symmetry.
The code is open source and it is distributed according to the GNU public license.
MulSKIPS is available on GitHub_, from where it can be downloaded as a tar file.

.. _GitHub: https://github.com/MulSKIPS/MulSKIPS

Deposition and evaporation of the substrate atoms are the active Monte Carlo events,
driving the stochastic evolution. In MulSKIPS, a dense super-lattice correctly
accommodates the original lattice of the ideal crystal along with a large class
of defective configurations. This feature makes the code unique in the range
of lattice Kinetic Monte Carlo codes currently available for sp3 materials.
Indeed, the code is able to simulate the evolution of both point and extended defects,
like stacking faults of different symmetries, antiphase boundaries and grain boundaries.
Moreover, MulSKIPS can simulate the morphological evolution during the growth process,
e.g. the epitaxial growth or etching of flat, structured, or patterned substrates,
as well as nanoparticles of various shapes.
In the case of surfaces, periodic boundary conditions are applied in the planes
orthogonal to the growth direction.

Features
========

Following we list all MulSKIPS functionalities:

* Epitaxial growth of nanoparticles, flat and patterned systems by:

  * physical vapor deposition [released] \ :cite:p:`LaMagna2019`,\ :cite:p:`Fisicaro2020`
  * chemical vapor deposition [released] 

* Melting / solidification processes at fixed temperature [released] \ :cite:p:`Calogero2022`,\ :cite:p:`Calogero2023`
* Nanosecond pulsed laser annealing processes [released] \ :cite:p:`Calogero2022`,\ :cite:p:`Calogero2023`
* Silicidation processes [under development]

References
==========

```{bibliography}
```