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

  * physical vapor deposition [released] [1]_ [2]_;
  * chemical vapor deposition [released] [3]_;

* Melting / solidification processes at fixed temperature [released] [4]_ [5]_
* Nanosecond pulsed laser annealing processes [released] [4]_ [5]_
* Silicidation processes [under development]

References
==========

.. [1] A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, I. Deretzis, F. La Via, and E. Smecca, "Simulation of the Growth Kinetics in Group IV Compound Semiconductors", physica status solidi (a) vol. 216, no. 10, p. 1800597, 2019, doi: 10.1002/pssa.201800597
.. [2] G. Fisicaro, C. Bongiorno, I. Deretzis, F. Giannazzo, F. La Via, F. Roccaforte, M. Zielinski, M. Zimbone, A. La Magna, "Genesis and Evolution of Extended Defects: The Role of Evolving Interface Instabilities in Cubic SiC", Applied Physics Reviews vol. 7, no. 2, p. 021402, Apr. 2020, doi: 10.1063/1.5132300
.. [3] https://hq.imm.cnr.it/content/super-lattice-kinetic-monte-carlo-method-simulate-cvd-epitaxy-si-based-materials
.. [4] G. Calogero, D. Raciti, P. Acosta-Alba, F. Cristiano, I. Deretzis, G. Fisicaro, K. Huet, S. Kerdilès, A. Sciuto and A. La Magna, "Multiscale modeling of ultrafast melting phenomena", npj Computational Materials 8, 36 (2022), doi: 10.1038/s41524-022-00720-y
.. [5] G. Calogero, D. Raciti, D. Ricciarelli, P. Acosta-Alba, F. Cristiano, R. Daubriac, R. Demoulin, I. Deretzis, G. Fisicaro, J.-M. Hartmann, S. Kerdilés, A. La Magna, "Atomistic insights into ultrafast SiGe nanoprocessing" J. Phys. Chem. C 2023, 127, 19867-19877, doi: 10.1021/acs.jpcc.3c05999
