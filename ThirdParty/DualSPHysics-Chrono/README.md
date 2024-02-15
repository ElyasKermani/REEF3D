# DualSPHysics-Chrono
This project contains the source code, dependencies and examples of **DualSPHysics v5.0.231** coupled to **Project Chrono v4.0.0**. This coupling is done via a communication interface (the so-called **DSPHChronoLib**).

## DualSPHysics
DualSPHysics ([https://dual.sphysics.org/](https://dual.sphysics.org/)) is based on the Smoothed Particle Hydrodynamics (SPH) method.
The code is developed to study free-surface flow phenomena where Eulerian methods can be difficult to apply, such as waves or impact of dam-breaks on off-shore structures. DualSPHysics is a set of C++, CUDA and Java codes designed to deal with real-life engineering problems. 

The source code is freely available in its public repository ([https://github.com/DualSPHysics/DualSPHysics](https://github.com/DualSPHysics/DualSPHysics)).

The list of developers of DualSPHysics can be seen on its website ([https://dual.sphysics.org/developers/](https://dual.sphysics.org/developers/)).

The package DualSPHysics-Chrono contains the source code of DualSPHysics v5.0.231, which is released under the GNU Lesser General Public License (LGPL) (see [https://github.com/DualSPHysics/DualSPHysics/blob/master/LICENSE](https://github.com/DualSPHysics/DualSPHysics/blob/master/LICENSE)).

## Project Chrono
Chrono ([https://www.projectchrono.org/](https://www.projectchrono.org/)) is a physics-based modelling and simulation infrastructure based on a platform-independent open-source design implemented in C++. A Project Chrono library can be embedded in a software project to simulate, for instance, wheeled and tracked vehicles operating on deformable terrains, robots, mechatronic systems, compliant mechanisms, and fluid solid interaction phenomena. Systems can be made of rigid and flexible/compliant parts with constraints, motors and contacts; parts can have three-dimensional shapes for collision detection. 

The source code of Chrono is freely available in its public repository ([https://github.com/projectchrono/chrono](https://github.com/projectchrono/chrono)).

The list of developers of Project Chrono can be seen on its website ([https://projectchrono.org/about/](https://projectchrono.org/about/)).

The DualSPHysics-Chrono package only contains the already compiled Chrono v4.0.0 executables. The source code of Chrono is released under a BSD-3-Clause license (see [https://github.com/projectchrono/chrono/blob/main/LICENSE](https://github.com/projectchrono/chrono/blob/main/LICENSE)).

## DSPHChronoLib
DSPHChronoLib is a source code developed in C++ to be used as a communication interface in charge of the information transfer between DualSPHysics and Chrono.
DSPHChronoLib is released under GNU Lesser General Public License (LGPL) and its Copyright can be seen in _DualSPHysics-Chrono/src_extra/DSPH-Chrono-Lib/LICENSE_.

**List of authors**

 Iván Martínez Estévez (ivan.martinez.estevez@uvigo.es). Universidade de Vigo, Spain.

 Dr José M. Domínguez (jmdominguez@uvigo.es). Universidade de Vigo, Spain.

 Dr Ricardo Canelas (ricardo.canelas@bentley.com). Bentley Systems, Lisbon, Portugal.

 Dr Bonaventura Tagliafierro (btagliafierro@gmail.com). University of Salerno, Italy.

 Dr Orlando García Feal (orlando@uvigo.es). Universidade de Vigo, Spain.

 Professor Alejandro J.C. Crespo (alexbexe@uvigo.es). Universidade de Vigo, Spain.

 Professor Moncho Gómez Gesteira (mggesteira@uvigo.es). Universidade de Vigo, Spain. 

## Documentation
**Note**: The complete documentation of this package is available in the following files:
- _Guide_DualSPHysics-Chrono.pdf_: Includes all the information to execute examples and compile the code.
- _Files_DualSPHysics-Chrono.pdf_: Includes the structure of directories and files of this package.
