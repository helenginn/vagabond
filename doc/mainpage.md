\mainpage Vagabond documentation v0.0
 
##  Introduction

Vagabond performs model refinement against X-ray crystallographic data by
treating protein structures as a bond-based system, rather than an atomic one,
and combines position and flexibility information into a single model refined
against the electron density in real space. In doing so, Vagabond allows us
to drop the number of parameters used in refinement by fully exploiting the
prior biological information and produce maps with less overfitting bias. It
also allows unusual atom probability distributions to be generated due to use
of ensembles.

It accepts input as MTZ (reflection) and PDB (atomic model) files, and produces
outputs (Vagabond .vbond model files, MTZ files to generate real space maps,
and averaged PDB files) along with various diagnostic results.

##  Constraints

Vagabond is in its infancy and will frequently encounter situations it does
not handle well. Please see the list of \link limits limitations. \endlink

##  Installation

Right now, you will need to get the latest version from Github and compile 
from source (https://www.github.com/helenginn/vagabond). Vagabond is not
available and will crash on Windows.

Vagabond uses the Meson build system (https://www.mesonbuild.com) and has
the following dependencies (don't forget to install the _devel_ versions):
* libpng (usually included)
* fftw3f
* Qt5 (GUI only)

Please refer to the full \link install installation instructions \endlink
for details on how to install Vagabond using meson.

##  Tutorial

A tutorial on Vagabond would be sensible, but maybe once I've figured out
how the GUI will be settled.

##  Manual

There will be a manual on Vagabond.

##  Developer documentation

Documentation on the code is available. Start by having a look at the 
Classes.

##  License

Vagabond is released under GPL v3.



