![Python Version](https://img.shields.io/badge/python-3.12-blue)
![OS](https://img.shields.io/badge/os-ubuntu%20%7C%20macos%20%7C%20windows-blue)
![License](https://img.shields.io/badge/license-MIT-green)
[![Documentation Status](https://readthedocs.org/projects/solidbeamfem/badge/?version=latest)](https://solidbeamfem.readthedocs.io/en/latest/)

# solidbeamFEM
Nonlinear finite element solver in Python for beams modeled as continua using brick elements in convected curvilinear coordinates. It alleviates geometry-induced locking modes (membrane, transverse shear, and curvature–thickness locking) using the Assumed Natural Strain (ANS) method. Locking modes can be alleviated selectively to identify which modes are responsible for the observed locking behavior. The solver also supports comparisons with alternative techniques such as h- and p-refinement and reduced integration rules. Additional information and example problems can be found in the ReadTheDocs documentation:

<div align="center">

**[Read the Docs Documentation](https://solidbeamfem.readthedocs.io)**

</div>

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Community Contributions](#community-contributions)


## Dependencies

This package depends only on:

- NumPy
- Matplotlib

The code has been tested with Python 3.12.


## Installation

Please refer to the documentation for detailed installation and usage instructions.

**[Read the Docs Documentation](https://solidbeamfem.readthedocs.io)**


## Community Contributions

Community contributions are welcome.

If you encounter any issues, please report them using the
[GitHub Issues tracker](https://github.com/jaafaralaswad/solidbeamfem/issues).

For efficiency, please include:

- A minimal working example that reproduces the issue
- The full error message or traceback