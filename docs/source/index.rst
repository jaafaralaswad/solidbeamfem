solidbeamFEM Documentation
==========================

Nonlinear finite element solver in Python for beams modeled as continua using
brick elements in convected curvilinear coordinates. It alleviates
geometry-induced locking modes (membrane, transverse shear, and
curvature–thickness locking) using the Assumed Natural Strain (ANS) method.
Locking modes can be alleviated selectively to identify which modes are
responsible for the observed locking behavior. The solver also supports
comparisons with alternative techniques such as h- and p-refinement and
reduced integration rules.

.. image:: https://img.shields.io/badge/GitHub-Repository-black?logo=github
   :target: https://github.com/jaafaralaswad/solidbeamfem
   :alt: GitHub Repository

.. toctree::
   :maxdepth: 2
   :hidden:

   weak_form
   ans_method
   examples


Dependencies
============

This package depends only on:

- NumPy
- Matplotlib


Package Installation
====================

Users are expected to clone the GitHub repository and edit the main driver
script directly to define the problem.

Install by cloning the GitHub repository
---------------------------------------

.. code-block:: bash

   git clone https://github.com/jaafaralaswad/solidbeamfem.git
   cd solidbeamfem
   pip install .


Locate and edit the main driver script
-------------------------------------

All problem parameters are defined in:

.. code-block:: text

   src/solidbeamfem/main.py

This file serves as the main entry point of the solver. Open ``main.py`` and
modify the parameters under the *Problem setup* section. Typical parameters
include:

- Beam geometry (width, height, length)
- Material properties (Young’s modulus, Poisson’s ratio)
- Loading type and magnitude
- Mesh discretization parameters
- Numerical integration points
- Assumed Natural Strain (ANS) locking-alleviation switches
- Solver tolerances and load increments
- Visualization options

After editing the parameters, save the file.


Usage
=====

Using ``solidbeamfem`` only requires defining the problem parameters.
All numerical details—including finite element assembly, nonlinear solution,
locking alleviation, and post-processing—are handled internally.

A typical workflow consists of:

1. Defining the geometry and material properties of the beam
2. Specifying the loading
3. Choosing the discretization and numerical parameters
4. Enabling optional locking-alleviation techniques
5. Running the solver and visualizing the results

Below is a minimal example illustrating a complete problem definition.


Problem definition
------------------

.. code-block:: python

   # =============================================================================
   # Problem setup
   # =============================================================================

   # Geometry (rectangular cross-section)
   width = 1.0       # beam width
   height = 0.5      # beam height
   length = 20.0     # beam length

   # Material (Saint Venant–Kirchhoff)
   E = 1.2e7         # Young's modulus
   nu = 0.0          # Poisson's ratio (nu = 0 avoids volumetric locking)

   # Loading
   load_type = "moment"   # "moment" or "shear"
   k_max = 0.5            # target nondimensional load

   # Discretization
   numel = 10             # elements along beam length
   ne_L = 2               # nodes per element in axial direction (>= 2)

   # Numerical integration
   ngp_c = 2              # Gauss points per cross-sectional direction
   ngp_l = ne_L           # Gauss points along beam length (<= 10)

   # Locking alleviation (Assumed Natural Strain)
   ANS_membrane = True    # membrane locking
   ANS_shear = True       # transverse shear locking
   ANS_curvature = True  # curvature–thickness locking

   # Incremental-iterative solver
   n_load_steps = 10      # load increments
   max_iterations = 20   # Newton–Raphson iterations per step
   tolerance = 1e-15     # convergence tolerance on energy norm

   # Visualization
   visualize_reference = True   # visualize reference configuration
   visualize_final = True       # visualize final configuration
   plot_displacement = True     # plot normalized tip displacement


Running the solver
------------------

Once the problem parameters are defined, the solver is executed by running
the main driver script.

.. code-block:: bash

   python -m solidbeamfem.main

or equivalently:

.. code-block:: bash

   python src/solidbeamfem/main.py

This will:

- Generate the finite element mesh
- Assemble and solve the nonlinear equilibrium equations
- Apply Assumed Natural Strain (ANS) corrections if enabled
- Produce the requested visualizations and plots

More advanced modifications—such as changing constitutive models or boundary
conditions—require manual changes to the corresponding implementation files
within the package.