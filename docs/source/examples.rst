Examples
========

Nuemrical Example 1
-------------------

The problem we solve here is a rectangular cantiliever beam subjected to bending moment applied at its tip.


.. figure:: images/cantilever_bending.jpg
   :width: 65%
   :align: center
   :alt: Cantilever beam under bending

We simulate a cantilever beam using the following parameters. The load level is
set to \( k_{\max} = 0.5 \), corresponding to an applied bending moment
\( M = \pi E I / L \). According to classical elasticity theory, this moment
induces a constant curvature that bends the beam into a half-circle in the
absence of shear deformation.

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
   numel = 20             # elements along beam length
   ne_L = 2               # nodes per element in axial direction (>= 2)

   # Numerical integration
   ngp_c = 2              # Gauss points per cross-sectional direction
   ngp_l = ne_L           # Gauss points along beam length (<= 10)

   # Locking alleviation (Assumed Natural Strain)
   ANS_membrane = False    # membrane locking
   ANS_shear = False       # transverse shear locking
   ANS_curvature = False  # curvature–thickness locking

   # Incremental-iterative solver
   n_load_steps = 10      # load increments
   max_iterations = 20   # Newton–Raphson iterations per step
   tolerance = 1e-15     # convergence tolerance on energy norm

   # Visualization
   visualize_reference = True   # visualize reference configuration
   visualize_final = True       # visualize final configuration
   plot_displacement = True     # plot normalized tip displacement


The solution is obviously locked.

.. figure:: images/locked_config_ex1.png
   :width: 65%
   :align: center
   :alt: Locked deformed configuration for bent beam


.. figure:: images/locked_plot_ex1.png
   :width: 65%
   :align: center
   :alt: Locked plot configuration for bent beam


Now, we activate the locking alleviation techniques

.. code-block:: python

   # Locking alleviation (Assumed Natural Strain)
   ANS_membrane = True    # membrane locking
   ANS_shear = True       # transverse shear locking
   ANS_curvature = True  # curvature–thickness locking


The solution now is in excellent agreement with analytical solution.

.. figure:: images/correct_config_ex1.png
   :width: 65%
   :align: center
   :alt: Correct deformed configuration for bent beam


.. figure:: images/correct_plot_ex1.png
   :width: 65%
   :align: center
   :alt: Correct plot configuration for bent beam


The influence of the different locking modes can be assessed by selectively
activating the corresponding alleviation mechanisms. For this problem, the
response is dominated by transverse shear locking.

.. figure:: images/comparison_example_1.png
   :width: 65%
   :align: center
   :alt: Locking modes comparison for example 1





Nuemrical Example 2
-------------------


.. code-block:: python

   # =============================================================================
   # Problem setup
   # =============================================================================

   # Geometry (rectangular cross-section)
   width = 1.0       # beam width
   height = 1.0      # beam height
   length = 300.0     # beam length

   # Material (Saint Venant–Kirchhoff)
   E = 12         # Young's modulus
   nu = 0.0          # Poisson's ratio (nu = 0 avoids volumetric locking)

   # Loading
   load_type = "shear"   # "moment" or "shear"
   k_max = 4.0            # target nondimensional load

   # Discretization
   numel = 3             # elements along beam length
   ne_L = 3               # nodes per element in axial direction (>= 2)

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
