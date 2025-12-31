Usage
=====

Using ``solidbeamfem`` only requires defining the problem.
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

Once the problem is defined, the solver can be executed directly.

.. code-block:: python

   import solidbeamfem

   solidbeamfem.main()

This will:

- Generate the finite element mesh
- Assemble and solve the nonlinear equilibrium equations
- Apply Assumed Natural Strain (ANS) corrections if enabled
- Produce the requested visualizations and plots

More advanced modifications—such as changing constitutive models or boundary conditions—require manual changes to the corresponding
implementation files within the package.