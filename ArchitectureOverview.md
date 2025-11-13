# Architecture Overview

The following document provides an overview of the architecture of the numerical simulation software. It outlines the main components, their responsibilities (i.e., functional requirements), and how they should interact with each other.

## Components
1. **Storage Layer**
    - **Responsibility**: Manage data storage and retrieval for various data structures used in the simulation.

    - **Key Classes**:

        - `Array2D`: A two-dimensional array class for storing scalar fields.
            - **Class Requirements**:
                - Provide read/write access to elements using row-major order.
                - An array should support both read and write operations (get, set).
                - An array should be initialized with a specified size and the size should be retrievable and remain constant.
                - Upon initialization, all elements should be set to zero.


2. **Discretization Layer**
    - **Responsibility**: Definition of the simulation domain for each system variable in a grid which represents the physical space. Each system variable may take values at different points in the cells of the grid (e.g., centered, on the edges, or on the corners), leading to the concept of staggered grids. This is useful for the discretization of partial differential equations governing the system, which is also part of this layer's responsibility.

    - **Key Classes**:
        - `FieldVariable`: Represents a physical variable defined on the simulation domain.
            - **Class Requirements**:
                - Each `FieldVariable` is represented by a 2D array of values corresponding to grid points, where the "scalar function" takes its values.

                - Provide methods to access and modify the variable's values at grid points.

                - Provide a method to interpolate values at arbitrary points within the domain using bilinear interpolation.

                - Support __initialization__ with a domain size specification (width, height),  grid spacing, and information about where in the grid-cell the variable takes its values: centered, up-down, or left-right. Upon initialization, all values should be set to zero.

                - **Note:** The `FieldVariable` class inherits from the `Array2D` class to utilize its storage capabilities. Its internal logic mainly consist on initializing the array with the correct size based on the domain size and grid spacing.


        - `StaggeredGrid`: Groups together multiple `FieldVariable` instances that are staggered with respect to each other. For our flow simulation, we have three main variables: pressure (p - centered), horizontal velocity (u - left-right), and vertical velocity (v - up-down), and two additional variables for intermediate calculations: F (left-right) and G (up-down).

            - **Class Requirements**:
                - For its initialization it receives optionally a set of FieldVariables (u, v, p, F, and G) to be staggered. Otherwise, it creates the field variables for the flow simulation with default parameters from the settings.

                - Upon initialization, if external variables are given, it validates that the variables can be staggered based on the size of the domain they represent, and the resolution (grid size).

                - Provide methods to access each of the staggered variables (u, v, p, F, G).

                - Based on the grid spacing, domain size, and system variables given, consider an additional cell padding around the domain to facilitate the application of __boundary conditions__.

                - **Note:** The class serves as an parent class to the SystemDiscretization class, where coordinated access to the variables is required to correctly compute the discretized equations.


        - `SystemDiscretization`: Manages the discretization of the entire system, coordinating multiple `FieldVariable` instances (is itself a `StaggeredGrid`, with additional methods for computations).

            - **Class Requirements**:
                - Provide an interface for classes which implement specific discretization schemes.

                - Should consider virtual methods for computing the operators required for the complete discretization of the Navier-Stokes differential equations:
                    - spatial derivatives:
                        - pressure gradients:
                            - p_x (left-right)
                            - p_y (up-down)

                        - diffusion terms:
                            - u_xx (left-right)
                            - u_yy (left-right)
                            - v_xx (up-down)
                            - v_yy (up-down)

                        - convection terms:
                            - u2_x: (left-right) - squared velocity component
                            - v2_y: (up-down) - squared velocity component
                            - (u*v)_y (left-right)
                            - (u*v)_x (up-down)

                        - continuity equation:
                            - u_x (centered)
                            - v_y (centered) 


                - Ensure that the discretization respects the staggered arrangement of the variables as specified in the point above.

                - Support initialization with a `StaggeredGrid` instance, allowing access to all field variables for discretization operations through inheritance.


        - `FiniteDifferencesDiscretization`: Implements the `SystemDiscretization` interface using finite difference methods for spatial discretization.

            - **Class Requirements**:
                - Provide concrete implementations for all virtual methods defined in the `SystemDiscretization` class, using finite difference schemes (e.g., central differences, upwind schemes) to compute the required spatial derivatives.
                
    


3. **Simulation Layer**
    - **Responsibility**: Manage the overall simulation process, including time-stepping, applying boundary conditions, and updating field variables based on the discretized equations. The main goal is the solution of the pressure and velocity fields over time by solving the discretized system of equations represented using the discretization layer.

    - **Key Classes**:

        - `PressureSolver`: Solves the pressure Poisson equation utilizing any appropriate numerical method (e.g., Successive Over-Relaxation, Gauss-Seidel, etc.).

            - **Class Requirements**:
                - Initialize the solver with a `SystemDiscretization` instance to access the necessary field variables and discretized operators.

                - Implement the system's differential equations in discretized form to set up the pressure Poisson equation.

                - Declare method(s) to apply/set boundary conditions specific to the pressure and velocity fields.

                - Declare a method which solves a single iteration of the pressure Poisson equation.

                - Declare a method to perform multiple iterations until convergence is reached or a maximum number of iterations is met.

                - Declare method to check if a convergence criteria is met.

                - Provide methods to retrieve the current state of the pressure and velocity fields.

                - **Note:** The actual implementation of the pressure solver algorithm will be done in subclasses that inherit from this class.


        - `GaussSeidelPressureSolver`: Inherits from `PressureSolver` and implements the Gauss-Seidel method for solving the pressure Poisson equation for a given time step.

            - **Class Requirements**:
                - Provide a concrete implementation of the pressure Poisson equation solver using the Gauss-Seidel iterative method.

                - Implement convergence criteria to determine when the solution has sufficiently converged.

                - Allow configuration of solver parameters such as relaxation factor and maximum iterations from the simulation settings

        - `SORPressureSolver`: Inherits from `PressureSolver` and implements the Successive Over-Relaxation (SOR) method for solving the pressure Poisson equation for a given time step.

            - **Class Requirements**:
                - Provide a concrete implementation of the pressure Poisson equation solver using the SOR iterative method.

                - Implement convergence criteria to determine when the solution has sufficiently converged.

                - Allow configuration of solver parameters such as relaxation factor and maximum iterations from the simulation settings
        
4. **Settings Layer**
    - **Responsibility**: Manage configuration settings for the simulation, including parameters such as grid size, maximum time step size, physical properties, and solver options.

    - **Key Classes**:

        - `Settings`: Encapsulates all configuration parameters required for the simulation.

            - **Class Requirements**:
                - Provide methods to load settings from a configuration .txt file.

                - Include all parameters mentioned in the given parameters.txt template, such as maximum time step size, reynolds number, and solver options (e.g., choice of pressure solver, convergence criteria).

                - Ensure that settings are validated upon loading to prevent invalid configurations.

                - Provide methods to retrieve individual settings as needed by other components of the simulation.

5. **Output Layer**
    - **Responsibility**: Handle the output of simulation results, including writing data to files in vtk format for visualization and analysis.

    - **Key Classes**:

        - `VTKWriter`: Manages the writing of simulation data to VTK files.

            - **Class Requirements**:
                - Provide methods to write the current state of field variables (pressure, velocity components) to VTK files.

                - Support configuration of output frequency and file naming conventions based on simulation settings.

                - Ensure that the output format is compatible with common visualization tools that support VTK files.


## Interaction Between Components

The following points summarize how the different components interact with each other. These reflects an overview of the logic in the main function of the simulation software:

- With the use of the **Settings Layer**, the simulation parameters are initialized and provided to other components.

- The **Discretization Layer** uses the **Storage Layer** to manage the data structures (2D arrays) for field variables and serves as the foundation for defining the simulation domain and discretizing the governing equations.

- In the **Simulation Layer**, the PressureSolver utilizes the operators defined in the **Discretization Layer** to explicitly compute the discretized version of the Navier-Stokes equations. This allows the PressureSolver to perform numerical computations and advance the simulation in time with the selected pressure solver.

- At the end of each time step, the **Output Layer** is invoked to write the current state of the simulation to VTK files for visualization and analysis.

## Conclusion

This architecture overview provides a high-level understanding of the main components of the numerical simulation software, their responsibilities, and how they interact with each other. Each component is designed to be modular, allowing for easy maintenance and potential future extensions or modifications.
