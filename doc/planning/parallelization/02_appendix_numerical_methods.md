# Appendix: Faster Numerical Methods for Parallelization

The standard Gauss-Seidel or SOR solvers are simple to implement but suffer from slow convergence, especially as the grid size increases. Their convergence rate degrades as $O(N^2)$. In a parallel setting, we need methods that are not only fast but also minimize communication overhead.

## 1. Red-Black Gauss-Seidel (SOR)

As mentioned in the main plan, Red-Black ordering allows for parallel execution of the Gauss-Seidel method.

*   **Pros**: Easy to implement, good parallel efficiency (communication only twice per iteration).
*   **Cons**: Convergence is still slow for large grids.

## 2. Conjugate Gradient (CG) Method

The Conjugate Gradient method is an iterative solver for systems of linear equations $Ax = b$ where $A$ is symmetric and positive-definite (SPD). The discretized Poisson equation for pressure is SPD.

### How it works
CG searches for the solution in "Krylov subspaces". It requires matrix-vector multiplications and dot products.

### Parallel Implementation
*   **Matrix-Vector Product**: Requires ghost layer exchange (neighbor communication).
*   **Dot Product**: Requires a global sum of local dot products.
    *   **MPI Command**: `MPI_Allreduce(local_sum, global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD)`.

### Performance
CG generally converges much faster than SOR for this type of problem, typically in $O(N^{1.5})$ or $O(N)$ with a good preconditioner.

## 3. Multigrid Methods (Geometric Multigrid)

Multigrid (MG) is one of the most efficient solvers for elliptic partial differential equations like the Poisson equation. It can achieve optimal $O(N)$ complexity (linear time).

### The Concept
Standard iterative solvers (like SOR) are good at smoothing out high-frequency errors (local noise) but very slow at reducing low-frequency errors (global offsets).
Multigrid solves this by using a hierarchy of grids:
1.  **Smoothing**: Run a few iterations of SOR on the fine grid to remove high-frequency error.
2.  **Restriction**: Transfer the residual (error) to a coarser grid (e.g., half the resolution).
3.  **Solve Coarse**: Solve for the error on the coarse grid (recursively).
4.  **Prolongation**: Interpolate the correction back to the fine grid.
5.  **Correction**: Update the fine grid solution.

### Parallel Implementation
*   Each level of the grid is partitioned among the processors.
*   Communication is required during smoothing (ghost layers) and during restriction/prolongation.
*   **Challenge**: As grids get very coarse, some processors might have very few or zero cells.

## Recommendation

For this project, start with **Red-Black SOR**. It is a direct adaptation of your existing code and teaches the fundamentals of parallel patterns.
If performance is still a bottleneck, **Conjugate Gradient** is the next logical step as it is easier to implement than Multigrid while offering significant speedups.
