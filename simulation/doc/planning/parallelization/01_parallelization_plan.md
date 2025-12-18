# Parallelization Plan using MPI

This document outlines the steps to parallelize the existing fluid dynamics solver using MPI (Message Passing Interface). The goal is to distribute the computational domain across multiple processes (ranks) to speed up the simulation.

## Overview

The parallelization strategy involves **Domain Decomposition**. The global grid is split into smaller rectangular sub-grids, one for each MPI process. Each process solves the equations on its local sub-grid. To maintain continuity across boundaries, processes exchange data with their neighbors using **Ghost Layers**.

## Step 1: MPI Initialization

Before any parallel work can begin, MPI must be initialized.

**Action Items:**
1.  Modify `src/main.cpp`.
2.  Call `MPI_Init(&argc, &argv)` at the beginning.
3.  Call `MPI_Finalize()` at the end.
4.  Get the rank and number of processes using `MPI_Comm_rank` and `MPI_Comm_size`.

```cpp
// src/main.cpp
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  // ... existing code ...

  MPI_Finalize();
  return 0;
}
```

## Step 2: The Partitioning Class

The `Partitioning` class is responsible for determining the layout of the parallel grid. It decides which process gets which part of the domain and who its neighbors are.

### 2.1 Decomposition Strategy
We need to split the `nRanks` processes into a 2D grid of processes: `nProcX` by `nProcY`.
*   **Hint**: You can use the MPI function `MPI_Dims_create(nRanks, 2, dims)` to let MPI find a good division.
*   Alternatively, you can compute it manually (e.g., if `nRanks` is 4, `2x2` is usually better than `4x1`).

### 2.2 Implementing `src/partitioning/partitioning.cpp`
You need to implement the methods defined in `partitioning.h`.

**Key Concepts:**
*   **Global Domain**: The total `nCells` from `Settings`.
*   **Local Domain**: The subset of cells this rank owns.
*   **Offset**: The global coordinate of the first local cell.

**Implementation Steps:**
1.  **`initialize(nCellsGlobal)`**:
    *   Determine `nProcX` and `nProcY` (e.g., using `MPI_Dims_create`).
    *   Calculate the position of the current rank in this process grid: `rankX` and `rankY`.
        *   `rankX = ownRankNo % nProcX`
        *   `rankY = ownRankNo / nProcX`
    *   Calculate `nCellsLocal`:
        *   Divide `nCellsGlobal` by `nProcX` and `nProcY`.
        *   **Careful**: If it doesn't divide evenly, the last ranks need to take the remainder.
    *   Calculate `nodeOffset`:
        *   Sum of `nCellsLocal` of all ranks to the left/bottom.
    *   Identify Neighbors:
        *   `leftNeighbour`: `rankX > 0 ? ownRankNo - 1 : MPI_PROC_NULL`
        *   `rightNeighbour`: `rankX < nProcX - 1 ? ownRankNo + 1 : MPI_PROC_NULL`
        *   `topNeighbour`: `rankY < nProcY - 1 ? ownRankNo + nProcX : MPI_PROC_NULL`
        *   `bottomNeighbour`: `rankY > 0 ? ownRankNo - nProcX : MPI_PROC_NULL`
        *   **Note**: `MPI_PROC_NULL` is a special value indicating no neighbor (boundary of the global domain).

## Step 3: Adapting the Staggered Grid

The `StaggeredGrid` (and `Discretization`) needs to work with the **local** size, but be aware of the **ghost layers**.

### 3.1 Ghost Layers
Currently, `StaggeredGrid` has `_uBoundaries`, `_vBoundaries`, etc., set to `{1,1,1,1}`. This is already correct!
*   **Physical Boundary**: A boundary that touches the edge of the global domain.
*   **Parallel Boundary**: A boundary that touches a neighbor rank.

In our code, we treat both the same way in terms of storage: we allocate an extra layer of cells (ghost layer) around the local domain.
*   If it's a **Physical Boundary**, we set values based on Dirichlet/Neumann BCs (as done in `setBoundaryConditions...`).
*   If it's a **Parallel Boundary**, we receive values from the neighbor.

**Action Items:**
*   Ensure `Discretization` is initialized with `nCellsLocal` from the `Partitioning` object, not `nCellsGlobal`.
*   The `meshWidth` (dx, dy) remains the same (calculated from `physicalSize / nCellsGlobal`).

## Step 4: Communication (Ghost Layer Exchange)

This is the most critical part. We need a mechanism to update the ghost layers.

### 4.1 The Exchange Function
Create a new class or function (e.g., in `Partitioning` or a new `ParallelHelper`) to handle data exchange.

**Function Signature:**
```cpp
void exchange(FieldVariable& f, const Partitioning& partitioning);
```

**Algorithm:**
1.  **Pack Data**: Copy the innermost layer of your local data into a send buffer.
    *   Left edge to send to Left Neighbor.
    *   Right edge to send to Right Neighbor.
    *   Top edge to send to Top Neighbor.
    *   Bottom edge to send to Bottom Neighbor.
2.  **Send/Recv**: Use MPI to exchange these buffers.
    *   **Hint**: Use `MPI_Sendrecv` to prevent deadlocks. It sends to one neighbor and receives from another simultaneously.
    *   Alternatively, use `MPI_Isend` (non-blocking send) and `MPI_Irecv` (non-blocking receive), followed by `MPI_Waitall`.
3.  **Unpack Data**: Copy the received data from the receive buffer into your ghost layers.

**Order of Exchange:**
To avoid corner issues, it's often easiest to:
1.  Exchange Left/Right.
2.  Unpack.
3.  Exchange Top/Bottom (including the ghost cells you just received from Left/Right).
4.  Unpack.
This ensures corner ghost cells are correctly populated.

## Step 5: Integration into Simulation Loop

Modify `src/simulation.cpp` to include the `Partitioning` object and call exchanges.

**Modifications:**
1.  **Constructor**:
    *   Create `Partitioning` object.
    *   Call `partitioning.initialize(settings->nCells)`.
    *   Initialize `Discretization` with `partitioning.nCellsLocal()`.
2.  **`runTimestep()`**:
    *   **Velocity Exchange**: Before `computeIntermediateVelocities`, exchange `u` and `v`.
        ```cpp
        exchange(discretization_->u(), partitioning_);
        exchange(discretization_->v(), partitioning_);
        ```
    *   **F/G Exchange**: Before `computeRHS`, exchange `f` and `g`.
        ```cpp
        exchange(discretization_->f(), partitioning_);
        exchange(discretization_->g(), partitioning_);
        ```
    *   **Pressure Solve**: The pressure solver needs to exchange `p` *inside* its iteration loop (see Step 6).

## Step 6: Parallel Pressure Solver

The standard Gauss-Seidel or SOR solver iterates through the grid `(i, j)` and uses the *most recent* values. In parallel, "most recent" might be on another processor.

### 6.1 The Issue with Standard GS/SOR
If you run standard GS in parallel, the values at the boundaries will lag behind (using values from the previous iteration until an exchange happens). This changes the convergence behavior and might even diverge.

### 6.2 Solution: Red-Black Gauss-Seidel
A common parallel approach is **Red-Black (Checkerboard) Gauss-Seidel**.
1.  **Color the grid**: Like a checkerboard, cells are Red `(i+j is even)` or Black `(i+j is odd)`.
2.  **Update Red cells**: They only depend on Black neighbors.
3.  **Exchange**: Communicate ghost layers for `p`.
4.  **Update Black cells**: They only depend on Red neighbors (which are now updated).
5.  **Exchange**: Communicate ghost layers for `p`.

**Action Items:**
*   Modify `PressureSolver` (or create `RedBlackSORPressureSolver`).
*   Split the loop into two passes: one for even `i+j`, one for odd `i+j`.
*   Insert `exchange(p)` between the passes.

## Step 7: Output Writers

You have `OutputWriterParaviewParallel`.
*   Pass the `Partitioning` object to it.
*   It uses `MPI_File_write_at` or similar (or VTK XML formats with `PStructuredGrid`) to write the data.
*   Since this is provided, just ensure you pass the correct `Partitioning` object and that `nCellsGlobal` and `nodeOffset` are correct.
