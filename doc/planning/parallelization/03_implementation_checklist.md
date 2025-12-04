# Implementation Checklist

Use this checklist to track your progress.

## Phase 1: Setup & Partitioning
- [ ] **MPI Init**: Add `MPI_Init` and `MPI_Finalize` to `src/main.cpp`.
- [ ] **Partitioning Class**: Implement `src/partitioning/partitioning.cpp`.
    - [ ] `initialize`: Compute `nProcX`, `nProcY`, `ownRankNo`, `nCellsLocal`, `nodeOffset`.
    - [ ] `initialize`: Determine neighbor ranks (left, right, top, bottom).
    - [ ] Implement getters (`nCellsLocal`, `nodeOffset`, etc.).
- [ ] **Verify Partitioning**: Write a small test in `main.cpp` to print each rank's local size and neighbors.

## Phase 2: Grid & Data Structures
- [ ] **Simulation Constructor**: Initialize `Partitioning` object.
- [ ] **Discretization**: Pass `nCellsLocal` (from Partitioning) instead of `nCellsGlobal` to `Discretization`/`StaggeredGrid`.
- [ ] **Verify Memory**: Ensure arrays are allocated with size `nCellsLocal + ghostLayers`. (This should happen automatically if `StaggeredGrid` uses the size passed to it).

## Phase 3: Communication
- [ ] **Exchange Function**: Implement a helper function to exchange ghost layers.
    - [ ] Pack data buffers (left, right, top, bottom).
    - [ ] `MPI_Sendrecv` (or `Isend`/`Irecv`) for left/right.
    - [ ] Unpack left/right.
    - [ ] `MPI_Sendrecv` for top/bottom.
    - [ ] Unpack top/bottom.
- [ ] **Verify Communication**: Create a test case where each rank fills its grid with `rankNo`. Exchange and check if ghost layers contain neighbor's `rankNo`.

## Phase 4: Integration
- [ ] **Velocity Exchange**: Add `exchange(u)` and `exchange(v)` in `runTimestep` (before `computeIntermediateVelocities`).
- [ ] **F/G Exchange**: Add `exchange(f)` and `exchange(g)` in `runTimestep` (before `computeRHS`).
- [ ] **Pressure Solver**:
    - [ ] Modify `PressureSolver` to take `Partitioning` object.
    - [ ] Implement Red-Black SOR (or just add `exchange(p)` inside the loop for standard SOR/GS, though convergence might be poor).
    - [ ] Add `exchange(p)` inside the solver loop.

## Phase 5: Output & Validation
- [ ] **Output**: Connect `OutputWriterParaviewParallel`.
- [ ] **Run**: Run with `mpirun -n 4 ./numsim`.
- [ ] **Visualize**: Open in Paraview (it handles parallel `.pvtu` or `.vtu` files).
- [ ] **Validate**: Compare results with the serial version (e.g., Cavity flow benchmark).
