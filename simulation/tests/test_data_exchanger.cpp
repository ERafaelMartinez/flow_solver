#include "../src/communication/exchanger.h"
#include "../src/partitioning/partitioning.h"
#include "../src/storage/field_variable.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

/**
 * DataExchanger Test
 *
 * Tests MPI ghost cell exchange with 4 processes on a 10x10 global grid.
 * Each process gets a 5x5 local grid + 2 ghost cells = 7x7 total.
 *
 * Partitioning Layout (rank indexing convention):
 *   {2, 3}  <- Top row (j=1)
 *   {0, 1}  <- Bottom row (j=0)
 *
 * Neighbor relationships:
 *   Rank 0 (Bottom-Left):  Right=1, Top=2
 *   Rank 1 (Bottom-Right): Left=0,  Top=3
 *   Rank 2 (Top-Left):     Right=3, Bottom=0
 *   Rank 3 (Top-Right):    Left=2,  Bottom=1
 *
 * Each rank fills its internal cells with: rank*1000 + i*100 + j
 * After exchange, we verify that ghost cells contain the correct neighbor data.
 */

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  if (nRanks != 4) {
    if (rank == 0)
      std::cerr << "This test requires exactly 4 ranks." << std::endl;
    MPI_Finalize();
    return 1;
  }

  // Define global grid
  std::array<int, 2> nCellsGlobal = {10, 10}; // 5x5 per rank
  auto partitioning = std::make_shared<Partitioning>();
  partitioning->initialize(nCellsGlobal);

  // Define FieldVariable
  std::array<int, 2> nCellsLocal = partitioning->nCellsLocal();
  std::array<int, 2> gridSize = {nCellsLocal[0] + 2,
                                 nCellsLocal[1] + 2}; // +2 for ghost layers
  std::array<double, 2> offset = {0, 0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable u(gridSize, offset, cellSize);
  u.setToZero();

  // Fill internal domain with rank-specific pattern
  // Value = rank * 1000 + i * 100 + j
  // i [1, nCellsLocal[0]], j [1, nCellsLocal[1]]
  for (int j = 1; j <= nCellsLocal[1]; j++) {
    for (int i = 1; i <= nCellsLocal[0]; i++) {
      u.at(i, j) = rank * 1000.0 + i * 100.0 + j;
    }
  }

  // Exchange
  DataExchanger exchanger(partitioning);
  exchanger.exchange(u);

  // Verify ghost cells for all ranks
  // Layout: {2, 3}  (top row)
  //         {0, 1}  (bottom row)
  bool pass = true;

  // Rank 0: Bottom-Left (col=0, row=0)
  // Neighbors: Right=1, Top=2, Left=none, Bottom=none
  if (rank == 0) {
    // Check Right Ghost (column size[0]-1) -> Should come from Rank 1, column 1
    // Rank 1 value at (1,j) = 1000 + 100 + j
    for (int j = 1; j <= nCellsLocal[1]; j++) {
      double val = u.at(u.size()[0] - 1, j);
      double expected = 1000.0 + 100.0 + j;
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 0] Right Ghost Error at j=" << j << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }

    // Check Top Ghost (row size[1]-1) -> Should come from Rank 2, row 1
    // Rank 2 value at (i,1) = 2000 + i*100 + 1
    for (int i = 1; i <= nCellsLocal[0]; i++) {
      double val = u.at(i, u.size()[1] - 1);
      double expected = 2000.0 + i * 100.0 + 1.0;
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 0] Top Ghost Error at i=" << i << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }
  }

  // Rank 1: Bottom-Right (col=1, row=0)
  // Neighbors: Left=0, Top=3, Right=none, Bottom=none
  else if (rank == 1) {
    // Check Left Ghost (column 0) -> Should come from Rank 0, column size[0]-2
    // Rank 0 value at (size[0]-2, j) = 0 + (size[0]-2)*100 + j = 400 + j (since
    // size=7, so col=5)
    for (int j = 1; j <= nCellsLocal[1]; j++) {
      double val = u.at(0, j);
      double expected = 0.0 + (nCellsLocal[0]) * 100.0 +
                        j; // column index 5 in Rank 0's internal
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 1] Left Ghost Error at j=" << j << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }

    // Check Top Ghost (row size[1]-1) -> Should come from Rank 3, row 1
    // Rank 3 value at (i,1) = 3000 + i*100 + 1
    for (int i = 1; i <= nCellsLocal[0]; i++) {
      double val = u.at(i, u.size()[1] - 1);
      double expected = 3000.0 + i * 100.0 + 1.0;
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 1] Top Ghost Error at i=" << i << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }
  }

  // Rank 2: Top-Left (col=0, row=1)
  // Neighbors: Right=3, Bottom=0, Left=none, Top=none
  else if (rank == 2) {
    // Check Right Ghost (column size[0]-1) -> Should come from Rank 3, column 1
    // Rank 3 value at (1,j) = 3000 + 100 + j
    for (int j = 1; j <= nCellsLocal[1]; j++) {
      double val = u.at(u.size()[0] - 1, j);
      double expected = 3000.0 + 100.0 + j;
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 2] Right Ghost Error at j=" << j << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }

    // Check Bottom Ghost (row 0) -> Should come from Rank 0, row size[1]-2
    // Rank 0 value at (i, size[1]-2) = 0 + i*100 + (size[1]-2) = i*100 + 5
    for (int i = 1; i <= nCellsLocal[0]; i++) {
      double val = u.at(i, 0);
      double expected =
          0.0 + i * 100.0 + nCellsLocal[1]; // row index 5 in Rank 0's internal
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 2] Bottom Ghost Error at i=" << i << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }
  }

  // Rank 3: Top-Right (col=1, row=1)
  // Neighbors: Left=2, Bottom=1, Right=none, Top=none
  else if (rank == 3) {
    // Check Left Ghost (column 0) -> Should come from Rank 2, column size[0]-2
    // Rank 2 value at (size[0]-2, j) = 2000 + (size[0]-2)*100 + j = 2000 + 500
    // + j
    for (int j = 1; j <= nCellsLocal[1]; j++) {
      double val = u.at(0, j);
      double expected = 2000.0 + nCellsLocal[0] * 100.0 + j;
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 3] Left Ghost Error at j=" << j << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }

    // Check Bottom Ghost (row 0) -> Should come from Rank 1, row size[1]-2
    // Rank 1 value at (i, size[1]-2) = 1000 + i*100 + (size[1]-2) = 1000 +
    // i*100 + 5
    for (int i = 1; i <= nCellsLocal[0]; i++) {
      double val = u.at(i, 0);
      double expected = 1000.0 + i * 100.0 + nCellsLocal[1];
      if (std::abs(val - expected) > 1e-6) {
        std::cout << "[Rank 3] Bottom Ghost Error at i=" << i << ": got " << val
                  << " expected " << expected << std::endl;
        pass = false;
      }
    }
  }

  if (pass)
    std::cout << "Rank " << rank << " Verification Passed." << std::endl;
  else
    std::cout << "Rank " << rank << " Verification FAILED." << std::endl;

  MPI_Finalize();
  return pass ? 0 : 1;
}
