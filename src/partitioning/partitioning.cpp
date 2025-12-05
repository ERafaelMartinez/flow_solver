#include "./partitioning.h"
#include <array>
#include <cmath>
#include <mpi.h>

void Partitioning::initialize(std::array<int, 2> nCellsGlobal) {
  // 1. We want to figure out what our indices are
  // 2. How many neighbours we have

  // rank determines which section of the grid we are looking at
  int processes_count = this->nRanks();
  int own_rank = this->ownRankNo();

  nCellsGlobal_ = nCellsGlobal;

  // TODO: Find a better way that would work with other numbers (real squars,
  // close to squares, prime)
  int nTilesX = std::sqrt(processes_count);
  int nTilesY = std::ceil(processes_count / nTilesX);

  // determine the position (indices) of the partition
  int subdomain_col = own_rank / nTilesX;
  int subdomain_row = own_rank % nTilesX;

  // TODO: we still have to deal with edge cases where the size of the subdomain
  // is not square.
  nCellsLocal_[0] = computeSubdomainSizeInAxis(nCellsGlobal[0], nTilesX);
  nCellsLocal_[1] = computeSubdomainSizeInAxis(nCellsGlobal[1], nTilesY);

  // Determine the offsets
  // TODO: Check if the assignment (col vs row) is correct and how other classes
  // expect to use this
  nodeOffset_[0] =
      computeSubdomainAxisOffset(nCellsGlobal_[0], nTilesX, subdomain_col);
  nodeOffset_[1] =
      computeSubdomainAxisOffset(nCellsGlobal_[1], nTilesY, subdomain_row);

  // determine the neighbour partitions
  neighbourLeft_ =
      getProcessAt(subdomain_col - 1, subdomain_row, nTilesX, nTilesY);
  neighbourRight_ =
      getProcessAt(subdomain_col + 1, subdomain_row, nTilesX, nTilesY);
  neighbourBottom_ =
      getProcessAt(subdomain_col, subdomain_row - 1, nTilesX, nTilesY);
  neighbourTop_ =
      getProcessAt(subdomain_col, subdomain_row + 1, nTilesX, nTilesY);
}

int Partitioning::getProcessAt(int i, int j, int numOfColumns,
                               int numOfRows) const {
  if (i < 0 || i > numOfColumns - 1 || j < 0 || j > numOfRows - 1) {
    return -1;
  }

  return i * numOfColumns + j;
}

int Partitioning::computeSubdomainSizeInAxis(int globalSizeInAxis,
                                             int axisPartitionCount) const {
  int processesCount = this->nRanks();
  int ownRank = this->ownRankNo();

  int sizePerSubdomain = globalSizeInAxis / axisPartitionCount;

  int sizeInAxis = globalSizeInAxis;
  // If we can divide the domain into multiple subdomains, then figure out the
  // correct local size
  if (sizePerSubdomain > 1) {
    if (ownRank < processesCount - 1) {
      sizeInAxis = sizePerSubdomain;
    } else {
      // For the last subdomain, we take the remaining of the global domain that
      // was not assigned to any previous subdomain
      sizeInAxis = globalSizeInAxis - (processesCount - 1) * sizePerSubdomain;
    }
  }

  return sizeInAxis;
}

int Partitioning::computeSubdomainAxisOffset(int globalSizeInAxis,
                                             int axisPartitionCount,
                                             int partitionIndexInAxis) const {

  // this is enough for now, since we currently now that only the right most and
  // top most partitions could differ in size than most. HACK: This though can
  // change depending on the method we use to compute how we want to partition
  // the domain
  int squareSizeForLeftAndBottomMostPartitions =
      globalSizeInAxis / axisPartitionCount;
  return partitionIndexInAxis * squareSizeForLeftAndBottomMostPartitions;
}

//! get the local number of cells in the own subdomain
std::array<int, 2> Partitioning::nCellsLocal() const { return nCellsLocal_; }

//! get the global number of cells in the whole computational domain
//! used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nCellsGlobal() const { return nCellsGlobal_; }

//! get the own MPI rank no
//! used in OutputWriterParaviewParallel and OutputWriterTextParallel
int Partitioning::ownRankNo() const {
  int own_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &own_rank);
  return own_rank;
}

//! number of MPI ranks
int Partitioning::nRanks() const {
  // Get the number of processes
  int processes_count;
  MPI_Comm_size(MPI_COMM_WORLD, &processes_count);
  return processes_count;
}

//! if the own partition has part of the bottom boundary of the whole domain
bool Partitioning::ownPartitionContainsBottomBoundary() const {
  return neighbourBottom_ == -1;
}

//! if the own partition has part of the top boundary of the whole domain
//! used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsTopBoundary() const {
  return neighbourTop_ == -1;
}

//! if the own partition has part of the left boundary of the whole domain
bool Partitioning::ownPartitionContainsLeftBoundary() const {
  return neighbourLeft_ == -1;
}

//! if the own partition has part of the right boundary of the whole domain
//! used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsRightBoundary() const {
  return neighbourRight_ == -1;
}

//! get the rank no of the left neighbouring rank
int Partitioning::leftNeighbourRankNo() const { return neighbourLeft_; }

//! get the rank no of the right neighbouring rank
int Partitioning::rightNeighbourRankNo() const { return neighbourRight_; }

//! get the rank no of the top neighbouring rank
int Partitioning::topNeighbourRankNo() const { return neighbourTop_; }

//! get the rank no of the bottom neighbouring rank
int Partitioning::bottomNeighbourRankNo() const { return neighbourBottom_; }

//! get the offset values for counting local nodes in x and y direction.
//! (i_local,j_local) + nodeOffset = (i_global,j_global)
//! used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nodeOffset() const { return nodeOffset_; }
