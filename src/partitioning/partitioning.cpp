#include "./partitioning.h"
#include <array>
#include <cmath>
#include <mpi.h>

void Partitioning::initialize(std::array<int, 2> nCellsGlobal) {
  // 1. We want to figure out what our indices are
  // 2. How many neighbours we have

  // rank determines which section of the grid we are looking at
  int processesCount = this->nRanks();
  int ownRank = this->ownRankNo();

  nCellsGlobal_ = nCellsGlobal;

  // determine the number of partitions/subdomains in x and y direction
  int nSubdomainsX = std::sqrt(processesCount);
  int nSubdomainsY = std::ceil(processesCount / nSubdomainsX);
  // TODO: Find a better way that would work with other numbers
  // (real squares, close to squares, prime)

  // determine the position (indices) of the partition
  int subdomainColIndex = ownRank / nSubdomainsY;
  int subdomainRowIndex = ownRank % nSubdomainsY;

  nCellsLocal_[0] = computeSubdomainSizeInAxis(nCellsGlobal[0], nSubdomainsX,
                                               subdomainColIndex);
  nCellsLocal_[1] = computeSubdomainSizeInAxis(nCellsGlobal[1], nSubdomainsY,
                                               subdomainRowIndex);

  // Determine the offsets
  // TODO: Check if the assignment (col vs row) is correct and how other classes
  // expect to use this
  nodeOffset_[0] = computeSubdomainAxisOffset(nCellsGlobal_[0], nSubdomainsX,
                                              subdomainColIndex);
  nodeOffset_[1] = computeSubdomainAxisOffset(nCellsGlobal_[1], nSubdomainsY,
                                              subdomainRowIndex);

  // determine the neighbour partitions
  neighbourLeft_ = getProcessAt(subdomainColIndex - 1, subdomainRowIndex,
                                nSubdomainsX, nSubdomainsY);
  neighbourRight_ = getProcessAt(subdomainColIndex + 1, subdomainRowIndex,
                                 nSubdomainsX, nSubdomainsY);
  neighbourBottom_ = getProcessAt(subdomainColIndex, subdomainRowIndex - 1,
                                  nSubdomainsX, nSubdomainsY);
  neighbourTop_ = getProcessAt(subdomainColIndex, subdomainRowIndex + 1,
                               nSubdomainsX, nSubdomainsY);
}

int Partitioning::getProcessAt(int i, int j, int numOfColumns,
                               int numOfRows) const {
  if (i < 0 || i > numOfColumns - 1 || j < 0 || j > numOfRows - 1) {
    return -1;
  }

  return i * numOfColumns + j;
}

int Partitioning::computeSubdomainSizeInAxis(int globalCellCountInAxis,
                                             int axisPartitionCount,
                                             int partitionIndexInAxis) const {
  int baseSize = globalCellCountInAxis / axisPartitionCount;
  int remainder = globalCellCountInAxis % axisPartitionCount;

  int sizeInAxis = baseSize;
  if (partitionIndexInAxis == axisPartitionCount - 1) {
    // The last partition gets the base size plus any remainder
    sizeInAxis += remainder;
  }
  return sizeInAxis;
}

int Partitioning::computeSubdomainAxisOffset(int globalCellCountInAxis,
                                             int axisPartitionCount,
                                             int partitionIndexInAxis) const {
  /*
   * calculate the offset of the partition in the axis based on the partition
   * index in the axis and the subdomain size distribution (made by the
   * computeSubdomainSizeInAxis method). This allocates the base size for each
   * partition and then adds the remainder to the last partition, thus for each
   * subdomain we get an integer multiple of the base size as offset.

   * NOTE: If computeSubdomainSizeInAxis changes, this method has to be updated
   * as well.
   */
  int baseSize = globalCellCountInAxis / axisPartitionCount;
  return partitionIndexInAxis * baseSize;
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
