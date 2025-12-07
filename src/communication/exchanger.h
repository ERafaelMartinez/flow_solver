#pragma once

#include "../partitioning/partitioning.h"
#include "../staggered_grid/staggered_grid.h"
#include "../storage/field_variable.h"
#include <array>
#include <assert.h>
#include <memory>
#include <mpi.h>
#include <vector>

class DataBuffer {
public:
  // constructor
  explicit DataBuffer(int nRegs) : nRegs_(nRegs) {
    // initialize data buffer internal array with zeros
    data_.resize(nRegs_, 0.0);
  }
  // destructor
  ~DataBuffer() = default;

  // Accessor for the internal array
  std::vector<double> &asArray() { return data_; }

  const std::vector<double> &asArray() const { return data_; }

  // Updator data buffer with new data
  void update(const std::vector<double> &newData) {
    // assert(newData.size() == nRegs_);
    data_ = newData;
  }

private:
  int nRegs_;
  std::vector<double> data_;
};

class DataExchanger {
public:
  // constructor
  DataExchanger(std::shared_ptr<Partitioning> partitioning);

  // exchange variables via MPI with neighbors
  void exchange(FieldVariable &fieldVar);
  
  // obtain the maximum velocity {maxU, maxV} from all ranks
  double getMaximumVelocity(std::array<double, 2> &maxVelocity);

  // obtain the maximum residual from all ranks
  double getResidual(double res);

private:
  std::shared_ptr<Partitioning> partitioning_;
  std::array<int, 4> neighborsRank_;

  // declare structured accessors for field variables and
  // their boundary rules
  std::array<std::array<int, 4>, 6> boundariesRules_;
  std::array<FieldVariable *, 6> fieldVars_;

  // provide new data buffers, one for each ghost boundary based on the
  // field variable's size
  std::array<std::shared_ptr<DataBuffer>, 4>
  getDataBuffers_(FieldVariable &fieldVar);

  // pack data from the field variables into the buffers
  void packDataBuffers_(std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers,
                        FieldVariable &fieldVar);

  // unpack data from the buffers into the field variables
  void
  unpackDataBuffers_(std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers,
                     FieldVariable &fieldVar);

  // receive data from neighbors and update buffer state
  void
  receiveDataBuffers_(std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers);

  // send data to neighbors from current buffer state
  void
  sendDataBuffers_(std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers);

  // declare MPI send/receive requests for non-blocking
  // communication
  std::vector<MPI_Request> sendRequests_;
  std::vector<MPI_Request> receiveRequests_;
};
