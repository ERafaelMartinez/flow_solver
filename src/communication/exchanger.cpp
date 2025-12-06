#include "exchanger.h"
#include "../partitioning/partitioning.h"
#include "../staggered_grid/staggered_grid.h"
#include "../storage/field_variable.h"
#include <array>

DataExchanger::DataExchanger(std::shared_ptr<Partitioning> partitioning) {
  partitioning_ = partitioning;
  neighborsRank_[0] = partitioning_->leftNeighbourRankNo();
  neighborsRank_[1] = partitioning_->rightNeighbourRankNo();
  neighborsRank_[2] = partitioning_->bottomNeighbourRankNo();
  neighborsRank_[3] = partitioning_->topNeighbourRankNo();
}

// get data buffers for the field variable: {left, right, bottom, top}
std::array<std::shared_ptr<DataBuffer>, 4>
DataExchanger::getDataBuffers_(FieldVariable &fieldVar) {
  std::array<std::shared_ptr<DataBuffer>, 4> dataBuffers;
  for (int i = 0; i < 4; ++i) {
    // allocate the buffers only for the case where neighbors exist
    if (neighborsRank_[i] != -1) {
      // determine if the boundary is vertical or horizontal
      if (i == 0 || i == 1) {
        dataBuffers[i] = std::make_shared<DataBuffer>(fieldVar.size()[1]);
      } else {
        dataBuffers[i] = std::make_shared<DataBuffer>(fieldVar.size()[0]);
      }
    } else {
      // in other case, set dummy buffer of size zero
      dataBuffers[i] = std::make_shared<DataBuffer>(0);
    }
  }
  return dataBuffers;
}

// pack data from the field variables into the buffers
void DataExchanger::packDataBuffers_(
    std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers,
    FieldVariable &fieldVar) {
  // internal bounds are the ones being sent to neighbors
  // assuming all staggered grids have ghost cells on all
  // boundaries, the following hard coded indices are
  // correct (look into staggered grid: boundaries as {1, 1, 1, 1})
  std::array<int, 2> rowRange = {0, fieldVar.size()[1]};
  std::array<int, 2> columnRange = {0, fieldVar.size()[0]};
  for (int i = 0; i < 4; ++i) {
    if (neighborsRank_[i] != -1) {
      switch (i) {
      case 0:
        dataBuffers[i]->update(
          // left internal-boundary column index is 1
          fieldVar.getColumnValues(1, rowRange)
        );
        break;
      case 1:
        dataBuffers[i]->update(
          // right internal-boundary column index is size[0] - 2
          fieldVar.getColumnValues(fieldVar.size()[0] - 2, rowRange)
        );
        break;
      case 2:
        dataBuffers[i]->update(
          // top internal-boundary row index is size[1] - 2
          fieldVar.getRowValues(fieldVar.size()[1] - 2, columnRange)
        );
        break;
      case 3:
        dataBuffers[i]->update(
          // bottom internal-boundary row index is 1
          fieldVar.getRowValues(1, columnRange)
        );
        break;
      }
    }
  }
}

// unpack data from the buffers into the field variables
void DataExchanger::unpackDataBuffers_(
    std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers,
    FieldVariable &fieldVar) {
  // ghost cell boundaries are the ones being received from neighbors
  // assuming all staggered grids have ghost cells on all
  // boundaries, the following hard coded indices are
  // correct
  std::array<int, 2> rowRange = {0, fieldVar.size()[1]};
  std::array<int, 2> columnRange = {0, fieldVar.size()[0]};
  for (int i = 0; i < 4; ++i) {
    if (neighborsRank_[i] != -1) {
      switch (i) {
      case 0:
        fieldVar.setColumnValues(
          // left ghost-column index is 0
          0, rowRange, dataBuffers[i]->asArray()
        );
        break;
      case 1:
        fieldVar.setColumnValues(
          // right ghost-column index is size[0] - 1
          fieldVar.size()[0] - 1, rowRange, dataBuffers[i]->asArray()
        );
        break;
      case 2:
        fieldVar.setRowValues(
          // top ghost-row index is size[1] - 1
          fieldVar.size()[1] - 1, columnRange, dataBuffers[i]->asArray()
        );
        break;
      case 3:
        fieldVar.setRowValues(
          // bottom ghost-row index is 0
          0, columnRange, dataBuffers[i]->asArray()
        );
        break;
      }
    }
  }
}

// send data to neighbors processes and store it in the buffers
void DataExchanger::sendDataBuffers_(
    std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers) {
  // declare a potential send request per boundary
  std::array<MPI_Request, 4> sendRequests;

  // iterate over all neighbor boundaries
  for (int i = 0; i < 4; ++i) {
    // ensure there is a neighbor process
    if (neighborsRank_[i] != -1) {
      // instantiate the send request
      sendRequests[i] = MPI_Request();
      // send data to neighbor process
      MPI_Isend(dataBuffers[i]->asArray().data(),
                dataBuffers[i]->asArray().size(), MPI_DOUBLE, neighborsRank_[i],
                0, MPI_COMM_WORLD, &sendRequests[i]);
      // add request to request pool in exchanger
      sendRequests_.push_back(sendRequests[i]);
    }
  }
}

// receive data from neighbors processes and store it in the buffers
void DataExchanger::receiveDataBuffers_(
    std::array<std::shared_ptr<DataBuffer>, 4> &dataBuffers) {
  // declare a potential receive request per boundary
  std::array<MPI_Request, 4> receiveRequests;

  // iterate over all neighbor boundaries
  for (int i = 0; i < 4; ++i) {
    // ensure there is a neighbor process
    if (neighborsRank_[i] != -1) {
      // instantiate the receive request
      receiveRequests[i] = MPI_Request();
      // receive data from neighbor process
      MPI_Irecv(dataBuffers[i]->asArray().data(),
                dataBuffers[i]->asArray().size(), MPI_DOUBLE, neighborsRank_[i],
                0, MPI_COMM_WORLD, &receiveRequests[i]);
      // add request to request pool in exchanger
      receiveRequests_.push_back(receiveRequests[i]);
    }
  }
}

void DataExchanger::exchange(FieldVariable &fieldVar) {
  // clear request pools
  sendRequests_.clear();
  receiveRequests_.clear();

  // get new data buffers
  std::array<std::shared_ptr<DataBuffer>, 4> sendBuffers =
      getDataBuffers_(fieldVar);
  std::array<std::shared_ptr<DataBuffer>, 4> receiveBuffers =
      getDataBuffers_(fieldVar);

  // pack data from the field variables into the buffers
  packDataBuffers_(sendBuffers, fieldVar);
  // send data to neighbors
  sendDataBuffers_(sendBuffers);
  // wait for all sends to finish
  MPI_Waitall(sendRequests_.size(), sendRequests_.data(), MPI_STATUSES_IGNORE);

  // receive data from neighbors
  receiveDataBuffers_(receiveBuffers);
  // wait for all receives to finish
  MPI_Waitall(receiveRequests_.size(), receiveRequests_.data(),
              MPI_STATUSES_IGNORE);
  // unpack data from the buffers into the field variables
  unpackDataBuffers_(receiveBuffers, fieldVar);
}

double DataExchanger::getMaximumTimeStepSize(double &timeStepSize) {
  // clear request pools
  sendRequests_.clear();
  receiveRequests_.clear();

  // all non-main ranks send their time step size to the main rank
  if (partitioning_->ownRankNo != 0) {
    std::shared_ptr<DataBuffer> sendBuffer = std::make_shared<DataBuffer>(1);
    sendBuffer->update(&timeStepSize);
    
    MPI_Request sendRequest;
    MPI_Isend(
      sendBuffer->asArray().data(), sendBuffer->asArray().size(),
      MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &sendRequest
    );
    sendRequests_.push_back(sendRequest);
  }

  // Then, the main rank receives the time step size from the other ranks
  // thus we need one receive buffer for each rank with non-main rank
  if (partitioning_->ownRankNo == 0) {
    // create receive buffers one per non-main rank
    std::vector< std::shared_ptr<DataBuffer> > receiveBuffers;
    receiveBuffers.resize(partitioning_->nRanks() - 1, std::make_shared<DataBuffer>(1));

    // create receive requests one per non-main rank
    receiveRequests_.resize(partitioning_->nRanks() - 1, MPI_Request());
    for (int i = 1; i < partitioning_->nRanks(); ++i) {
      MPI_Irecv(
        receiveBuffers[i - 1]->asArray().data(),
        receiveBuffers[i - 1]->asArray().size(),
        MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
        &receiveRequests_[i - 1]
      );
    }
    // wait for all receives to finish
    MPI_Waitall(
      receiveRequests_.size(), receiveRequests_.data(), MPI_STATUSES_IGNORE
    );
    // extract maximum time step size from the buffers into the time step size
    double maxTimeStepSize = timeStepSize;
    for (int i = 1; i < partitioning_->nRanks(); ++i) {
      maxTimeStepSize = std::max(maxTimeStepSize, receiveBuffers[i]->asArray()[0]);
    }

    // communicate maximum time step size to all ranks
    MPI_Bcast(&maxTimeStepSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // wait for broadcast to finish
    MPI_Waitall(
      receiveRequests_.size(), receiveRequests_.data(), MPI_STATUSES_IGNORE
    );

    return maxTimeStepSize;
  } else {
    // the other ranks read the maximum time step size from the main rank
    double maxTimeStepSize;
    MPI_Bcast(&maxTimeStepSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return maxTimeStepSize;
  }

}
