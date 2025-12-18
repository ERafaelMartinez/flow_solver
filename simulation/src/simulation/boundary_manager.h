#pragma once

#include "../discretization/discretization.h"
#include "../partitioning/partitioning.h"
#include "../settings/settings.h"
#include <memory>

/**
 * @brief Manages boundary conditions for velocity and intermediate velocity
 * fields.
 *
 * This class encapsulates all boundary condition logic for the simulation,
 * providing fine-grained control over individual boundaries which is essential
 * for parallel domain decomposition.
 */
class BoundaryManager {
public:
  /**
   * @brief Construct a new Boundary Manager
   *
   * @param discretization Shared pointer to the discretization
   * @param partitioning Shared pointer to the partitioning
   * @param settings Pointer to simulation settings
   */
  BoundaryManager(std::shared_ptr<Discretization> discretization,
                  std::shared_ptr<Partitioning> partitioning,
                  Settings *settings);

  /**
   * @brief Apply boundary conditions to velocity fields u and v
   */
  void setBoundaryConditionsVelocity();

  /**
   * @brief Apply boundary conditions to intermediate velocity fields F and G
   */
  void setBoundaryConditionsFG();

private:
  std::shared_ptr<Discretization> discretization_;
  std::shared_ptr<Partitioning> partitioning_;
  Settings *settings_;

  // Helper methods for applying boundary conditions to u and v
  // Split by boundary side and velocity component for parallel control
  void applyTopBoundaryU();
  void applyBottomBoundaryU();
  void applyTopBoundaryV();
  void applyBottomBoundaryV();
  void applyLeftBoundaryU();
  void applyRightBoundaryU();
  void applyLeftBoundaryV();
  void applyRightBoundaryV();

  // Helper methods for applying boundary conditions to F and G
  void applyTopBoundaryF();
  void applyBottomBoundaryF();
  void applyTopBoundaryG();
  void applyBottomBoundaryG();
  void applyLeftBoundaryF();
  void applyRightBoundaryF();
  void applyLeftBoundaryG();
  void applyRightBoundaryG();
};
