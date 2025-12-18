#include "../src/storage/field_variable.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// Test helper function
bool doubleEquals(double a, double b, double epsilon = 1e-9) {
  return std::abs(a - b) < epsilon;
}

// Test 1: Basic getColumnValues and setColumnValues with range starting at 0
bool testColumnValuesRangeZero() {
  std::cout << "Test 1: Column values with range [0, height)..." << std::endl;

  std::array<int, 2> gridSize = {5, 7}; // 5 columns, 7 rows
  std::array<double, 2> offset = {0.0, 0.0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable field(gridSize, offset, cellSize);

  // Initialize field with known values
  for (int j = 0; j < gridSize[1]; ++j) {
    for (int i = 0; i < gridSize[0]; ++i) {
      field.at(i, j) = i * 100.0 + j;
    }
  }

  // Get column 2 values
  std::array<int, 2> rowRange = {0, 7};
  std::vector<double> columnValues = field.getColumnValues(2, rowRange);

  // Verify we got the right values
  if (columnValues.size() != 7) {
    std::cout << "  FAILED: Expected 7 values, got " << columnValues.size()
              << std::endl;
    return false;
  }

  for (int j = 0; j < 7; ++j) {
    double expected = 2 * 100.0 + j; // column 2, row j
    if (!doubleEquals(columnValues[j], expected)) {
      std::cout << "  FAILED: columnValues[" << j << "] = " << columnValues[j]
                << ", expected " << expected << std::endl;
      return false;
    }
  }

  // Now set column 3 with new values
  std::vector<double> newValues = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
  field.setColumnValues(3, rowRange, newValues);

  // Verify the values were set correctly
  for (int j = 0; j < 7; ++j) {
    if (!doubleEquals(field.at(3, j), newValues[j])) {
      std::cout << "  FAILED: field.at(3, " << j << ") = " << field.at(3, j)
                << ", expected " << newValues[j] << std::endl;
      return false;
    }
  }

  std::cout << "  PASSED" << std::endl;
  return true;
}

// Test 2: setColumnValues with range starting at non-zero offset (THE CRITICAL
// TEST)
bool testColumnValuesRangeOffset() {
  std::cout << "Test 2: Column values with range [2, 5)..." << std::endl;

  std::array<int, 2> gridSize = {5, 7};
  std::array<double, 2> offset = {0.0, 0.0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable field(gridSize, offset, cellSize);
  field.setToZero();

  // Set column 1, rows 2-4 (range [2, 5))
  std::array<int, 2> rowRange = {2, 5};
  std::vector<double> values = {100.0, 200.0,
                                300.0}; // 3 values for rows 2, 3, 4

  field.setColumnValues(1, rowRange, values);

  // Verify: rows 0, 1 should be 0, rows 2-4 should have our values, rows 5-6
  // should be 0
  for (int j = 0; j < 7; ++j) {
    double expected;
    if (j < 2 || j >= 5) {
      expected = 0.0;
    } else {
      expected = values[j - 2]; // This is the critical indexing!
    }

    if (!doubleEquals(field.at(1, j), expected)) {
      std::cout << "  FAILED: field.at(1, " << j << ") = " << field.at(1, j)
                << ", expected " << expected << std::endl;
      return false;
    }
  }

  std::cout << "  PASSED" << std::endl;
  return true;
}

// Test 3: Basic getRowValues and setRowValues with range starting at 0
bool testRowValuesRangeZero() {
  std::cout << "Test 3: Row values with range [0, width)..." << std::endl;

  std::array<int, 2> gridSize = {8, 5}; // 8 columns, 5 rows
  std::array<double, 2> offset = {0.0, 0.0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable field(gridSize, offset, cellSize);

  // Initialize field
  for (int j = 0; j < gridSize[1]; ++j) {
    for (int i = 0; i < gridSize[0]; ++i) {
      field.at(i, j) = i * 100.0 + j;
    }
  }

  // Get row 3 values
  std::array<int, 2> columnRange = {0, 8};
  std::vector<double> rowValues = field.getRowValues(3, columnRange);

  // Verify
  if (rowValues.size() != 8) {
    std::cout << "  FAILED: Expected 8 values, got " << rowValues.size()
              << std::endl;
    return false;
  }

  for (int i = 0; i < 8; ++i) {
    double expected = i * 100.0 + 3; // column i, row 3
    if (!doubleEquals(rowValues[i], expected)) {
      std::cout << "  FAILED: rowValues[" << i << "] = " << rowValues[i]
                << ", expected " << expected << std::endl;
      return false;
    }
  }

  // Set row 2 with new values
  std::vector<double> newValues = {11.0, 22.0, 33.0, 44.0,
                                   55.0, 66.0, 77.0, 88.0};
  field.setRowValues(2, columnRange, newValues);

  // Verify
  for (int i = 0; i < 8; ++i) {
    if (!doubleEquals(field.at(i, 2), newValues[i])) {
      std::cout << "  FAILED: field.at(" << i << ", 2) = " << field.at(i, 2)
                << ", expected " << newValues[i] << std::endl;
      return false;
    }
  }

  std::cout << "  PASSED" << std::endl;
  return true;
}

// Test 4: setRowValues with range starting at non-zero offset (THE CRITICAL
// TEST)
bool testRowValuesRangeOffset() {
  std::cout << "Test 4: Row values with range [3, 6)..." << std::endl;

  std::array<int, 2> gridSize = {8, 5};
  std::array<double, 2> offset = {0.0, 0.0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable field(gridSize, offset, cellSize);
  field.setToZero();

  // Set row 1, columns 3-5 (range [3, 6))
  std::array<int, 2> columnRange = {3, 6};
  std::vector<double> values = {111.0, 222.0,
                                333.0}; // 3 values for columns 3, 4, 5

  field.setRowValues(1, columnRange, values);

  // Verify: columns 0-2 should be 0, columns 3-5 should have our values,
  // columns 6-7 should be 0
  for (int i = 0; i < 8; ++i) {
    double expected;
    if (i < 3 || i >= 6) {
      expected = 0.0;
    } else {
      expected = values[i - 3]; // This is the critical indexing!
    }

    if (!doubleEquals(field.at(i, 1), expected)) {
      std::cout << "  FAILED: field.at(" << i << ", 1) = " << field.at(i, 1)
                << ", expected " << expected << std::endl;
      return false;
    }
  }

  std::cout << "  PASSED" << std::endl;
  return true;
}

// Test 5: Simulate ghost cell exchange pattern
bool testGhostCellPattern() {
  std::cout << "Test 5: Ghost cell exchange pattern..." << std::endl;

  // Simulate a 5x5 internal grid with ghost cells (7x7 total)
  std::array<int, 2> gridSize = {7, 7};
  std::array<double, 2> offset = {0.0, 0.0};
  std::array<double, 2> cellSize = {1.0, 1.0};

  FieldVariable field(gridSize, offset, cellSize);
  field.setToZero();

  // Fill internal cells (indices 1-5) with a pattern
  for (int j = 1; j <= 5; ++j) {
    for (int i = 1; i <= 5; ++i) {
      field.at(i, j) = i * 10.0 + j;
    }
  }

  // Extract left internal boundary (column 1)
  std::array<int, 2> rowRange = {0, 7};
  std::vector<double> leftBoundary = field.getColumnValues(1, rowRange);

  // Simulate receiving this as right ghost cells from left neighbor
  // Set it to column 0 (left ghost)
  field.setColumnValues(0, rowRange, leftBoundary);

  // Verify column 0 now has the same values as column 1
  for (int j = 0; j < 7; ++j) {
    if (!doubleEquals(field.at(0, j), field.at(1, j))) {
      std::cout << "  FAILED: Ghost cell mismatch at row " << j << std::endl;
      return false;
    }
  }

  std::cout << "  PASSED" << std::endl;
  return true;
}

int main() {
  std::cout << "=== FieldVariable Unit Tests ===" << std::endl << std::endl;

  bool allPassed = true;

  allPassed &= testColumnValuesRangeZero();
  allPassed &= testColumnValuesRangeOffset();
  allPassed &= testRowValuesRangeZero();
  allPassed &= testRowValuesRangeOffset();
  allPassed &= testGhostCellPattern();

  std::cout << std::endl;
  if (allPassed) {
    std::cout << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
  } else {
    std::cout << "=== SOME TESTS FAILED ===" << std::endl;
    return 1;
  }
}
