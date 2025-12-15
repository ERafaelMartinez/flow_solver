# Makefile for numsim tests
# This Makefile compiles and runs unit tests for the numsim project

CXX = mpic++
CXXFLAGS = -std=c++14 -Wall -Wextra -I../src
LDFLAGS = 

# Source directories
SRC_DIR = ../src
STORAGE_DIR = $(SRC_DIR)/storage
COMM_DIR = $(SRC_DIR)/communication
PART_DIR = $(SRC_DIR)/partitioning

# Build directory
BUILD_DIR = ../build/tests

# Object files needed for tests
STORAGE_OBJS = $(BUILD_DIR)/array2d.o $(BUILD_DIR)/field_variable.o
COMM_OBJS = $(BUILD_DIR)/exchanger.o
PART_OBJS = $(BUILD_DIR)/partitioning.o

# Test executables
TEST_FIELD_VARIABLE = $(BUILD_DIR)/test_field_variable
TEST_DATA_EXCHANGER = $(BUILD_DIR)/test_data_exchanger

.PHONY: all clean test test_field_variable test_data_exchanger

all: $(TEST_FIELD_VARIABLE) $(TEST_DATA_EXCHANGER)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile storage objects
$(BUILD_DIR)/array2d.o: $(STORAGE_DIR)/array2d.cpp $(STORAGE_DIR)/array2d.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/field_variable.o: $(STORAGE_DIR)/field_variable.cpp $(STORAGE_DIR)/field_variable.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile communication objects
$(BUILD_DIR)/exchanger.o: $(COMM_DIR)/exchanger.cpp $(COMM_DIR)/exchanger.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile partitioning objects
$(BUILD_DIR)/partitioning.o: $(PART_DIR)/partitioning.cpp $(PART_DIR)/partitioning.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build test_field_variable (no MPI needed)
$(TEST_FIELD_VARIABLE): test_field_variable.cpp $(STORAGE_OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Build test_data_exchanger (MPI needed)
$(TEST_DATA_EXCHANGER): test_data_exchanger.cpp $(STORAGE_OBJS) $(COMM_OBJS) $(PART_OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Run all tests
test: test_field_variable test_data_exchanger

# Run field variable tests (serial)
test_field_variable: $(TEST_FIELD_VARIABLE)
	@echo "=== Running FieldVariable Tests ==="
	@$(TEST_FIELD_VARIABLE)
	@echo ""

# Run data exchanger tests (parallel with 4 processes)
test_data_exchanger: $(TEST_DATA_EXCHANGER)
	@echo "=== Running DataExchanger Tests (4 processes) ==="
	@mpirun -np 4 $(TEST_DATA_EXCHANGER)
	@echo ""

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)

# Help target
help:
	@echo "Available targets:"
	@echo "  all                  - Build all tests"
	@echo "  test                 - Run all tests"
	@echo "  test_field_variable  - Build and run FieldVariable tests"
	@echo "  test_data_exchanger  - Build and run DataExchanger tests (requires 4 MPI processes)"
	@echo "  clean                - Remove build artifacts"
	@echo "  help                 - Show this help message"
