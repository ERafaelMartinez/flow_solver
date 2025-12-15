#!/bin/bash

# Strong Scaling Analysis Script
# This script profiles the execution time of numsim_parallel with varying processor counts
# to analyze strong scaling behavior (fixed problem size, varying processors)

# Configuration
BUILD_DIR="./build"
EXECUTABLE="${BUILD_DIR}/numsim_parallel"
PARAMETER_FILE="./parameters/parameters-5.txt"
OUTPUT_FILE="scaling_results_s5.csv"
LOG_DIR="./scaling_logs"

# Processor counts to test (modify as needed)
PROC_COUNTS=(1 2 3 4)

# Number of runs per configuration for averaging
NUM_RUNS=3

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1" >&2
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" >&2
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" >&2
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Function to check if executable exists
check_executable() {
    if [ ! -f "$EXECUTABLE" ]; then
        print_error "Executable not found: $EXECUTABLE"
        print_info "Please build the project first using 'make release' or 'make debug'"
        exit 1
    fi
}

# Function to check if parameter file exists
check_parameter_file() {
    if [ ! -f "$PARAMETER_FILE" ]; then
        print_error "Parameter file not found: $PARAMETER_FILE"
        exit 1
    fi
}

# Function to create log directory
create_log_dir() {
    mkdir -p "$LOG_DIR"
    print_info "Created log directory: $LOG_DIR"
}

# Function to run a single test
run_test() {
    local nprocs=$1
    local run_num=$2
    local log_file="${LOG_DIR}/run_p${nprocs}_r${run_num}.log"
    
    print_info "Running with $nprocs processor(s), run $run_num/$NUM_RUNS..."
    
    # Clean output directory
    rm -rf ${BUILD_DIR}/out
    
    # Run the program and measure time
    local start_time=$(date +%s.%N)
    
    cd "$BUILD_DIR" && mpiexec -n $nprocs ./numsim_parallel ../$PARAMETER_FILE > "../${log_file}" 2>&1
    local exit_code=$?
    
    cd ..
    
    local end_time=$(date +%s.%N)
    local elapsed=$(awk "BEGIN {print $end_time - $start_time}")
    
    if [ $exit_code -ne 0 ]; then
        print_error "Run failed with exit code $exit_code (see $log_file)"
        echo "-1"
        return 1
    fi
    
    echo "$elapsed"
    return 0
}

# Function to calculate statistics
calculate_stats() {
    local times=("$@")
    local count=${#times[@]}
    
    # Use awk to calculate all statistics at once
    local stats=$(printf '%s\n' "${times[@]}" | awk '
    {
        sum += $1
        values[NR] = $1
        if (NR == 1) {
            min = $1
            max = $1
        } else {
            if ($1 < min) min = $1
            if ($1 > max) max = $1
        }
    }
    END {
        mean = sum / NR
        for (i = 1; i <= NR; i++) {
            variance_sum += (values[i] - mean) ^ 2
        }
        stddev = sqrt(variance_sum / NR)
        printf "%.6f %.6f %.6f %.6f", mean, stddev, min, max
    }
    ')
    
    echo "$stats"
}

# Main execution
main() {
    print_info "=== Strong Scaling Analysis ==="
    print_info "Executable: $EXECUTABLE"
    print_info "Parameter file: $PARAMETER_FILE"
    print_info "Processor counts: ${PROC_COUNTS[*]}"
    print_info "Runs per configuration: $NUM_RUNS"
    echo ""
    
    # Checks
    check_executable
    check_parameter_file
    create_log_dir
    
    # Initialize results file
    echo "NumProcs,MeanTime(s),StdDev(s),MinTime(s),MaxTime(s),Speedup,Efficiency(%)" > "$OUTPUT_FILE"
    
    # Store baseline time (serial execution)
    local baseline_time=0
    
    # Run tests for each processor count
    for nprocs in "${PROC_COUNTS[@]}"; do
        print_info "========================================="
        print_info "Testing with $nprocs processor(s)"
        print_info "========================================="
        
        local times=()
        local failed=0
        
        # Run multiple times for statistics
        for run in $(seq 1 $NUM_RUNS); do
            local elapsed=$(run_test $nprocs $run)
            
            if [ "$elapsed" == "-1" ]; then
                failed=1
                break
            fi
            
            times+=($elapsed)
            print_success "Run $run completed in ${elapsed}s"
        done
        
        if [ $failed -eq 1 ]; then
            print_error "Skipping processor count $nprocs due to failures"
            continue
        fi
        
        # Calculate statistics
        read mean stddev min max <<< $(calculate_stats "${times[@]}")
        
        # Calculate speedup and efficiency
        if [ $nprocs -eq 1 ] || [ $nprocs -eq ${PROC_COUNTS[0]} ]; then
            baseline_time=$mean
            speedup=1.0000
            efficiency=100.00
        else
            speedup=$(awk "BEGIN {printf \"%.4f\", $baseline_time / $mean}")
            efficiency=$(awk "BEGIN {printf \"%.2f\", ($speedup / $nprocs) * 100}")
        fi
        
        # Write results
        echo "$nprocs,$mean,$stddev,$min,$max,$speedup,$efficiency" >> "$OUTPUT_FILE"
        
        print_success "Results for $nprocs processor(s):"
        print_info "  Mean time: ${mean}s"
        print_info "  Std dev:   ${stddev}s"
        print_info "  Min time:  ${min}s"
        print_info "  Max time:  ${max}s"
        print_info "  Speedup:   ${speedup}x"
        print_info "  Efficiency: ${efficiency}%"
        echo ""
    done
    
    print_success "=== Scaling Analysis Complete ==="
    print_info "Results saved to: $OUTPUT_FILE"
    print_info "Logs saved to: $LOG_DIR/"
    echo ""
    
    # Display summary table using awk
    print_info "Summary Table:"
    awk -F',' '
    NR == 1 {
        printf "%-10s %-15s %-15s %-15s %-15s %-10s %-15s\n", $1, $2, $3, $4, $5, $6, $7
        printf "%-10s %-15s %-15s %-15s %-15s %-10s %-15s\n", "----------", "---------------", "---------------", "---------------", "---------------", "----------", "---------------"
    }
    NR > 1 {
        printf "%-10s %-15s %-15s %-15s %-15s %-10s %-15s\n", $1, $2, $3, $4, $5, $6, $7
    }
    ' "$OUTPUT_FILE"
}

# Run main function
main
