#!/bin/bash

# Script to build and run tests inside Docker container

echo "🧪 Building and running numsim tests..."
echo "========================================"

# Build tests using the test Makefile
echo "📁 Building tests..."
if cd tests && make all; then
    echo "✅ Tests built successfully!"
    echo ""
    
    # Run FieldVariable tests
    echo "🧪 Running FieldVariable tests..."
    echo "-----------------------------------"
    make test_field_variable
    FIELD_TEST_RESULT=$?
    echo ""
    
    # Run DataExchanger tests
    echo "🧪 Running DataExchanger tests..."
    echo "-----------------------------------"
    make test_data_exchanger
    EXCHANGER_TEST_RESULT=$?
    echo ""
    
    # Summary
    echo "========================================"
    echo "📊 Test Summary:"
    if [ $FIELD_TEST_RESULT -eq 0 ]; then
        echo "  ✅ FieldVariable tests: PASSED"
    else
        echo "  ❌ FieldVariable tests: FAILED"
    fi
    
    if [ $EXCHANGER_TEST_RESULT -eq 0 ]; then
        echo "  ✅ DataExchanger tests: PASSED"
    else
        echo "  ❌ DataExchanger tests: FAILED"
    fi
    echo "========================================"
    
    # Exit with failure if any test failed
    if [ $FIELD_TEST_RESULT -ne 0 ] || [ $EXCHANGER_TEST_RESULT -ne 0 ]; then
        exit 1
    fi
else
    echo "❌ Test build failed"
    exit 1
fi
