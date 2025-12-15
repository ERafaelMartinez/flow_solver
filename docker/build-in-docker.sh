#!/bin/bash

# Script to build the numsim project inside Docker container

echo "🔨 Building numsim project..."
echo "================================"

# Build using the project's Makefile
echo "📁 Using project Makefile for build..."
if make debug; then
    echo "✅ Build successful!"
    echo ""
    echo "🎉 numsim executable created successfully!"
    echo "📋 To run the simulation:"
    echo "   cd build && ./numsim ../parameters.txt"
    echo ""
    echo "📁 Build directory contents:"
    ls -la build/
else
    echo "❌ Build failed"
    exit 1
fi