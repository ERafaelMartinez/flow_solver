#!/bin/bash

# Script to clean up, rebuild, and relaunch the Docker image

echo "🧹 Cleaning up existing Docker containers and images..."
./docker-cleanup.sh

if [ $? -ne 0 ]; then
    echo "❌ Cleanup failed. Aborting."
    exit 1
fi

echo "🔨 Rebuilding the Docker image..."
./docker-build.sh

if [ $? -ne 0 ]; then
    echo "❌ Build failed. Aborting."
    exit 1
fi

echo "🚀 Relaunching the Docker container..."
./docker-enter.sh

if [ $? -eq 0 ]; then
    echo "✅ Docker container is ready!"
else
    echo "❌ Failed to enter the Docker container."
    exit 1
fi