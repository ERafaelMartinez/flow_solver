#!/bin/bash

# Cleanup script for Docker environment

echo "🧹 Cleaning up Docker environment..."

# Change to docker directory
cd "$(dirname "$0")"

# Stop and remove containers
echo "🛑 Stopping containers..."
docker-compose down

# Optional: Remove the built image (uncomment if you want to rebuild from scratch)
# echo "🗑️  Removing Docker image..."
# docker-compose down --rmi all

# Optional: Remove volumes (uncomment if you want to clean build cache)
# echo "🗑️  Removing volumes..."
# docker-compose down --volumes

echo "✅ Cleanup complete!"
echo "💡 To start fresh, run: docker/docker-build.sh"