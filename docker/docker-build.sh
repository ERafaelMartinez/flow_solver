#!/bin/bash

# Build the Docker image and start the development environment
echo "🐳 Building Docker development environment..."

# Change to docker directory
echo "$(dirname "$0")"
cd "$(dirname "$0")"

docker-compose build

# Start the container
echo "🚀 Starting development container..."
docker-compose up -d

# Check if container started successfully
if [ $? -eq 0 ]; then
    echo "✅ Container started successfully!"
    echo "📝 You can now enter the container with:"
    echo "   docker/docker-enter.sh"
    echo ""
    echo "🔧 Or run the build script inside:"
    echo "   docker/docker-enter.sh -c './build-in-docker.sh'"
    echo ""
    echo "📁 All Docker files are now in the docker/ folder"
else
    echo "❌ Failed to start container"
    exit 1
fi