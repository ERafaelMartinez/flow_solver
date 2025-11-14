#!/bin/bash

# Quick script to enter the Docker development environment

echo "🐳 Entering Docker development container..."

# Change to docker directory
cd "$(dirname "$0")"

# Check if container is running
if [ "$(docker-compose ps -q numsim-dev)" ]; then
    # If user passed arguments, execute them in container, otherwise start interactive bash
    if [ $# -eq 0 ]; then
        docker-compose exec numsim-dev /bin/bash
    else
        docker-compose exec numsim-dev "$@"
    fi
else
    echo "⚠️  Container not running. Starting it first..."
    docker-compose up -d
    echo "✅ Container started. Entering now..."
    if [ $# -eq 0 ]; then
        docker-compose exec numsim-dev /bin/bash
    else
        docker-compose exec numsim-dev "$@"
    fi
fi