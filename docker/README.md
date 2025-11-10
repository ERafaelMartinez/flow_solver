# Docker Development Environment

This folder contains all Docker-related files for the numsim project.

## Quick Start

```bash
# Build and start the development environment
docker/docker-build.sh

# Enter the container
docker/docker-enter.sh

# Build the project (inside container)
./build-in-docker.sh

# Run simulation (inside container)
./numsim parameters.txt

# Cleanup when done (from host)
docker/docker-cleanup.sh
```

## Files

- `Dockerfile` - Ubuntu 22.04 environment with VTK and development tools
- `docker-compose.yml` - Service configuration with volume mounting
- `docker-build.sh` - Build image and start container
- `docker-enter.sh` - Enter running container
- `docker-cleanup.sh` - Stop and cleanup containers

## Features

- ✅ VTK 9.1 pre-installed
- ✅ G++ compiler and CMake
- ✅ Real-time file sharing between host and container
- ✅ Persistent build cache
- ✅ Clean project organization

## WORKFLOW WITH DOCKER

### 🚀 Complete Workflow

#### **Step 1: Build and Start the Docker Environment**
```bash
# From your project root directory (/Users/Rafael/Workspace/2025ws-numsim/)
docker/docker-build.sh
```

This script will:
- Build the Docker image with Ubuntu + VTK + g++
- Start the container with your project files mounted
- Give you instructions for next steps

#### **Step 2: Enter the Container**
```bash
# From your project root (while container is running)
docker/docker-enter.sh
```

You'll now be inside the container with a shell prompt like:
```
root@container:/workspace#
```

#### **Step 3: Build Your Project (Inside Container)**
```bash
# You're now inside the container
./build-in-docker.sh
```

This will:
- Configure the project with CMake
- Compile your numsim executable
- Show you the build results

#### **Step 4: Run Your Simulation (Inside Container)**
```bash
# After successful build, still inside container
cd build
./numsim ../parameters.txt
```

This runs your fluid dynamics simulation with the parameters from `parameters.txt`.

#### **Step 5: View Results (On Your Mac)**
```bash
# Exit the container (press Ctrl+D or type 'exit')
exit

# Back on your Mac, check the output files
ls -la build/
# Your simulation results will be here, ready for visualization
```

#### **Step 6: Cleanup When Done**
```bash
# From your project root on Mac
docker/docker-cleanup.sh
```

### 🔄 Daily Development Cycle

Once you have everything set up, your typical workflow becomes:

```bash
# 1. Start development session
docker/docker-enter.sh

# 2. Build and test (inside container)
./build-in-docker.sh && cd build && ./numsim ../parameters.txt

# 3. Exit container when done
exit

# 4. View/analyze results on Mac using your preferred tools
```

### 📋 Quick Reference Commands

| Action | Command |
|--------|---------|
| **First time setup** | `docker/docker-build.sh` |
| **Enter container** | `docker/docker-enter.sh` |
| **Build project** | `./build-in-docker.sh` (inside container) |
| **Run simulation** | `./numsim ../parameters.txt` (inside container) |
| **One-liner build+run** | `docker/docker-enter.sh -c "./build-in-docker.sh && cd build && ./numsim ../parameters.txt"` |
| **Cleanup everything** | `docker/docker-cleanup.sh` |

### 🎯 Key Benefits of This Workflow

✅ **Edit on Mac**: Use VS Code, your favorite editor  
✅ **Build in Linux**: VTK and all dependencies work perfectly  
✅ **Results on Mac**: Output files immediately available for analysis  
✅ **Clean separation**: Docker complexity hidden in `docker/` folder  
✅ **No data loss**: Your files always stay on your Mac  

The beauty of this setup is that after the initial `docker/docker-build.sh`, you mainly just use `docker/docker-enter.sh` to jump in and out of your development environment!