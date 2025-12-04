FROM ubuntu:22.04

# Set non-interactive mode to avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install essential development tools
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    cmake \
    git \
    wget \
    curl \
    libvtk9-dev \
    libvtk9.1 \
    vtk9 \
    python3-vtk9 \
    pkg-config \
    qtbase5-dev \
    qttools5-dev \
    qttools5-dev-tools \
    libqt5opengl5-dev \
    libqt5widgets5 \
    libqt5core5a \
    libqt5gui5 \
    libopengl-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    libgl1-mesa-dev \
    libgl1-mesa-glx \
    libglew-dev \
    openjdk-11-jdk \
    mpich \
    && rm -rf /var/lib/apt/lists/*

# Ensure libjawt.so is correctly linked
RUN echo "/usr/lib/jvm/java-11-openjdk-amd64/lib" > /etc/ld.so.conf.d/java.conf \
    && ldconfig

# Ensure libjvm.so is correctly linked
RUN echo "/usr/lib/jvm/java-11-openjdk-amd64/lib/server" >> /etc/ld.so.conf.d/java.conf \
    && ldconfig

# Set environment variables for VTK
ENV VTK_DIR=/usr/lib/x86_64-linux-gnu/cmake/vtk-9.1

# Create a working directory for your project
WORKDIR /workspace

# Create a build directory (files will be mounted via volume)
RUN mkdir -p build

# Default command
CMD ["/bin/bash"]
