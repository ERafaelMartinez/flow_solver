# Numsim WS2025

## How to start

1. Make sure you have the following tools installed:

- clang or g++
- make
- cmake
- vtk

2. Build the project for the first time

execute the following command to make a clean (debug) build of the project in the root dir

```sh
make pristine
```

> **note**: you need to always call `make pristine` if any of the CMakeLists.txt files are changed
> **note**: you can also use `make pristine_release` for a release build


3. re-build the project after changes using the following command in the root dir

```sh
make
```

> **note**: you can also use `make release` for a release build

4. Run the simulation by executing the binary in the build folder

```sh
cd build/
./numsim
```

5. Open the output using Paraview

You can then open the output files in `build/out/` using ParaView
