# VertexCoverSolver

## Requirements

- Development tools (compiler, dependencies, ...)
  
  - Ubuntu/Debian:
    
    ```bash
    sudo apt update
    sudo apt install build-essential
    ```
  - Ubuntu/Debian:
    
    ```bash
    sudo apt install cmake
    ```

## Compilation

In the root directory, run the following commands

Linux Makefiles:

```bash
cmake . -B build -DCMAKE_BUILD_TYPE=RelWithDbInfo 
cmake --build build
```

The compiled executables are located in the `build/bin` subdirectory.

## Run the benchmark script

e.g.

```bash
./../vc-data-students/benchmark-fast.sh "./build/bin/VertexCoverSolverClique"
```

## Choosing other Bounds

In the ArrayGraph.cpp file, find the function ```int ArrayGraph::getLowerBoundVC()```. There, one can select one of four implemented bounds: Clique Bound, Cycle Bound, LP Bound.