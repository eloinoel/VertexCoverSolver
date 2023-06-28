# VertexCoverSolver

## Requirements

- Development tools (compiler, dependencies, ...)

  - Ubuntu/Debian:
    ```bash
    sudo apt update
    sudo apt install build-essential
    ```

- CMake >= 3.14

  - Ubuntu/Debian:
    ```bash
    sudo apt install cmake
    ```

## Compilation

cmake . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# from the logs folder, run

./../../vc-data-students/benchmark-fast.sh "./../build/bin/VertexCoverSolver"

# test one case

build/bin/VertexCoverSolver < ../vc-data-students/1-random/000002_000000000012.dimacs

# rebuild and test one case from logs folder

cd ..; cmake --build build; cd logs; ../build/bin/VertexCoverSolver < ./../../vc-data-students/1-random/000002_000000000012.dimacs

# rebuild and execute benchmark in logs folder

cd ..; cmake --build build; cd logs; ./../../vc-data-students/benchmark-fast.sh "./../build/bin/VertexCoverSolver" 30 3

## SSH remote testing machines

# connect

ssh algeng-ss23-team6@aba01.akt.tu-berlin.de
password: 9b5MEDypeZPf

# upload data to their remote machine

scp -r ./vc-data-students/ algeng-ss23-team6@aba01.akt.tu-berlin.de:./local/
scp ./VertexCoverSolver algeng-ss23-team6@aba01.akt.tu-berlin.de:./local/src

# execute gprof profiling

cmake -D DCMAKE_CXX_FLAGS="-pg" DCMAKE_EXE_LINKER_FLAGS="-pg" DCMAKE_SHARED_LINKER_FLAGS="-pg" . -B build

cd ..; cmake --build build; cd logs; ../build/bin/VertexCoverSolver < ./../../vc-data-students/1-random/000600_000000003600.dimacs

gprof ../build/bin/VertexCoverSolver gmon.out > results.txt

# gprof call graph
1. copy results.txt into gprof2dot.py directory
2. python3 gprof2dot.py -w results.txt > call_graph.dot
view dot file in viewer --> xdot call_graph.dot
convert dot file to png --> dot -Tpng call_graph.dot -o call_graph.png
