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

# execute valgrind in logs folder
cd ..; cmake --build build; cd logs; valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind_out.txt ./../build/bin/VertexCoverSolver < ./../../vc-data-students/1-random/000600_000000003600.dimacs

## SSH remote testing machines

# connect

ssh algeng-ss23-team6@aba01.akt.tu-berlin.de
password: 9b5MEDypeZPf

# upload data to their remote machine

scp -r ./vc-data-students/ algeng-ss23-team6@aba01.akt.tu-berlin.de:./local/
scp ./VertexCoverSolver algeng-ss23-team6@aba01.akt.tu-berlin.de:./local/src

# execute gprof profiling

cmake . -B build -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg

cd ..; cmake --build build; cd logs; ../build/bin/VertexCoverSolver < ./../../vc-data-students/1-random/000600_000000003600.dimacs

gprof ../build/bin/VertexCoverSolver gmon.out > results.txt
