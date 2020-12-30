# `gfold`

A C++ implementation of the G-FOLD algorithm.

G-FOLD is an acronym for Guidance for Fuel-Optimal Large
Diverts. It is a guidance algorithm that attempts to find
the best trajectory for powered descent, such as for rocket
boosters and planetary landers. It uses convex optimization
to get as close to the landing target as possible, and then
optimizes for fuel usage to select the most optimal
trajectory.

This project attempts to implement the G-FOLD algorithm in
C++ as a demo that plots the some of the numerical examples
in the original paper using matplotlib-cpp.

# Demo

![gfold.png](https://i.postimg.cc/VvbrrRYD/gfold.png)

# Building

#### Prerequisites

  * CMake
  * Eigen  
  * Catch2
  * Python3 Matplotlib

#### Build

``` shell
git clone https://github.com/AgentTroll/gfold.git
cd gfold
mkdir build && cd build
cmake .. && make
./gfold
```

# Credits

Built with [CLion](https://www.jetbrains.com/clion/)

Utilizes:

  * [Eigen](https://eigen.tuxfamily.org/)
  * [Epigraph](https://github.com/EmbersArc/Epigraph)
  * [matplotlib-cpp](https://github.com/lava/matplotlib-cpp)

# References

"Lossless Convexification of Non-Convex Control Bound and
Pointing Constraints of the Soft Landing Optimal Control
Problem ." B. Acikmese, J. M. Carson III, and L. Blackmore.
IEEE Transactions on Control Systems Technology, Volume 21,
Issue 6 (2013).

Link: http://larsblackmore.com/iee_tcst13.pdf
