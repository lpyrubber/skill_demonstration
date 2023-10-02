# skill_demonstration
This repository is a hodgepodge of my previous works in Compuational Physics.

* Particle in Cell method (1-D plasma simulation for demonstrating Landu damping) 

* Finite Element method (Material mechanics)

* Mesh Generator (Delauney triangle)

* Conjegate Gradient (Linear Equation Solver)

* Shallow Water Equation (1-D shallow water shock wave with finite volume method)
  * AVX, MPI, OpenMP, OpenCL, CUDA 
  * MPI + CUDA for multi GPU case
  * OpenMP + AVX for further speed up (thread level + instruction level parallelism)
  * OpenMP + MPI

* Euler Equation (1-D shock wave with finite volume method + AUSM scheme)
  * AVX, MPI, OpenMP, OpenCL, CUDA
  
* validation platform
  * 13 inch Macbook Pro (2020 version)
  * Self-constructed PC
    * CPU: AMD Ryzen 3700x (8 physical cores, 16 cores with hyper threads)
    * DRAM: DDR4 16 GB 
    * GPU: Nvidia RTX 3080   10GP GDDR6
          Nvidia GTX 1050ti 4GB  GDDR5
    * Motherboard: ASUS TUF GAMING X570-PLUS
    * OS: Window/Ubuntu 20.04