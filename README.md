# skill_demonstration
This repository is a hodgepodge of my previous works in Computational Physics.

* Particle in Cell Method (1-D plasma simulation for mimicing Landau damping) 

* Mesh Generator (Delauney Triangle)
  * Bowyer-Watson algorithm
  wiki:https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm

* Conjugate Gradient Method (Linear Equation Solver for 1-D heat diffusion problem)
  * including CUDA version

* Finite Element Method (Material Mechanics)
  * Using above Mesh Generator to greate triangular mesh
  * Using serial Cojugate Gradiend code as solver

* Shallow Water Equation (1-D shallow water shock wave with Finite Volume Method)
  * including AVX, MPI, OpenMP, OpenCL,and CUDA version 
  * Some hybrid version:
    * MPI + CUDA for multi GPU usage
    * OpenMP + AVX for further speedup (thread level + instruction level parallelism)
    * OpenMP + MPI

* Euler Equation (1-D shock wave with finite volume method + AUSM scheme)
  * including AVX, MPI, OpenMP, OpenCL, and CUDA version
  
-------------------------
* Validation Platform:
  * 13 inch Macbook Pro (2020 intel CPU version)
  * Self-constructed PC
    * CPU: AMD Ryzen 3700x (8 physical cores, 16 cores with hyper threads)
    * DRAM: DDR4 16 GB 
    * GPU: Nvidia RTX 3080   10GP GDDR6
          Nvidia GTX 1050ti 4GB  GDDR5
    * Motherboard: ASUS TUF GAMING X570-PLUS
    * OS: Window/Ubuntu 20.04