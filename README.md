# skill demonstration
This repository is a hodgepodge of my previous works in Computational Physics.
* Particle in Cell Method (1-D plasma simulation for mimicing Landau damping) 
* Mesh Generator (Delauney Triangle)
  * Bowyer-Watson algorithm
  wiki:https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm

* Conjugate Gradient Method (Linear equation solver for 1-D heat diffusion problem)
  * Including CUDA version

* Finite Element Method (2-D elastic solid mechanics)
  * Using Mesh Generator to create triangular mesh
  * Using serial Conjugate Gradient code as solver

* Shallow Water Equation (1-D shallow water shock wave with finite volume method)
  * Including AVX, MPI, OpenMP, OpenCL,and CUDA version 
  * Hybrid version:
    * MPI + CUDA for multi GPU usage
    * OpenMP + AVX for further speedup (thread level + instruction level parallelism)
    * OpenMP + MPI

* Euler Equation (1-D shock wave with finite volume method + AUSM scheme)
  * Including AVX, MPI, OpenMP, OpenCL, and CUDA version
  
-------------------------
Validation Platform:
* 13 inch Macbook Pro (2020 intel CPU version)
* Self-constructed PC
  * CPU: AMD Ryzen 3700X (8 physical cores, 16 cores with hyper threads)
  * DRAM: DDR4 16GB 
  * GPU: 
    * Nvidia RTX 3080, 10GB GDDR6
    * Nvidia GTX 1050ti, 4GB  GDDR5
  * Motherboard: ASUS TUF GAMING X570-PLUS
  * OS: Window 10 / Ubuntu 20.04