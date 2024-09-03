
# Parallel Marching Squares using Pthreads

## Project Overview
This project is an implementation of the Marching Squares algorithm designed to generate contour lines for topographic maps. The main goal of the project is to optimize an existing sequential implementation by parallelizing it using Pthreads in C/C++. The parallel implementation is expected to produce the same output as the sequential version but with improved execution times.
## Marching Squares Algorithm
The Marching Squares algorithm is a computer graphics algorithm that generates contour lines (isolines) for a two-dimensional scalar field. It operates on a grid of cells, each containing four data points (corners), and determines the contour segment inside each cell based on the values at the corners.
## Sequential Implementation
The sequential implementation of the algorithm processes each grid cell one by one to determine the contour segments and then writes the result to an output file. This approach is straightforward but can be slow for large grids.
## Parallel Implementation
To improve the performance, the grid is divided into smaller blocks, and each block is processed independently in parallel using Pthreads. By distributing the workload across multiple threads, the execution time can be significantly reduced, especially for large data sets.

## Run the program:

```bash
- ./tema1_par < Input file > <Output file> <Number of threads>
```

- Input file containing the grid data.
- Output file where the generated contours will be saved.
- Number of threads to be used for parallel execution.
