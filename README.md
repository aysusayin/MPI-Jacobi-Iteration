# MPI-Jacobi-Iteration
## This is a homework for CMPE 478: Parallel Processing lecture in Bogazici University
### Solving Poisson Equation with Jacobi iteration
In this homework, we needed to solve systems of equations that arise from the
discretization of Poisson’s equation on a cube domain. The poissons equation is the following:
```
u_xx + u_yy + u_zz = f(x,y,z)
```
The domain is [0,1] and the value of the boundary points are given. The program takes argument N as input and divide the domain into a (n+1) x (n+1) + (n+1) cube grid. Using the finite difference method, the approximation of each grid point is found as:

```
u(x_i, y_j,z_k) = 1/6 * [u(x_i-1, y_j,z_k) + u(x_i+1, y_j,z_k) + u(x_i, y_j-1,z_k) + u(x_i, y_j+1,z_k) + u(x_i, y_j,z_k-1) + u(x_i, y_j,z_k+1)] - 1/(6n^2) * f(x_i, y_j,z_k)
```

In order to solve the equation, the boundary values must be known. In this work the value of the boundary points are calculated with exact() function and it is given as f = x*y*z.

To exploit MPI, I divided the domain into subcubes. Each processor works on a subcube of the domain and exchanges the values in the boundaries. I assigned the subcubes to processors as in the example. (The blue numbers are the processor ids.) Also added ghost cells to do the computations easily.

Example:

![alt text](https://i.ibb.co/dPBWNrZ/cube.jpg)

Then jacobi iteration is done in each processor. The local error is calculated by taking the difference of previous and the current values. Then all local errors are summed up using reduction method. When the global error is less then the predetermined error threshold, iteration stops. 


#### How to Run?
To compile the program: 

```
mpicxx main.cpp -o solver

```

Then run the program like this:

```
 mpiexec -n [Number of Processors] ./solver [N]
```

Number of processors is a perfect cube number and N+1 is divisible by the cubroot of number of processors. Program outputs the final error. If you don't want to display anything on the screen, put --silent flag before N.

#### Performance Analysis
To see the analysis, run `./analyze.sh`. You can run this file without making the project since it makes it.
 
 NOTE: "analyze.sh" program generates a file named "MPI_Analysis.csv" and "tmp_file". So
 please be careful if you have a file with the same name since the program overwrites it.
 
There is an example csv file that contains the results generated on my computer called "ExampleAnalysis.csv".

My CPU information:
```
Architecture:        x86_64
CPU op-mode(s):      32-bit, 64-bit
Byte Order:          Little Endian
CPU(s):              4
On-line CPU(s) list: 0-3
Thread(s) per core:  2
Core(s) per socket:  2
Socket(s):           1
NUMA node(s):        1
Vendor ID:           GenuineIntel
CPU family:          6
Model:               78
Model name:          Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz
Stepping:            3
CPU MHz:             3059.423
CPU max MHz:         3100.0000
CPU min MHz:         400.0000
BogoMIPS:            5184.00
Virtualization:      VT-x
L1d cache:           32K
L1i cache:           32K
L2 cache:            256K
L3 cache:            4096K
NUMA node0 CPU(s):   0-3
```
