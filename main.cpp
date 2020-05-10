#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <bits/stdc++.h>
#include <time.h>

using namespace std;

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int n = argv[];
    if(world_rank == 0){
        
    }
    int[][][][];


    MPI_Finalize();
    return 0;
}