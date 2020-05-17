#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <cmath>

using namespace std;

int n;
int s;
int cube_edge_size;
int world_rank;
float global_err;

struct Point {
    float x;
    float y;
    float z;
};

float exact(Point p) {
    return p.x * p.y * p.z;
}

float func(Point p) {
    return 0;
}

Point posInMainCubeDomain(int x, int y, int z) {
    Point p;
    p.x = (cube_edge_size * (int) (world_rank / s) + x) / n;
    p.y = (cube_edge_size * (int) (world_rank % s) + y) / n;
    p.z = (cube_edge_size * (int) (world_rank / (s * s)) + z) / n;
    return p;
}

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Datatype MPI_xy_plane;
    MPI_Datatype MPI_yz_plane;
    MPI_Datatype MPI_xz_plane;
    MPI_Type_vector(cube_edge_size * cube_edge_size,1, cube_edge_size, MPI_FLOAT, &MPI_xy_plane) ;
    MPI_Type_vector(cube_edge_size, cube_edge_size,cube_edge_size*cube_edge_size, MPI_FLOAT, &MPI_xz_plane) ;
    MPI_Type_vector(1, cube_edge_size*cube_edge_size, 0, MPI_FLOAT, &MPI_yz_plane) ;
    MPI_Request req[12];
    MPI_Status stats[12];

    n = atof(argv[1]);
    s = cbrt(world_size);
    cube_edge_size = (n / s) + 2; // +2 is for ghost points
    float arr[2][cube_edge_size][cube_edge_size][cube_edge_size];
    float error_threshold = 10E-4;
    int prev = 0;
    int curr = 1;
    global_err = 10E4;

    // initialize sub cube
    for (int i = 0; i < cube_edge_size; i++) {
        for (int j = 0; i < cube_edge_size; i++) {
            for (int k = 0; i < cube_edge_size; i++) {
                Point p = posInMainCubeDomain(i-1, j-1, k-1);
                if (p.x == 0 or p.x == 1 or p.y == 0 or p.y == 1 or p.z == 0 or p.z == 1) {
                    arr[prev][i][j][k] = exact(p);
                } else {
                    arr[prev][i][j][k] = 0;
                }
            }
        }
    }


//jacobi iteration
    while (global_err > error_threshold) {
        float local_err = 0;

        // send values to other processors - non blocking
        if(world_rank % s != 0){
            // not on the left surface of the cube
            MPI_Isend(&arr[curr][0][1][0], 1, MPI_xz_plane, world_rank - 1, 2, MPI_COMM_WORLD, &req[2]);
        }
        if(world_rank % s != s-1){
            // not on the right surface of the cube
            MPI_Isend(&arr[curr][0][cube_edge_size-2][0], 1, MPI_xz_plane, world_rank + 1, 5, MPI_COMM_WORLD, &req[5]);
        }
        if(world_rank % (s*s) >= s){
            // not on the up surface of the cube
            MPI_Isend(&arr[curr][1][0][0], 1, MPI_yz_plane, world_rank - s, 1, MPI_COMM_WORLD, &req[1]);
        }
        if(world_rank % (s*s) < (s-1) * s){
            // not on the bottom surface of the cube
            MPI_Isend(&arr[curr][cube_edge_size-2][0][0], 1, MPI_yz_plane, world_rank + s, 4, MPI_COMM_WORLD, &req[4]);
        }
        if(world_rank >= (s*s)){
            // not on the front surface of the cube
            MPI_Isend(&arr[curr][0][0][1], 1, MPI_xy_plane, world_rank-(s*s), 0, MPI_COMM_WORLD, &req[0]);
        }
        if(world_rank < world_size - (s*s)){
            // not on the back surface of the cube
            MPI_Isend(&arr[curr][0][0][cube_edge_size-2], 1, MPI_xy_plane, world_rank+(s*s), 3, MPI_COMM_WORLD, &req[3]);
        }


        // get values from other processors - blocking
        if(world_rank % s != 0){
            // not on the left surface of the cube
            MPI_Irecv(&arr[prev][0][cube_edge_size-1][0], 1, MPI_xz_plane, world_rank - 1, 5, MPI_COMM_WORLD, &req[11]);
        }
        if(world_rank % s != s-1){
            // not on the right surface of the cube
            MPI_Irecv(&arr[prev][0][0][0], 1, MPI_xz_plane, world_rank + 1, 2, MPI_COMM_WORLD, &req[8]);
        }
        if(world_rank % (s*s) >= s){
            // not on the up surface of the cube
            MPI_Irecv(&arr[prev][cube_edge_size-1][0][0], 1, MPI_yz_plane, world_rank - s, 4, MPI_COMM_WORLD, &req[10]);
        }
        if(world_rank % (s*s) < (s-1) * s){
            // not on the bottom surface of the cube
            MPI_Irecv(&arr[prev][0][0][0], 1, MPI_yz_plane, world_rank + s, 1, MPI_COMM_WORLD, &req[7]);
        }
        if(world_rank >= (s*s)){
            // not on the front surface of the cube
            MPI_Irecv(&arr[prev][0][0][cube_edge_size-1], 1, MPI_xy_plane, world_rank-(s*s), 3, MPI_COMM_WORLD, &req[9]);
        }
        if(world_rank < world_size - (s*s)){
            // not on the back surface of the cube
            MPI_Irecv(&arr[prev][0][0][0], 1, MPI_xy_plane, world_rank+(s*s), 0, MPI_COMM_WORLD, &req[6]);
        }

        MPI_Waitall(12,req,stats);

        for (int i = 1; i < cube_edge_size-1; i++) {
            for (int j = 1; i < cube_edge_size-1; i++) {
                for (int k = 1; i < cube_edge_size-1; i++) {
                    // calculate new values
                    Point p = posInMainCubeDomain(i-1, j-1, k-1);
                    arr[curr][i][j][k] = 1 / 6 *
                                         (arr[prev][i - 1][j][k] + arr[prev][i + 1][j][k] + arr[prev][i][j - 1][k] +
                                          arr[prev][i][j + 1][k] + arr[prev][i][j][k - 1] + arr[prev][i][j][k + 1]) -
                                         (1 / (6 * (n - 1) * (n - 1))) * func(p);
                    // calculate error
                    local_err += abs(arr[curr][i][j][k] - arr[prev][i][j][k]);
                }
            }
        }
        prev = 1 - prev;
        curr = 1 - curr;
        MPI_Allreduce(&local_err, &global_err, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

    if (world_rank == 0) {
        printf("%f", global_err);
    }


    MPI_Finalize();

    return 0;
}