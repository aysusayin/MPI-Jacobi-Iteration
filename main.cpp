#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <fstream>

using namespace std;

int N;
int s;
int cube_edge_size;
int world_rank;
const double ERROR_THRESHOLD = 10E-10;

struct Point {
  double x;
  double y;
  double z;
};

double exact(Point p) { return p.x * p.y * p.z; }

double func(Point p) { return 0; }

Point posInMainCubeDomain(int x, int y, int z) {
  Point p;
  p.x =
      (1.0 * (cube_edge_size - 2) * ((int)(world_rank % (s * s)) / s) + x) / N;
  p.y = (1.0 * (cube_edge_size - 2) * ((int)(world_rank % s)) + y) / N;
  p.z = (1.0 * (cube_edge_size - 2) * ((int)(world_rank / (s * s))) + z) / N;
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
  N = atof(argv[1]);
  s = cbrt(world_size);
  int prev = 0;
  int curr = 1;
  double global_err = 10E4;
  cube_edge_size = (N + 1) / s + 2;  // +2 is for ghost cells
  double arr[2][cube_edge_size][cube_edge_size][cube_edge_size];

  // initialize sub cube
  for (int i = 0; i < cube_edge_size; i++) {
    for (int j = 0; j < cube_edge_size; j++) {
      for (int k = 0; k < cube_edge_size; k++) {
        arr[prev][i][j][k] = 0;
        arr[curr][i][j][k] = 0;
        if (!(i == 0 or i == cube_edge_size - 1 or j == 0 or
              j == cube_edge_size - 1 or k == 0 or k == cube_edge_size - 1)) {
          Point p = posInMainCubeDomain(i - 1, j - 1, k - 1);
          if (p.x == 0 or p.x == 1 or p.y == 0 or p.y == 1 or p.z == 0 or
              p.z == 1) {
            // Boundary points
            double e = exact(p);
            arr[prev][i][j][k] = e;
            arr[curr][i][j][k] = e;
          }
        }
      }
    }
  }

  // Define datatypes
  MPI_Type_vector(cube_edge_size * cube_edge_size, 1, cube_edge_size,
                  MPI_DOUBLE, &MPI_xy_plane);
  MPI_Type_vector(cube_edge_size, cube_edge_size,
                  cube_edge_size * cube_edge_size, MPI_DOUBLE, &MPI_xz_plane);
  MPI_Type_vector(1, cube_edge_size * cube_edge_size, 0, MPI_DOUBLE,
                  &MPI_yz_plane);

  MPI_Type_commit(&MPI_xy_plane);
  MPI_Type_commit(&MPI_yz_plane);
  MPI_Type_commit(&MPI_xz_plane);

  MPI_Request req[12];
  for (int i = 0; i < 12; i++) {
    req[i] = MPI_REQUEST_NULL;
  }
  MPI_Status stats[12];

  // jacobi iteration
  while (global_err > ERROR_THRESHOLD) {
    double local_err = 0;
    prev = 1 - prev;
    curr = 1 - prev;
    // send values to other processors - non blocking
    if (world_rank % s != 0) {
      // not on the left surface of the cube
      MPI_Isend(&arr[curr][0][1][0], 1, MPI_xz_plane, world_rank - 1, 2,
                MPI_COMM_WORLD, &req[2]);
    }
    if (world_rank % s != s - 1) {
      // not on the right surface of the cube
      MPI_Isend(&arr[curr][0][cube_edge_size - 2][0], 1, MPI_xz_plane,
                world_rank + 1, 5, MPI_COMM_WORLD, &req[5]);
    }
    if (world_rank % (s * s) >= s) {
      // not on the up surface of the cube
      MPI_Isend(&arr[curr][1][0][0], 1, MPI_yz_plane, world_rank - s, 1,
                MPI_COMM_WORLD, &req[1]);
    }
    if (world_rank % (s * s) < (s - 1) * s) {
      // not on the bottom surface of the cube
      MPI_Isend(&arr[curr][cube_edge_size - 2][0][0], 1, MPI_yz_plane,
                world_rank + s, 4, MPI_COMM_WORLD, &req[4]);
    }
    if (world_rank >= (s * s)) {
      // not on the front surface of the cube
      MPI_Isend(&arr[curr][0][0][1], 1, MPI_xy_plane, world_rank - (s * s), 0,
                MPI_COMM_WORLD, &req[0]);
    }
    if (world_rank < world_size - (s * s)) {
      // not on the back surface of the cube
      MPI_Isend(&arr[curr][0][0][cube_edge_size - 2], 1, MPI_xy_plane,
                world_rank + (s * s), 3, MPI_COMM_WORLD, &req[3]);
    }

    // get values from other processors - blocking
    if (world_rank % s != 0) {
      // not on the left surface of the cube
      MPI_Irecv(&arr[prev][0][cube_edge_size - 1][0], 1, MPI_xz_plane,
                world_rank - 1, 5, MPI_COMM_WORLD, &req[11]);
    }
    if (world_rank % s != s - 1) {
      // not on the right surface of the cube
      MPI_Irecv(&arr[prev][0][0][0], 1, MPI_xz_plane, world_rank + 1, 2,
                MPI_COMM_WORLD, &req[8]);
    }
    if (world_rank % (s * s) >= s) {
      // not on the up surface of the cube
      MPI_Irecv(&arr[prev][cube_edge_size - 1][0][0], 1, MPI_yz_plane,
                world_rank - s, 4, MPI_COMM_WORLD, &req[10]);
    }
    if (world_rank % (s * s) < (s - 1) * s) {
      // not on the bottom surface of the cube
      MPI_Irecv(&arr[prev][0][0][0], 1, MPI_yz_plane, world_rank + s, 1,
                MPI_COMM_WORLD, &req[7]);
    }
    if (world_rank >= (s * s)) {
      // not on the front surface of the cube
      MPI_Irecv(&arr[prev][0][0][cube_edge_size - 1], 1, MPI_xy_plane,
                world_rank - (s * s), 3, MPI_COMM_WORLD, &req[9]);
    }
    if (world_rank < world_size - (s * s)) {
      // not on the back surface of the cube
      MPI_Irecv(&arr[prev][0][0][0], 1, MPI_xy_plane, world_rank + (s * s), 0,
                MPI_COMM_WORLD, &req[6]);
    }

    MPI_Waitall(12, req, stats);

    for (int i = 1; i < cube_edge_size - 1; i++) {
      for (int j = 1; j < cube_edge_size - 1; j++) {
        for (int k = 1; k < cube_edge_size - 1; k++) {
          // calculate new values
          Point p = posInMainCubeDomain(i - 1, j - 1, k - 1);
          if (!(p.x == 0 or p.x == 1 or p.y == 0 or p.y == 1 or p.z == 0 or
                p.z == 1)) {
            // Don't update boundry points
            arr[curr][i][j][k] =
                1.0 / 6 *
                    (arr[prev][i - 1][j][k] + arr[prev][i + 1][j][k] +
                     arr[prev][i][j - 1][k] + arr[prev][i][j + 1][k] +
                     arr[prev][i][j][k - 1] + arr[prev][i][j][k + 1]) -
                (1 / (6 * N * N)) * func(p);
            // calculate error
            local_err += abs(arr[curr][i][j][k] - arr[prev][i][j][k]);
          }
        }
      }
    }
    MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }

  if (world_rank == 0) {
    printf("global error:");
    cout << global_err << "\n";
  }

  MPI_Finalize();

  return 0;
}
