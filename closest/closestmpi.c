#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <mpi.h>

typedef struct {
    int x;
    int y;
} Point;

double calculate_distance(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

Point* generate_points(int n) {
    Point* points = (Point*)malloc(n * sizeof(Point));

    for (int i = 0; i < n; ++i) {
        again:
        {
            points[i].x = rand() % (n + 1);
            points[i].y = rand() % (n + 1);
        }
        if (points[i].x == points[i].y) {
            goto again;
        }
    }

    return points;
}

void write_points_to_file(Point* points, int n, int rank) {
    char filename[20];
    sprintf(filename, "points_%d.txt", rank);
    FILE* file = fopen(filename, "w");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; ++i) {
        fprintf(file, "%d %d\n", points[i].x, points[i].y);
    }

    fclose(file);
}

void find_shortest_distance(Point* points, int n, int* closest_pair_index1, int* closest_pair_index2) {
    double min_distance = DBL_MAX;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double distance = calculate_distance(points[i], points[j]);
            if (distance < min_distance) {
                min_distance = distance;
                *closest_pair_index1 = i;
                *closest_pair_index2 = j;
            }
        }
    }
}

int main(int argc, char** argv) {
    srand(time(NULL));

    MPI_Init(&argc, &argv);

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n;
    if (rank == 0) {
        printf("Enter the number of points (n): ");
        scanf("%d", &n);
    }

    // Broadcast n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Point* generated_points = generate_points(n);

    // Scatter the points to all processes
    Point* local_points = (Point*)malloc((n / world_size) * sizeof(Point));
    MPI_Scatter(generated_points, n / world_size, MPI_2INT, local_points, n / world_size, MPI_2INT, 0, MPI_COMM_WORLD);

    write_points_to_file(local_points, n / world_size, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    int closest_pair_index1, closest_pair_index2;
    find_shortest_distance(local_points, n / world_size, &closest_pair_index1, &closest_pair_index2);

    // Gather the closest pair indices from all processes to rank 0
    int* all_indices = (int*)malloc(2 * world_size * sizeof(int));
    MPI_Gather(&closest_pair_index1, 1, MPI_INT, all_indices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&closest_pair_index2, 1, MPI_INT, all_indices + world_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double shortest_distance;
    int final_closest_pair_index1, final_closest_pair_index2;

    if (rank == 0) {
        int overall_closest_pair_index1, overall_closest_pair_index2;
        shortest_distance = DBL_MAX;

        for (int i = 0; i < world_size; ++i) {
            double distance = calculate_distance(generated_points[all_indices[i]],
                                                 generated_points[all_indices[i + world_size]]);
            if (distance < shortest_distance) {
                shortest_distance = distance;
                overall_closest_pair_index1 = all_indices[i];
                overall_closest_pair_index2 = all_indices[i + world_size];
            }
        }

        final_closest_pair_index1 = overall_closest_pair_index1;
        final_closest_pair_index2 = overall_closest_pair_index2;
    }

    // Broadcast the final closest pair indices to all processes
    MPI_Bcast(&final_closest_pair_index1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_closest_pair_index2, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double cpu_time_used = end_time - start_time;

    if (rank == 0) {
        printf("Shortest distance between points: %.2f\n", shortest_distance);
        printf("Closest pair of points: (%d, %d) and (%d, %d)\n",
               generated_points[final_closest_pair_index1].x, generated_points[final_closest_pair_index1].y,
               generated_points[final_closest_pair_index2].x, generated_points[final_closest_pair_index2].y);
        printf("Time taken to calculate the shortest distance: %.6f seconds\n", cpu_time_used);
    }

    free(generated_points);
    free(local_points);
    free(all_indices);

    MPI_Finalize();

    return 0;
}
