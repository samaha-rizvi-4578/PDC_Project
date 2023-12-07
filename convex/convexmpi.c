#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

typedef struct {
    int x;
    int y;
} Point;

double calculate_distance(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

Point* generate_points(int n, int rank) {
    Point* points = (Point*)malloc(n * sizeof(Point));

    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            points[i].x = rand() % n + 1;
            points[i].y = rand() % n + 1;
        }
    }

    // Broadcast generated points to all processes
    MPI_Bcast(points, n * 2, MPI_INT, 0, MPI_COMM_WORLD);

    return points;
}

// Check if three points make a clockwise turn
int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; // colinear
    return (val > 0) ? 1 : 2; // clock or counterclockwise
}

// Check if point q is on the left side of the line formed by p and r
int on_left_side(Point p, Point q, Point r) {
    return orientation(p, q, r) == 2;
}

int count_convex_hull_points(Point* points, int n, int rank) {
    if (n < 3) {
        printf("Convex hull not possible with less than 3 points.\n");
        return 0;
    }

    int convexHullCount = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int i = 0; i < n - 2; ++i) {
        for (int j = i + 1; j < n - 1; ++j) {
            for (int k = j + 1; k < n; ++k) {
                int is_convex = 1;

                for (int m = 0; m < n; ++m) {
                    if (m != i && m != j && m != k) {
                        if (on_left_side(points[i], points[j], points[k]) !=
                            on_left_side(points[i], points[j], points[m])) {
                            is_convex = 0;
                            break;
                        }
                    }
                }

                if (is_convex) {
                    convexHullCount++;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double cpu_time_used = end_time - start_time;

    // Communication Time
    MPI_Barrier(MPI_COMM_WORLD);
    double comm_start_time = MPI_Wtime();

    // Scatter convexHullCount to all processes
    int* allConvexHullCounts = (int*)malloc(world_size * sizeof(int));
    MPI_Gather(&convexHullCount, 1, MPI_INT, allConvexHullCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double comm_end_time = MPI_Wtime();
    double comm_time_used = comm_end_time - comm_start_time;

    // Print results only from process 0
    if (rank == 0) {
        printf("Number of points in the convex hull: %d\n", convexHullCount);
        printf("Computation time: %.6f seconds\n", cpu_time_used);
        printf("Communication time: %.6f seconds\n", comm_time_used);

        // Calculate total time
        double total_time = comm_time_used + cpu_time_used;
        printf("Total time: %.6f seconds\n", total_time);
    }

    free(allConvexHullCounts);

    return convexHullCount;
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    MPI_Init(&argc, &argv);

    int world_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numPoints;

    if (rank == 0) {
        printf("Enter the number of points: ");
        scanf("%d", &numPoints);
    }

    MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Point* generated_points = generate_points(numPoints, rank);

    int convexHullCount = count_convex_hull_points(generated_points, numPoints, rank);

    free(generated_points);

    MPI_Finalize();

    return 0;
}
