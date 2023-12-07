#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

typedef struct {
    int x, y;
} Point;

int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; // colinear
    return (val > 0) ? 1 : 2; // clock or counterclockwise
}

void convexHull(Point points[], int n, int rank) {
    if (n < 3) return;

    int hullSize = 0;
    Point *hull = malloc(n * sizeof(Point));
    if (hull == NULL) {
        perror("Memory allocation error");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int leftmost = 0;
    for (int i = 1; i < n; i++)
        if (points[i].x < points[leftmost].x)
            leftmost = i;

    int p = leftmost, q;
    do {
        hull[hullSize++] = points[p];

        q = (p + 1) % n;
        for (int i = 0; i < n; i++) {
            if (orientation(points[p], points[i], points[q]) == 2)
                q = i;
        }

        p = q;
    } while (p != leftmost);

    printf("Process %d: Convex Hull Points:\n", rank);
    for (int i = 0; i < hullSize; i++) {
        printf("(%d, %d)\n", hull[i].x, hull[i].y);
    }

    free(hull);
}

int main(int argc, char *argv[]) {
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

    // Generate random points
    Point *points = malloc(numPoints * sizeof(Point));
    if (points == NULL) {
        perror("Memory allocation error");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (rank == 0) {
        for (int i = 0; i < numPoints; i++) {
            points[i].x = rand() % 100;
            points[i].y = rand() % 100;
        }
    }

    // Broadcast the generated points to all processes
    MPI_Bcast(points, numPoints * 2, MPI_INT, 0, MPI_COMM_WORLD);


    // Measure time for convex hull calculation
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    convexHull(points, numPoints, rank);
    double end = MPI_Wtime();

   // if (rank == 0) {
        // Calculate the time taken
     //   printf("Time taken to calculate convex hull: %f seconds\n", end - start);
   // }
	  printf("Time taken to calculate convex hull: %f seconds\n", end - start);
    free(points); // Don't forget to free the allocated memory

    MPI_Finalize();

    return 0;
}
