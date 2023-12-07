#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int x, y;
} Point;

int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; // colinear
    return (val > 0) ? 1 : 2; // clock or counterclockwise
}

void convexHull(Point points[], int n) {
    if (n < 3) return;

    int hullSize = 0;
    Point *hull = malloc(n * sizeof(Point));
    if (hull == NULL) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
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

    printf("Convex Hull Points:\n");
    for (int i = 0; i < hullSize; i++) {
        printf("(%d, %d)\n", hull[i].x, hull[i].y);
    }

    free(hull);  // Don't forget to free the allocated memory
}

int main() {
    int numPoints;

    printf("Enter the number of points: ");
    scanf("%d", &numPoints);

    // Generate random points
    Point *points = malloc(numPoints * sizeof(Point));
    if (points == NULL) {
        perror("Memory allocation error");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < numPoints; i++) {
        points[i].x = rand() % 100;
        points[i].y = rand() % 100;
    }



    // Measure time for convex hull calculation
    clock_t start = clock();
    convexHull(points, numPoints);
    clock_t end = clock();

    // Calculate the time taken
    double timeTaken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken to calculate convex hull: %f seconds\n", timeTaken);

    free(points);  // Don't forget to free the allocated memory

    return 0;
}
