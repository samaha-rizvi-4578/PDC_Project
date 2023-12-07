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

void convexHull(Point points[], int n, FILE *outputFile) {
    if (n < 3) return;

    int hullSize = 0;
    Point hull[n];

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
    // Write convex hull points to the output file
    fprintf(outputFile, "Convex Hull:\n");
    for (int i = 0; i < hullSize; i++) {
        fprintf(outputFile, "(%d, %d)\n", hull[i].x, hull[i].y);
        printf("(%d, %d)\n", hull[i].x, hull[i].y);
    }
}

int main() {
    FILE *inputFile, *outputFile;
    inputFile = fopen("points.txt", "w");
    if (inputFile == NULL) {
        perror("Error opening file");
        return 1;
    }

    int numPoints;

    printf("Enter the number of points: ");
    scanf("%d", &numPoints);

    // Generate random points and write to the "points.txt" file
    fprintf(inputFile, "%d\n", numPoints);
    for (int i = 0; i < numPoints; i++) {
        Point p;
        p.x = rand() % 100;
        p.y = rand() % 100;
        fprintf(inputFile, "%d %d\n", p.x, p.y);
    }

    fclose(inputFile);

    // Open "points.txt" for reading
    inputFile = fopen("points.txt", "r");
    if (inputFile == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read points from "points.txt"
    fscanf(inputFile, "%d", &numPoints);
    Point points[numPoints];
    for (int i = 0; i < numPoints; i++)
        fscanf(inputFile, "%d %d", &points[i].x, &points[i].y);

    fclose(inputFile);



    // Open "convex_hull.txt" for writing
    outputFile = fopen("convex_hull.txt", "w");
    if (outputFile == NULL) {
        perror("Error creating file");
        return 1;
    }

    // Measure time for convex hull calculation
    clock_t start = clock();
    convexHull(points, numPoints, outputFile);
    clock_t end = clock();

    // Calculate the time taken
    double timeTaken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken to calculate convex hull: %f seconds\n", timeTaken);

    fclose(outputFile);

    return 0;
}
