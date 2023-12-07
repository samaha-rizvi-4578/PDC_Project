#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

Point* generate_points(int n) {
    Point* points = (Point*)malloc(n * sizeof(Point));

    for (int i = 0; i < n; ++i) {
        points[i].x = rand() % n + 1;
        points[i].y = rand() % n + 1;
    }

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

int count_convex_hull_points(Point* points, int n) {
    if (n < 3) {
        printf("Convex hull not possible with less than 3 points.\n");
        return 0;
    }

    int convexHullCount = 0;

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

    return convexHullCount;
}

int main() {
    srand(time(NULL));

    int n;
    printf("Enter the number of points (n): ");
    scanf("%d", &n);

    Point* generated_points = generate_points(n);

    clock_t start_time = clock();

    int convexHullCount = count_convex_hull_points(generated_points, n);

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Number of points in the convex hull: %d\n", convexHullCount);
    printf("Time taken to calculate convex hull: %.6f seconds\n", cpu_time_used);

    free(generated_points);

    return 0;
}
