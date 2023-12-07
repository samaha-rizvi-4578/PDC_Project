//closest  opemp:
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <omp.h>

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
        if(points[i].x == points[i].y)
        {
            goto again;
        }
    }

    return points;
}

void write_points_to_file(Point* points, int n) {
    FILE* file = fopen("points.txt", "w");

    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < n; ++i) {
        #pragma omp critical
        {
            fprintf(file, "%d %d\n", points[i].x, points[i].y);
        }
    }

    fclose(file);
}

void find_shortest_distance(Point* points, int n, int* closest_pair_index1, int* closest_pair_index2) {
    double min_distance = DBL_MAX;

    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double distance = calculate_distance(points[i], points[j]);
            #pragma omp critical
            {
                if (distance < min_distance) {
                    min_distance = distance;
                    *closest_pair_index1 = i;
                    *closest_pair_index2 = j;
                }
            }
        }
    }
}

int main() {
    srand(time(NULL));

    int n;
    printf("Enter the number of points (n): ");
    scanf("%d", &n);

    Point* generated_points = generate_points(n);

    write_points_to_file(generated_points, n);

    clock_t start_time = clock();

    int closest_pair_index1, closest_pair_index2;
    find_shortest_distance(generated_points, n, &closest_pair_index1, &closest_pair_index2);
    double shortest_distance = calculate_distance(generated_points[closest_pair_index1], generated_points[closest_pair_index2]);

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Shortest distance between points: %.2f\n", shortest_distance);
    printf("Closest pair of points: (%d, %d) and (%d, %d)\n",
           generated_points[closest_pair_index1].x, generated_points[closest_pair_index1].y,
           generated_points[closest_pair_index2].x, generated_points[closest_pair_index2].y);
    printf("Time taken to calculate the shortest distance: %.6f seconds\n", cpu_time_used);

    free(generated_points);

    return 0;
}