#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define Nx 50
#define Ny 50
#define STEPS 500000

#define dt 0.1
#define du 0.03
#define dv 0.08
#define au 0.08
#define av 0.045
#define bu -0.6
#define bv 0.0
#define cu 0.03
#define cv -0.015
#define Du 0.035
#define Dv 1.5
#define Fmax 0.2
#define Gmax 0.5

// Grid arrays
double *u, *v;
int *x_menos, *x_mas, *y_menos, *y_mas;

void initialize_grid() {
    srand(42);  // Fixed seed for reproducibility
    for (int i = 0; i < Nx * Ny; i++) {
        u[i] = (double)rand() / RAND_MAX;
        v[i] = (double)rand() / RAND_MAX;
    }
}

void save_txt_output(const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Could not open output file");
        return;
    }
    for (int i = 0; i < Nx * Ny; i++) {
        fprintf(f, "%.4f\n", u[i]);
    }
    fclose(f);
    printf("Saved final state to %s\n", filename);
}

void compute_neighbors() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int n = i + j * Nx;
            x_menos[n] = ((i - 1 + Nx) % Nx) + j * Nx;
            x_mas[n] = ((i + 1) % Nx) + j * Nx;
            y_menos[n] = (i + 1) % Nx + ((j - 1 + Ny) % Ny) * Nx;
            y_mas[n] = (i - 1 + Nx) % Nx + ((j + 1) % Ny) * Nx;
        }
    }
}

void evolucion() {
    double *aux_u = malloc(Nx * Ny * sizeof(double));
    double *aux_v = malloc(Nx * Ny * sizeof(double));

    for (int n = 0; n < Nx * Ny; n++) {
        double F = fmax(-Fmax, fmin((au - du) * u[n] + bu * v[n] + cu, Fmax));
        double G = fmax(-Gmax, fmin(av * u[n] + (bv - dv) * v[n] + cv, Gmax));

        int y_vecino = (n % 2 == 0) ? y_menos[n] : y_mas[n];

        aux_u[n] = u[n] + dt * (F + Du * (u[x_menos[n]] + u[x_mas[n]] + u[y_vecino] - 3 * u[n]));
        aux_v[n] = v[n] + dt * (G + Dv * (v[x_menos[n]] + v[x_mas[n]] + v[y_vecino] - 3 * v[n]));
    }

    for (int n = 0; n < Nx * Ny; n++) {
        u[n] = fmax(0.0, fmin(1.0, aux_u[n]));
        v[n] = fmax(0.0, fmin(1.0, aux_v[n]));
    }

    free(aux_u);
    free(aux_v);
}

void save_image() {
    unsigned char *image = malloc(Nx * Ny * sizeof(unsigned char));
    for (int i = 0; i < Nx * Ny; i++) {
        image[i] = (unsigned char)(u[i] * 255);
    }
    stbi_write_png("pattern_triangular.png", Nx, Ny, 1, image, Nx);
    printf("Saved final state as pattern_triangular.png\n");
    free(image);
}

int main() {
    u = malloc(Nx * Ny * sizeof(double));
    v = malloc(Nx * Ny * sizeof(double));
    x_menos = malloc(Nx * Ny * sizeof(int));
    x_mas = malloc(Nx * Ny * sizeof(int));
    y_menos = malloc(Nx * Ny * sizeof(int));
    y_mas = malloc(Nx * Ny * sizeof(int));

    if (!u || !v || !x_menos || !x_mas || !y_menos || !y_mas) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    initialize_grid();
    compute_neighbors();

    for (int t = 0; t < STEPS; t++) {
        evolucion();
        if ((t + 1) % 1000 == 0)
            printf("Step: %d/%d\n", t + 1, STEPS);
    }

    save_txt_output("tri_result.txt");

    free(u); free(v);
    free(x_menos); free(x_mas); free(y_menos); free(y_mas);

    return 0;
}
