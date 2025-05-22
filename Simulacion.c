#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "GPR.h"
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
#define bu -0.55
#define bv 0.0
#define cu 0.03
#define cv -0.015
#define Du 0.02
#define Dv 1.5
#define Fmax 0.2
#define Gmax 0.5

// Grid arrays
double *u, *v;
int *x_menos, *x_mas, *y_menos, *y_mas;

// Function to initialize the grid with random values
void initialize_grid() {
    srand(time(NULL));
    for (int i = 0; i < Nx * Ny; i++) {
        u[i] = (double)rand() / RAND_MAX;
        v[i] = (double)rand() / RAND_MAX;
    }
}

// Function to compute neighbor indices with periodic boundaries
void compute_neighbors() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int n = i + j * Nx;
            x_menos[n] = ((i - 1 + Nx) % Nx) + j * Nx;
            x_mas[n] = ((i + 1) % Nx) + j * Nx;
            y_menos[n] = i + ((j - 1 + Ny) % Ny) * Nx;
            y_mas[n] = i + ((j + 1) % Ny) * Nx;
        }
    }
}

// Evolution function
void evolucion() {
    double *aux_u = (double *)malloc(Nx * Ny * sizeof(double));
    double *aux_v = (double *)malloc(Nx * Ny * sizeof(double));

    for (int n = 0; n < Nx * Ny; n++) {
        double F = fmax(-Fmax, fmin((au - du) * u[n] + bu * v[n] + cu, Fmax));
        double G = fmax(-Gmax, fmin(av * u[n] + (bv - dv) * v[n] + cv, Gmax));

        aux_u[n] = u[n] + dt * (F + Du * (u[x_menos[n]] + u[x_mas[n]] + u[y_menos[n]] + u[y_mas[n]] - 4 * u[n]));
        aux_v[n] = v[n] + dt * (G + Dv * (v[x_menos[n]] + v[x_mas[n]] + v[y_menos[n]] + v[y_mas[n]] - 4 * v[n]));
    }

    for (int n = 0; n < Nx * Ny; n++) {
        u[n] = fmax(0.0, fmin(1.0, aux_u[n]));
        v[n] = fmax(0.0, fmin(1.0, aux_v[n]));
    }

    free(aux_u);
    free(aux_v);
}

// Function to save the final result as an image
void save_image() {
    unsigned char *image = (unsigned char *)malloc(Nx * Ny * sizeof(unsigned char));

    // Normalize values to 0-255 grayscale
    for (int i = 0; i < Nx * Ny; i++) {
        image[i] = (unsigned char)(u[i] * 255);
    }

    stbi_write_png("result.png", Nx, Ny, 1, image, Nx);
    printf("Image saved as result.png\n");

    free(image);
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

int main() {
    // Allocate memory
    u = (double *)malloc(Nx * Ny * sizeof(double));
    v = (double *)malloc(Nx * Ny * sizeof(double));
    x_menos = (int *)malloc(Nx * Ny * sizeof(int));
    x_mas = (int *)malloc(Nx * Ny * sizeof(int));
    y_menos = (int *)malloc(Nx * Ny * sizeof(int));
    y_mas = (int *)malloc(Nx * Ny * sizeof(int));

    if (!u || !v || !x_menos || !x_mas || !y_menos || !y_mas) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    initialize_grid();
    compute_neighbors();

    for (int t = 0; t < STEPS; t++) {
        evolucion();
        if ((t + 1) % 1000 == 0) {
            printf("Tiempo: %d/%d\n", (int)((t + 1)*dt), (int)(STEPS*dt));
        }

        if ((t + 1) % 100000 == 0)    save_image();
    }

    save_txt_output("squared_result.txt");
    save_image();

    // Free allocated memory
    free(u);
    free(v);
    free(x_menos);
    free(x_mas);
    free(y_menos);
    free(y_mas);

    return 0;
}
