#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define Nx 25
#define Ny 100
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

double *u, *v;
int neighbors[6][Nx * Ny];

void initialize_grid() {
    srand(42);
    for (int i = 0; i < Nx * Ny; i++) {
        u[i] = (double)rand() / RAND_MAX;
        v[i] = (double)rand() / RAND_MAX;
    }
}

void compute_neighbors() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            int idx = x + y * Nx;
            int dxs[6], dys[6];
            if (y % 2 == 1) { // odd row
                int dx_tmp[6] = {0, 0, 0, 1, 0, 1};
                int dy_tmp[6] = {-2, 2, -1, -1, 1, 1};
                for (int i = 0; i < 6; i++) {
                    dxs[i] = dx_tmp[i];
                    dys[i] = dy_tmp[i];
                }
            } else { // even row
                int dx_tmp[6] = {0, 0, -1, 0, -1, 0};
                int dy_tmp[6] = {-2, 2, -1, -1, 1, 1};
                for (int i = 0; i < 6; i++) {
                    dxs[i] = dx_tmp[i];
                    dys[i] = dy_tmp[i];
                }
            }

            for (int i = 0; i < 6; i++) {
                int nx = (x + dxs[i] + Nx) % Nx;
                int ny = (y + dys[i] + Ny) % Ny;
                neighbors[i][idx] = nx + ny * Nx;
            }
        }
    }
}

void evolucion() {
    double *aux_u = malloc(Nx * Ny * sizeof(double));
    double *aux_v = malloc(Nx * Ny * sizeof(double));

    for (int n = 0; n < Nx * Ny; n++) {
        double F = fmax(-Fmax, fmin((au - du) * u[n] + bu * v[n] + cu, Fmax));
        double G = fmax(-Gmax, fmin(av * u[n] + (bv - dv) * v[n] + cv, Gmax));

        double u_sum = 0.0, v_sum = 0.0;
        for (int i = 0; i < 6; i++) {
            u_sum += u[neighbors[i][n]];
            v_sum += v[neighbors[i][n]];
        }

        aux_u[n] = u[n] + dt * (F + Du * (u_sum - 6.0 * u[n]));
        aux_v[n] = v[n] + dt * (G + Dv * (v_sum - 6.0 * v[n]));
    }

    for (int n = 0; n < Nx * Ny; n++) {
        u[n] = fmax(0.0, fmin(1.0, aux_u[n]));
        v[n] = fmax(0.0, fmin(1.0, aux_v[n]));
    }

    free(aux_u);
    free(aux_v);
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

void save_image() {
    unsigned char *image = malloc(Nx * Ny);
    for (int i = 0; i < Nx * Ny; i++) {
        image[i] = (unsigned char)(u[i] * 255);
    }
    stbi_write_png("pattern_hexagonal.png", Nx, Ny, 1, image, Nx);
    free(image);
    printf("Saved final state as pattern_hexagonal.png\n");
}

int main() {
    u = malloc(Nx * Ny * sizeof(double));
    v = malloc(Nx * Ny * sizeof(double));

    if (!u || !v) {
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

    save_txt_output("hex_result.txt");
    save_image();

    free(u);
    free(v);
    return 0;
}
