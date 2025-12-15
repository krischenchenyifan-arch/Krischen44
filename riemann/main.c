#include <stdio.h>
#include <stdlib.h>
#include "riemann.h"

// Thermodynamic constants
#define GAMMA 1.4
#define R 1.0
#define CFL 0.25
#define FINAL_TIME 0.2
#define CV (R/(GAMMA - 1.0))
#define L 1.0
#define NX 200
#define DX (L/NX)

int main() {

    // Set up a simple shock tube problem
    const int NO_CELLS = 200;
    const int NO_INTERFACES = NO_CELLS + 1;
    // 1D Euler equations have 3 primitives: density, velocity, temperature
    double p0[NO_CELLS]; // Density
    double p1[NO_CELLS]; // Velocity
    double p2[NO_CELLS]; // Temperature

    // 1D Euler equations have 3 conservatives: density, momentum, energy
    double u0[NO_CELLS]; // Density
    double u1[NO_CELLS]; // Momentum
    double u2[NO_CELLS]; // Energy

    // Set up the fluxes of conserved quantities at each interface
    double f0[NO_INTERFACES]; // Density flux
    double f1[NO_INTERFACES]; // Momentum flux
    double f2[NO_INTERFACES]; // Energy flux

    /*
    Now we are ready to solve:

    du/dt + df/dx = 0

    */

    // Initialize the problem
    for (int i = 0; i < NO_CELLS; i++) {
        if (i < NO_CELLS / 2) {
            p0[i] = 10.0;   // Density
            p1[i] = 0.0;   // Velocity
            p2[i] = 1.0;   // Temperature
        } else {
            p0[i] = 1.0;
            p1[i] = 0.0;
            p2[i] = 1.0;
        }

        // Compute the conserved variables from the primitive variables
        u0[i] = p0[i];
        u1[i] = p0[i] * p1[i];
        u2[i] = p0[i] * (CV * p2[i] + 0.5 * p1[i] * p1[i]);
    }

    return 0;
}