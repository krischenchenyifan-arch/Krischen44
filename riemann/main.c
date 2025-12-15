#include "riemann.h"

#define DEBUG 0

// Thermodynamic constants
#define GAMMA 1.4
#define R 1.0
#define CFL 0.25
#define FINAL_TIME 0.2
#define CV (R/(GAMMA - 1.0))
#define L 1.0
#define NX 200
#define DX (L/NX)
#define MAX_ITERATIONS 1000

// We are doing a 1D shock tube problem
// This means all of the cell interface normals are exactly the same
double nx = 1.0; double ny = 0.0; double nz = 0.0;
double px = 0.0; double py = 1.0; double pz = 0.0;
double qx = 0.0; double qy = 0.0; double qz = 1.0;


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

    float DT = 1e-6; // Initial time step; this will be updated dynamically

    /*
    Now we are ready to solve:

    du/dt + df/dx = 0

    */

    // Initialize the problem
    for (int i = 0; i < NO_CELLS; i++) {
        if (i < 0.5*NO_CELLS) {
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

    // Start some time stepping
    float time = 0.0;
    for (int timestep = 0; timestep < MAX_ITERATIONS; timestep++) {
        // Compute the maximum CFL condition
        double MAX_CFL = CPU_Compute_Max_CFL(p0, p1, p2, DX, DT, NO_CELLS);
        // Update the time step based on the CFL condition
        DT = DT * ( CFL / MAX_CFL );   // IF CFL = 0.25 but MAX_CFL = 0.5, then DT_new = DT_old * (0.25/0.5) = 0.5 * DT_old 

        // Compute fluxes at each interface
        printf("Computing fluxes at timestep %d with DT = %e\n", timestep, DT);
        for (int i = 0; i < NO_INTERFACES; i++) {
            double QL_rho, QL_vx, QL_vy, QL_vz, QL_cRT;
            double QR_rho, QR_vx, QR_vy, QR_vz, QR_cRT;

            if (i == 0) {
                // Use a reflective boundary condition on the left
                QL_rho = p0[0];
                QL_vx = -p1[0];
                QL_vy = 0.0;
                QL_vz = 0.0;
                QL_cRT = sqrt(R * p2[0]);
                QR_rho = p0[0];
                QR_vx = p1[0];
                QR_vy = 0.0;
                QR_vz = 0.0;
                QR_cRT = sqrt(R * p2[0]);
            } else if (i == NO_INTERFACES - 1) {
                // Use a reflective boundary condition on the right
                QL_rho = p0[NO_CELLS - 1];
                QL_vx = p1[NO_CELLS - 1];
                QL_vy = 0.0;
                QL_vz = 0.0;
                QL_cRT = sqrt(R * p2[NO_CELLS - 1]);
                QR_rho = p0[NO_CELLS - 1];
                QR_vx = -p1[NO_CELLS - 1];
                QR_vy = 0.0;
                QR_vz = 0.0;
                QR_cRT = sqrt(R * p2[NO_CELLS - 1]);
            } else {
                // Normal case: take values from adjacent cells
                QL_rho = p0[i-1];
                QL_vx = p1[i-1];
                QL_vy = 0.0;
                QL_vz = 0.0;
                QL_cRT = sqrt(R * p2[i-1]);
                QR_rho = p0[i];
                QR_vx = p1[i];
                QR_vy = 0.0;
                QR_vz = 0.0;
                QR_cRT = sqrt(R * p2[i]);
            }

            // Compute the flux at this interface
            double flux[5]; // The solver is a 3D solver, though we'll only use 1D here
            double interface_p[5];
            if (DEBUG) {
                printf("Computing flux at interface %d between cells %d and %d\n", i, i - 1, i);
                printf("  Left state: rho = %g, vx = %g, T = %g\n", QL_rho, QL_vx, (QL_cRT * QL_cRT) / R);
                printf("  Right state: rho = %g, vx = %g, T = %g\n", QR_rho, QR_vx, (QR_cRT * QR_cRT) / R);
                printf("-------------------------------------------\n");
            }
            CPU_Calc_Flux(flux, interface_p,
                QL_rho, QL_vx, QL_vy, QL_vz, QL_cRT,
                QR_rho, QR_vx, QR_vy, QR_vz, QR_cRT, R, GAMMA,
                nx, ny, nz,
                px, py, pz,
                qx, qy, qz, 0); // wall_flag = 0; we handle the wall situation above
            
            // Store the fluxes of conserved variables
            f0[i] = flux[0]; // Density flux
            f1[i] = flux[1]; // Momentum flux
            f2[i] = flux[4]; // Energy flux (note: flux[2] is y-momentum, flux[3] is z-momentum)

        } // End loop over interfaces

        // Update the conserved variables in each cell using the fluxes
        for (int i = 0; i < NO_CELLS; i++) {
            if (DEBUG) printf("Upating cell %d: left mom flux = %g, right mom flux = %g\n", i, f1[i], f1[i + 1]);
            u0[i] = u0[i] - (DT / DX) * (f0[i + 1] - f0[i]);
            u1[i] = u1[i] - (DT / DX) * (f1[i + 1] - f1[i]);
            u2[i] = u2[i] - (DT / DX) * (f2[i + 1] - f2[i]);
        }           

        // Update the primitive variables from the conserved variables
        for (int i = 0; i < NO_CELLS; i++) {
            p0[i] = u0[i]; // Density
            p1[i] = u1[i] / u0[i]; // Velocity
            double kinetic_energy = 0.5 * p1[i] * p1[i];
            double total_energy_per_mass = u2[i] / u0[i];
            double internal_energy_per_mass = total_energy_per_mass - kinetic_energy;
            p2[i] = internal_energy_per_mass / CV; // Temperature   
        }      

        time += DT;
        if (time >= FINAL_TIME) {
            printf("Reached final time at timestep %d: time = %g\n", timestep, time);
            break;
        }

    }  // End time stepping

    // Save the results
    FILE *output_file = fopen("shock_tube_results.dat", "w");
    for (int i = 0; i < NO_CELLS; i++) {
        float x = (i + 0.5) * DX; // Cell center
        fprintf(output_file, "%g\t%g\t%g\t%g\n", x, p0[i], p1[i], p2[i]);
    }
    fclose(output_file);

    return 0;
}