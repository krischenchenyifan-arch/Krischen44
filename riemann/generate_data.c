/*
    Generate Data 
    Prof. Matt Smith, NCKU
    This code does something quite simple:

    The SHLL 

*/

#include "riemann.h"

#define DEBUG 0

// Thermodynamic constants
#define GAMMA 1.4
#define R 1.0
#define CV (R/(GAMMA - 1.0))

#define NO_SAMPLES 10000

#define MIN_MACH -1.0
#define MAX_MACH 1.0
#define MIN_RHO 1e-6
#define MAX_RHO 2.0
#define MIN_T 1e-6
#define MAX_T 2.0
#define R 1.0
#define GAMMA 1.4

// All generated fluxes are 1D
// This means all of the cell interface normals are exactly the same
double nx = 1.0; double ny = 0.0; double nz = 0.0;
double px = 0.0; double py = 1.0; double pz = 0.0;
double qx = 0.0; double qy = 0.0; double qz = 1.0;

int main() {

    // Now to write the speed data to file
    FILE *fp = fopen("samples.dat", "w");

    for (int sample = 0; sample < NO_SAMPLES; sample++) {
        if (DEBUG) printf("Generating sample %d/%d\n", sample + 1, NO_SAMPLES);
        float QL_rho, QL_vx, QL_cRT;
        float QR_rho, QR_vx, QR_cRT;
        
        // Randomly generate left state
        QL_rho = MIN_RHO + ((float)rand() / RAND_MAX) * (MAX_RHO - MIN_RHO);
        float QL_Mach = MIN_MACH + ((float)rand() / RAND_MAX) * (MAX_MACH - MIN_MACH);
        QL_cRT = sqrt(R*(MIN_T + ((float)rand() / RAND_MAX) * (MAX_T - MIN_T)));
        QL_vx = QL_Mach * QL_cRT * sqrt(GAMMA);
        // Randomly generate left state
        QR_rho = MIN_RHO + ((float)rand() / RAND_MAX) * (MAX_RHO - MIN_RHO);
        float QR_Mach = MIN_MACH + ((float)rand() / RAND_MAX) * (MAX_MACH - MIN_MACH);
        QR_cRT = sqrt(R*(MIN_T + ((float)rand() / RAND_MAX) * (MAX_T - MIN_T)));
        QR_vx = QR_Mach * QR_cRT * sqrt(GAMMA);

        // Compute the flux at this interface
        double flux[5]; // The solver is a 3D solver, though we'll only use 1D here
        double interface_p[5];
        CPU_Calc_Flux(flux, interface_p,
            QL_rho, QL_vx, 0.0, 0.0, QL_cRT,
            QR_rho, QR_vx, 0.0, 0.0, QR_cRT, R, GAMMA,
            nx, ny, nz,
            px, py, pz,
            qx, qy, qz, 0);

        if (DEBUG) printf("  Left state: rho = %g, vx = %g, T = %g\n", QL_rho, QL_vx, (QL_cRT * QL_cRT) / R);
        if (DEBUG) printf("  Right state: rho = %g, vx = %g, T = %g\n", QR_rho, QR_vx, (QR_cRT * QR_cRT) / R);
        if (DEBUG) printf("  Computed fluxes: mass = %g, mom = %g, energy = %g\n", flux[0], flux[1], flux[4]);

        // Now to compute the characteristic speeds for verification
        double FLUXES_L[3];
        double FLUXES_R[3];
        double U_L[3];
        double U_R[3];

        // Compute conserved quantities
        U_L[0] = QL_rho;
        U_L[1] = QL_rho * QL_vx;
        U_L[2] = QL_rho * (((CV*QL_cRT*QL_cRT)/R) + 0.5 * QL_vx * QL_vx);
        U_R[0] = QR_rho;
        U_R[1] = QR_rho * QR_vx;
        U_R[2] = QR_rho * (((CV*QR_cRT*QR_cRT)/R) + 0.5 * QR_vx * QR_vx);

        FLUXES_L[0] = QL_rho * QL_vx; // Mass flux
        FLUXES_L[1] = QL_rho * QL_vx * QL_vx + QL_rho * (QL_cRT * QL_cRT); // Momentum flux
        FLUXES_L[2] = QL_vx*( U_L[2] + QL_rho * (QL_cRT * QL_cRT) ); // Energy flux

        FLUXES_R[0] = QR_rho * QR_vx; // Mass flux
        FLUXES_R[1] = QR_rho * QR_vx * QR_vx + QR_rho * (QR_cRT * QR_cRT); // Momentum flux
        FLUXES_R[2] = QR_vx*( U_R[2] + QR_rho * (QR_cRT * QR_cRT) ); // Energy flux

        // Compute the SHLL Split fluxes
        double FP[3], FM[3];
        // FP FIRST (L)
        double Z1 = 0.5*(QL_Mach + 1.0);
        FP[0] = FLUXES_L[0]*Z1;
        FP[1] = FLUXES_L[1]*Z1;
        FP[2] = FLUXES_L[2]*Z1;
        // FM FIRST (R)
        float Z3 = 0.5*(QR_Mach - 1.0);
        FM[0] = -FLUXES_R[0]*Z3;
        FM[1] = -FLUXES_R[1]*Z3;
        FM[2] = -FLUXES_R[2]*Z3;

        // double FLUX_DIFF[3];
        //FLUX_DIFF[0] = flux[0] - (0.5*(FLUXES_L[0]+FLUXES_R[0]) + 0.5*QL_Mach*FLUXES_L[0] - 0.5*QR_Mach*FLUXES_R[0]);
        //FLUX_DIFF[1] = flux[1] - (0.5*(FLUXES_L[1]+FLUXES_R[1]) + 0.5*QL_Mach*FLUXES_L[1] - 0.5*QR_Mach*FLUXES_R[1]);
        //FLUX_DIFF[2] = flux[4] - (0.5*(FLUXES_L[2]+FLUXES_R[2]) + 0.5*QL_Mach*FLUXES_L[2] - 0.5*QR_Mach*FLUXES_R[2]);

        // Compute the flux difference between the riemann flux and the SHLL flux (Without the dissipation term)
        double FLUX_DIFF[3];
        FLUX_DIFF[0] = flux[0] - (FP[0] + FM[0]);
        FLUX_DIFF[1] = flux[1] - (FP[1] + FM[1]);
        FLUX_DIFF[2] = flux[4] - (FP[2] + FM[2]);  // That was a nasty bug!


        // Now to compute the speed
        double SPEEDS[3];
        //SPEEDS[0] = FLUX_DIFF[0] / (0.5*((1.0-QR_Mach*QR_Mach)*U_R[0] - (1.0-QL_Mach*QL_Mach)*U_L[0]));
        //SPEEDS[1] = FLUX_DIFF[1] / (0.5*((1.0-QR_Mach*QR_Mach)*U_R[1] - (1.0-QL_Mach*QL_Mach)*U_L[1]));
        //SPEEDS[2] = FLUX_DIFF[2] / (0.5*((1.0-QR_Mach*QR_Mach)*U_R[2] - (1.0-QL_Mach*QL_Mach)*U_L[2]));
        SPEEDS[0] = FLUX_DIFF[0] / (0.5*(U_R[0] - U_L[0]));
        SPEEDS[1] = FLUX_DIFF[1] / (0.5*(U_R[1] - U_L[1]));
        SPEEDS[2] = FLUX_DIFF[2] / (0.5*(U_R[2] - U_L[2]));
        
        if (DEBUG) printf("  Estimated speeds: %g, %g, %g\n", SPEEDS[0], SPEEDS[1], SPEEDS[2]);

        // Now try to compute the final fluxes using these speeds as a test
        double FLUX_TEST[3];
        //FLUX_TEST[0] = FP[0] + FM[0] + SPEEDS[0]*(0.5*((1.0-QR_Mach*QR_Mach)*U_R[0] - (1.0-QL_Mach*QL_Mach)*U_L[0]));
        //FLUX_TEST[1] = FP[1] + FM[1] + SPEEDS[1]*(0.5*((1.0-QR_Mach*QR_Mach)*U_R[1] - (1.0-QL_Mach*QL_Mach)*U_L[1]));
        //FLUX_TEST[2] = FP[2] + FM[2] + SPEEDS[2]*(0.5*((1.0-QR_Mach*QR_Mach)*U_R[2] - (1.0-QL_Mach*QL_Mach)*U_L[2]));
        FLUX_TEST[0] = FP[0] + FM[0] + SPEEDS[0]*(0.5*(U_R[0] - U_L[0]));
        FLUX_TEST[1] = FP[1] + FM[1] + SPEEDS[1]*(0.5*(U_R[1] - U_L[1]));
        FLUX_TEST[2] = FP[2] + FM[2] + SPEEDS[2]*(0.5*(U_R[2] - U_L[2]));

        if (DEBUG) printf("  Reconstructed fluxes: mass = %g, mom = %g, energy = %g\n", FLUX_TEST[0], FLUX_TEST[1], FLUX_TEST[2]);
        if (DEBUG) printf("------------------------------------------------------\n");

        // Now to write to file
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
            QL_rho, QL_vx, QL_cRT,
            QR_rho, QR_vx, QR_cRT,
            SPEEDS[0], SPEEDS[1], SPEEDS[4]);
    }
    fclose(fp);

}

