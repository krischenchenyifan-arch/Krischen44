/* riemann.h */

// Compute the analytical riemann flux for the euler equations
void CPU_Calc_Flux(double *flux, double *interface_p,
    double QL_rho, double QL_ux, double QL_vy, double QL_vz, double QL_cRT,
    double QR_rho, double QR_ux, double QR_vy, double QR_vz, double QR_cRT, double R, double GAMMA,
    double nx, double ny, double nz,
    double px, double py, double pz,
    double qx, double qy, double qz, int wall_flag);