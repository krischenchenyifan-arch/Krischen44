#include "riemann.h"

// Return the maximum CFL number across all cells
double CPU_Compute_Max_CFL(double *p0, double *p1, double *p2, float DX, float DT, float NO_CELLS) {
	float MAX_CFL = -1.0;
	for (int cell = 0; cell < NO_CELLS; cell++) {
		// Extract rho, u and T out from the arrays
		double rho = p0[cell];
		double u = p1[cell];
		double T = p2[cell];
		// Check the temperature
		if (T < 0) {
			printf("Error: Negative temperature in cell %d: T = %f\n. Aborting.", cell, T);
			exit(1);
		}
		double a = sqrt(1.4 * 1.0 * T); // GAMMA = 1.4, R = 1.0
		double CFL = (fabs(u) + a) * (DT / DX);
		if (CFL > MAX_CFL) {
			MAX_CFL = CFL;
		}
	}

	// Check the return value
	if (MAX_CFL < 0) {
		printf("Error: MAX_CFL is negative!. Aborting. \n");
		exit(1);
	} else {
		return MAX_CFL;
	}
}



void CPU_Calc_Flux(double *flux, double *interface_p,
    double QL_rho, double QL_ux, double QL_vy, double QL_vz, double QL_cRT,
    double QR_rho, double QR_ux, double QR_vy, double QR_vz, double QR_cRT, double R, double GAMMA,
    double nx, double ny, double nz,
    double px, double py, double pz,
    double qx, double qy, double qz, int wall_flag) {
	
	// -------- Comments by Peter Jacobs & Matthew Smith -----------------------------
	//  (i) put into the implosion (rectefm) prgram. FORTRAN90 version, June 2005
	//  (ii) Modified for use in MOCVD simulations in January 2013.
	//
	// Inviscid flux calculation based on the Riemann solver in...
	// P. A. Jacobs (1992)
	// Approximate Riemann Solver for Hypervelocity Flows.
	// A.I.A.A. Journal Vol. 30(10):2558--2561
	//
	// Version...
	// -------
	// 1.0   : 09-Mar-93 : basic isentropic solution with linear
	//                     interpolation through the expansion fans
	// 1.01  : 18-Mar-93 : add debug, constants double precision
	// 2.0   : 20-May-93 : strong-shock correction added
	// 3.0   : 18-Aug-93 : strong shock correction replaced by
	//                     general-strength shock correction
	// 3.1   : 29-Oct-93 : Added a couple more Newton steps in the 
	//                     shock-correction stage and modified the
	//                     exception handling for negative Pstar.
	//                     This should fix the non-zero mass flux
	//                     problem at boundaries and it should also
	//                     make the solver more robust (as if it wasn't
	//                     already robust).
	//
	//
	// P. A. Jacobs
	// Department of Mechanical Engineering
	// The University of Queensland
	// St Lucia, Qld 4072
	// p.jacobs@uq.edu.au
	//
	// modified for calling by MNM's 2&3d codes which allow for
	// chemical species and reactions (which may be frozen, as in this case).
	// M. N. M. May 1993.
	//
	//     Converted all references to double precision constants 1.0d0 etc
	//     to 1.0 etc.  This is better on the CRAY, I think.  Compilers will use
	//     the appropriate (8 byte) version of the constant.
	//     M. N. M. 13 Sep 1993.
	//---------------------------------------------------------------
	// 
	// may 2005, as part of the pfm TDEFM stuff
	// modified for fortran 90, chemical species removed, reduced to 2d
	// cmp input replaced by cRT = sqrt(RT) input, as for pfm
	// does NOT return interface state anymore
	// 
	//     Purpose...
	//     -------
	//     Given the flow states either side of the interface,
	//     compute the interface flow state and the fluxes across
	//     the interface.
	// 
	//     Input...
	//     -----
	//     real(dbl), intent(in) :: QL_rho,QL_ux,QL_vy,QL_cRT
	//     QL_rho  : Left density
	//     QL_ux   : Left x-velocity
	//     QL_vy   : Left y-velocity
	//     QL_cRT  : Left thermal velocity sqrt(RT)
	//     QR_rho, QR_ux, QR_vy, QR_cRT : Right state
	//     real(dbl), intent(in) :: QR_rho,QR_ux,QR_vy,QR_cRT
	// 
	//     nx, ny  : unit normals for the interface
	//     px, py  : unit normals for the interface
	//     real(dbl), intent(in) :: nx,ny,px,py
	//     Output...
	//     ------
	//     flxmss  : mass flux across the interface
	//              (mass/unit-area/unit-time)
	//     flxmnx  : flux of x-momentum across the interface
	//     flxmny  : flux of y-momentum
	//     flxeng  : flux of energy

	// DECLARE CONSTANTS
	double QL_u, QL_v, QL_w, QR_u, QR_v, QR_w;
	double QR_RT, QL_RT, QR_T, QL_T, QR_E, QL_E, QR_p, QL_p, QR_a, QL_a, QL_e, QR_e;
	double sqrR, sqrL, gm1, gp1, base,expon,pwr,z;
	double uLbar, uRbar,vacuum, term1, term2, F, geff;
	double ustar, pstar, eLstar, eRstar, rhoLstar, rhoRstar, TLstar, TRstar, aLstar, aRstar;
	double dFdpstar, delp, temporary, wspeedL, wspeedR;
	double VA, VB;
	double QI_rho, QI_u, QI_p, QI_e, QI_a, QI_T, frac, QI_v, QI_w, flxnmn, flxpmn, flxqmn, E_tot;
	double RHOMIN = 1.0e-15; 
	double EMIN = 1.0e-15;
	double AMIN = 1.0e-15;
	double BIGRAT = 1.5;
	double PMIN = EMIN*RHOMIN*(GAMMA - 1.0); 
	double CV = R/(GAMMA - 1.0);
	double TMIN = EMIN/CV; 
	double mflx, pxflx, pyflx, pzflx, eflx;
	int option;
	
	// Left hand normals
	QL_u = nx*QL_ux + ny*QL_vy + nz*QL_vz; 
	QL_v = px*QL_ux + py*QL_vy + pz*QL_vz;
	QL_w = qx*QL_ux + qy*QL_vy + qz*QL_vz;
	
	// Right hand normals
	QR_u = nx*QR_ux + ny*QR_vy + nz*QR_vz;
	QR_v = px*QR_ux + py*QR_vy + pz*QR_vz;
	QR_w = qx*QR_ux + qy*QR_vy + qz*QR_vz;
	
	// Switch normal components in case one is a wall (Reflective conditions)
	if (wall_flag == 1) {
		QR_u = -1.0*QR_u;
	} else if (wall_flag == -1) {
		QL_u = -1.0*QL_u;
	}

	// Compute the initial state properties.
	QR_RT = QR_cRT*QR_cRT;
	QL_RT = QL_cRT*QL_cRT;
	QL_T = QL_RT/R;
	QR_T = QR_RT/R;
	QL_e = QL_T*CV; 
	QR_e = QR_T*CV;
	QL_p = QL_rho*QL_RT;
	QR_p = QR_rho*QR_RT;
	QL_a = sqrt(GAMMA*QL_RT);
	QR_a = sqrt(GAMMA*QR_RT);
	sqrL = sqrt(QL_rho);
	sqrR = sqrt(QR_rho);

	geff = GAMMA; // Effective GAMMA

	gm1 = geff - 1.0; // Pretty obvious
	gp1 = geff + 1.0;

	//     -----------------------------------------------------
	//     STAGE 1: Explicit solution using two isentropic waves.
	//              This gives pstar and ustar
	//     -----------------------------------------------------
	//
	//     Intermediate variable. 

	base = QL_p/QR_p;
	expon = gm1/(2.0*geff);
	pwr = pow(base,expon);
	z = (QR_a/QL_a)*pwr;

	//   Riemann invariants. 
	uLbar = QL_u + (2.0/gm1*QL_a);
	uRbar = QR_u - (2.0/gm1*QR_a);

	vacuum = 0; // Assume no vacuum between states 
	if ((uLbar - uRbar) <= 0.0 ) {
	    // We have a situation in which a (near) vacuum is formed
	    //  between the left and right states.
	    vacuum = 1; // A near vacuum is present
	    ustar = 0.0;
	    pstar = PMIN; 
	    eLstar = EMIN;
	    eRstar = EMIN;
	    rhoLstar = RHOMIN;
	    rhoRstar = RHOMIN;
	    TLstar = TMIN;
	    TRstar = TMIN;
	    aLstar = AMIN;
	    aRstar = AMIN;
	}
	if (vacuum == 0) {
	    // Positive-pressure solution.
	    ustar = (uLbar*z + uRbar)/(1.0 + z);
	    base = (0.5*gm1*(uLbar - uRbar)/(QL_a*(1.0 + z)));
	    expon = 2.0*geff/gm1;
	    pwr = pow(base,expon);
	    pstar = QL_p*pwr;
	}

	//     ------------------------------------
	//     STAGE 2: The strong-shock correction 
	//     ------------------------------------

	if (vacuum == 0)  {
	    if ((pstar > (BIGRAT*QL_p)) && (pstar > (BIGRAT*QR_p))) {
		// Analytical solution to strong shock using shock relations
		// if both of the pressure jumps are large enough.
		// Take 4 Newton steps to get accurate estimates of pstar and ustar
		// --- Iteration 1 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p+gm1/(2.0*geff)); // check later
		term2 = sqrt(gp1/(2.0 * geff)*pstar/QR_p+gm1/(2.0*geff));
		F = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1- QR_u- QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    * ( pstar / QL_p + (3.0 * geff - 1.0) / gp1 ) 
		    / (term1 * term1 * term1)
		    - gp1 * QR_a / (4.0 * geff * geff * QR_p)
		    * ( pstar / QR_p + (3.0 * geff - 1.0) / gp1 )
		    / (term2 * term2 * term2);
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		
		//       --- Iteration 2 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//        --- Iteration 3 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1 * pstar;
		}
		//           --- Iteration 4 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p 
		    + gm1/(2.0*geff));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1)
		    /(term1*term1*term1)
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		pstar = pstar - F/dFdpstar;
		if (pstar < PMIN) {
		    pstar = PMIN;
		}
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		
		//       --- Calculate Velocity ---
		ustar = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1;
	    } else if (pstar > (BIGRAT*QR_p)) { 
		//           Treat the right-moving wave as a shock, the
		//           left-moving wave as an isentropic wave, and take
		//           four Newton steps to improve the guess for pstar, ustar.
		//           --- Iteration 1 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 2 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//          --- Iteration 3 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 4 ---
		term1 = pow((pstar/QL_p),(gm1/(2.0*geff)));
		term2 = sqrt(gp1/(2.0*geff)*pstar/QR_p
		    + gm1/(2.0*geff));
		F = QL_u
		    - 2.0*QL_a/gm1*(term1-1.0)
		    - QR_u
		    - QR_a/geff*(pstar/QR_p-1.0)/term2;
		dFdpstar = -QL_a/(geff*pstar)*term1
		    - gp1*QR_a/(4.0*geff*geff*QR_p)
		    *(pstar/QR_p+(3.0*geff-1.0)/gp1)
		    /(term2*term2*term2);
		delp = F / dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Velocity ---
		ustar = QR_u + QR_a/geff*(pstar/QR_p-1.0)/term2;
	    } else if (pstar > (BIGRAT*QL_p)) {
		//           Treat the left-moving wave as a shock, the
		//           right-moving wave as an isentropic wave, and take
		//           four Newton steps to improve the guess for pstar, ustar.
		//           --- Iteration 1 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 2 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//           --- Iteration 3 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}
		//          --- Iteration 4 ---
		term1 = sqrt(gp1/(2.0*geff)*pstar/QL_p
		    + gm1/(2.0*geff));
		term2 = pow((pstar/QR_p),(gm1/(2.0*geff)));
		F = QL_u
		    - QL_a/geff*(pstar/QL_p-1.0)/term1
		    - QR_u
		    - 2.0*QR_a/gm1*(term2-1.0);
		dFdpstar = -gp1*QL_a/(4.0*geff*geff*QL_p)
		    *(pstar/QL_p+(3.0*geff-1.0)/gp1 )
		    /(term1*term1*term1)
		    -QR_a/(geff*pstar)*term2;
		delp = F/dFdpstar;
		if (delp < pstar) {
		    pstar = pstar - delp;
		} else {
		    pstar = 0.1*pstar;
		}

		//           --- Velocity ---
		ustar = QL_u-QL_a/geff*(pstar/QL_p-1.0)/term1;
	    }
	}

	//     ---------------------------------------------
	//     STAGE 2a: Back out the other quantities using
	//               thermodynamic relations.
	//     ---------------------------------------------

	if (vacuum == 0 ) {
	    if (pstar > QL_p) {
		//   Back out values using the shock relations.
		//   Density -- from the Rankine-Hugoniot relations 
		
		rhoLstar = QL_rho*
		    (gp1*pstar+gm1*QL_p)/
		    (gp1*QL_p+gm1*pstar);
		
		//    Specific energy -- from the Equation of state 
		eLstar = pstar/(gm1*rhoLstar);
		
		//    Local speed of sound -- Perfect gas version.
		aLstar = sqrt(geff*gm1*eLstar);
	    } else {
		
		//        Use the isentropic-wave relations. 
		//        Local speed of sound -- Riemann invariants. 
		aLstar = (uLbar-ustar)*0.5*gm1;
		
		//        Specific energy -- sound speed 
		eLstar = aLstar*aLstar/(geff*gm1);
		
		//       Density -- equation of state
		rhoLstar = pstar/(gm1*eLstar);
	    }
	    if (pstar > QR_p) {
		//           Back out values using the shock relations.
		//           Density -- from the Rankine-Hugoniot relations 
		rhoRstar = QR_rho*
		    (gp1*pstar+gm1*QR_p)/
		    (gp1*QR_p+gm1*pstar);
		
		//           Specific energy -- from the Equation of state 
		eRstar = pstar/(gm1*rhoRstar);
		
		//           Local speed of sound -- Perfect gas version. 
		aRstar = sqrt(geff*gm1*eRstar);
	    } else {
		
		//           Use the isentropic-wave relations. 
		//           Local speed of sound -- Riemann invariants. 
		aRstar = (ustar-uRbar)*0.5*gm1;
		
		//           Specific energy -- sound speed
		eRstar = aRstar*aRstar/(geff*gm1);
		
		//           Density -- equation of state
		rhoRstar = pstar/(gm1*eRstar);
	    }
	    
	    //    Temperatures -- equation of state also.
	    TLstar = QL_T*(pstar*QL_rho)/(QL_p*rhoLstar);
	    TRstar = QR_T*(pstar*QR_rho)/(QR_p*rhoRstar);
	}

	//     ***********
	//     Wave speeds. 
	//     ***********

	if ((pstar > QL_p) && (vacuum == 0)) {
	    //     Left wave is a shock. 
	    temporary = 0.5*gp1*QL_p/QL_rho*(pstar/QL_p+gm1/gp1);
	    wspeedL = QL_u - sqrt(temporary);
	} else {
	    //     Left wave is an expansion fan.
	    wspeedL = QL_u - QL_a;
	}
	//getch();
	if ((pstar > QR_p) && (vacuum == 0)) {
	    //     Right wave is a shock. 
	    temporary = 0.5*gp1*QR_p/QR_rho*(pstar/QR_p+gm1/gp1);
	    wspeedR = QR_u+sqrt(temporary);
	} else {
	    //    Right wave is an expansion fan. 
	    wspeedR = QR_u + QR_a;
	}
	//getch();

	//     ************************************
	//     Decide which way the waves are going. 
	//     ************************************
	option = 0;


	if (ustar > 0.0) {
	    //        The left wave and the contact discontinuity determine
	    //        the cell interface quantities. 
	    if ((pstar - QL_p) >= 0.0) {
		//           The left wave is a shock. 
		if (wspeedL >= 0.0) {
		    //              All waves have gone into the right cell.
		    //              The values are taken from the left cell only. 
		    option = 1;
		} else {
		    //              The values are taken from behind the left shock.
		    option = 2;
		}
	    } else {
		//           The left wave is an expansion.
		VA = QL_u - QL_a;
		VB = ustar - aLstar;
		if (VA >= 0.0) {
		    //              All waves have gone into the right cell.
		    //              The values are taken from the left cell only.
		    option = 3;
		} else if (VB > 0.0) {
		    //              The left rarefaction straddles the cell interface.
		    //              Interpolate velocity inside the left rarefaction. 
		    option = 4;
		} else {
		    //              The values come from behind the left rarefaction. 
		    option = 5;
		}
	    }
	} else {
	    //        The right wave and the contact discontinuity determine
	    //        the cell interface quantities. 
	    if ((pstar - QR_p) >= 0.0) {
		//          The right wave is a shock. 
		if (wspeedR < 0.0) {
		    //              All waves have gone into the left cell.
		    //              The values are taken from the right cell only.
		    option = 11;
		} else {
		    //              The values are taken from behind the right shock.
		    option = 12;
		}
	    } else {
		//           The right wave is an expansion. 
		VA = QR_u + QR_a;
		VB = ustar + aRstar;
		if (VA <= 0.0) {
		    //              All waves have gone into the left cell.
		    //              The values are taken from the right cell only. 
		    option = 13;
		} else if (VB < 0.0) {
		    //              The right rarefaction straddles the cell interface.
		    //              Interpolate velocity inside the right rarefaction. 
		    option = 14;
		} else {
		    //              The values come from behind the right rarefaction.
		    option = 15;
		}
	    }
	}

	//     ************************************************
	//     *  Now, copy or interpolate the relevant data. *
	//     ************************************************

	if ((option == 1) || (option == 3)) {
	    //        All waves have gone into the right cell.
	    //        The values are taken from the left cell only. 
	    QI_rho = QL_rho;
	    QI_u   = QL_u;
	    QI_p   = QL_p;
	    QI_e   = QL_e;
	    QI_a   = QL_a;
	    QI_T   = QL_T;
	}

	if ((option == 2) || (option == 5)) {
	    //        The values are taken from behind the left wave.
	    QI_rho = rhoLstar;
	    QI_u   = ustar;
	    QI_p   = pstar;
	    QI_e   = eLstar;
	    QI_a   = aLstar;
	    QI_T   = TLstar;
	}

	if (option == 4) {
	    //        The left wave is an expansion. 
	    VA = QL_u - QL_a;
	    VB = ustar - aLstar;
	    //        The left rarefaction straddles the cell interface.
	    //        Interpolate velocity inside the left rarefaction. 
	    frac = (-VA)/(VB-VA);
	    QI_u = QL_u - frac * (QL_u - ustar);
	    //        Take the easy way out and linearly interpolate. 
	    QI_a   = QL_a-frac*(QL_a-aLstar);
	    QI_rho = QL_rho-frac*(QL_rho-rhoLstar);
	    QI_p   = QL_p-frac*(QL_p-pstar);
	    QI_e   = QL_e-frac*(QL_e-eLstar);
	    QI_T   = QL_T-frac*(QL_T-TLstar);
	}

	//     For the options below ...
	//     The right wave and the contact discontinuity determine
	//     the cell interface quantities. 

	if ((option == 11) || (option == 13)) {
	    //        All waves have gone into the left cell.
	    //        The values are taken from the right cell only. 
	    QI_rho = QR_rho;
	    QI_u   = QR_u;
	    QI_p   = QR_p;
	    QI_e   = QR_e;
	    QI_a   = QR_a;
	    QI_T   = QR_T;
	}

	if ((option == 12) || (option == 15)) {
	    //        The values are taken from behind the right wave. 
	    QI_rho = rhoRstar;
	    QI_u   = ustar;
	    QI_p   = pstar;
	    QI_e   = eRstar;
	    QI_a   = aRstar;
	    QI_T   = TRstar;
	}

	if (option == 14) { 
	    //        The right wave is an expansion. 
	    VA = QR_u + QR_a;
	    VB = ustar + aRstar;
	    //        The right rarefaction straddles the cell interface.
	    //        Interpolate velocity inside the right rarefaction. 
	    frac = (-VB)/(VA-VB);
	    QI_u = ustar+frac*(QR_u-ustar);
	    //        Take the easy way out and linearly interpolate. 
	    QI_a   = aRstar   + frac * (QR_a - aRstar);
	    QI_rho = rhoRstar + frac * (QR_rho - rhoRstar);
	    QI_p   = pstar   + frac * (QR_p - pstar);
	    QI_e   = eRstar   + frac * (QR_e - eRstar);
	    QI_T   = TRstar   + frac * (QR_T - TRstar);
	}


	//     ******************
	//     Passive Quantities.
	//     ******************
	//
	//     We assume that the transverse velocity is unaffected by
	//     the normal interactions.  We only need to select the
	//     correct value.

	if (QI_u < 0.0) {
	    QI_v = QR_v;
	    QI_w = QR_w;
	} else {
	    QI_v = QL_v;
	    QI_w = QL_w;
	}

	//     ******************
	//     Combine the fluxes
	//     ******************

	//    Mass/unit-area/unit-time
	mflx = QI_rho*QI_u; // Mss flux

	//    Normal momentum
	flxnmn = QI_rho*QI_u*QI_u + QI_p;

	//    Tangential momentums
	flxpmn = mflx*QI_v;
	flxqmn = mflx*QI_w;

	//    Energy Flux
	E_tot = QI_e + 0.5 * (QI_u*QI_u + QI_v*QI_v + QI_w*QI_w);
	eflx = QI_rho*E_tot*QI_u + QI_p*QI_u; // Energy flux

	//     Convert back to global axes
	pxflx = flxnmn*nx + flxpmn*px + flxqmn*qx; // Momentum fluxes
	pyflx = flxnmn*ny + flxpmn*py + flxqmn*qy;
	pzflx = flxnmn*nz + flxpmn*pz + flxqmn*qz;

	// Final Flux calculations
	flux[0] = mflx;
	flux[1] = pxflx;
	flux[2] = pyflx;
	flux[3] = pzflx;
	flux[4] = eflx;

	// States now
	interface_p[0] = QI_rho;
	interface_p[1] = QI_u;
	interface_p[2] = QI_v;
	interface_p[3] = QI_w;
	interface_p[4] = QI_T;	
}
