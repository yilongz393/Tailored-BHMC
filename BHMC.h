//Created on 08/2021 by Yilong Zhou

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <iomanip>

using namespace std;

// BASIN-HOPPING MONTE CARLO OPTIMIZATION
class BHMC{
public:
	//Parameters
	const int NP_I;                 // Number of NP of type I
	const int NP_II;                // Number of NP of type II
	const double RADI_I;            // Radius of macro-particles I (nm)
	const double RADI_II;           // Radius of macro-particles II (nm)
	const int NATOM;    			// Total Number of NPs in the system

	const double Chi;               // Ratio between surface tension difference of NP_I and NP_II with both matrix
	const double eff;               // Factor accounting for effective radius
	
	const double v1;                // Ratio between surface tension difference of NPs and interfacial tension Gamma12
	const double v2;                // Nondimensionalize the Hamaker constant by pi*RI^2*Gamma12

	const double PI;                // Value of pi  // if float pi=3.1415927
	//const double AHM;               // Hamaker constant
	const double SIGMA;             // atom size parameter (nm)
	const double HEIGHT;
	const double DENS;              // System density in atoms/A^3  // Make L>5*sigma

	const double LSIMBOX;

	//POTENTIAL, MAGNITUDE OF GRADIENT VECTORS
	double EPOT;
	double EPOTR;             // Real potential energy without modification

	double PPOT;              //POTENTIAL ENERGY DUE TO PARTICLE-PARTICLE INTERACTION
	double IPOT;              //POTENTIAL ENERGY DUE TO PARTICLE-INTERFACE INTERACTION

	// ARRAYS FOR PRE-LOCATION,FORCE
	double* RX;
	double* RY;
	double* RZ;
	double* FX;
	double* FY;
	double* FZ;

	// INITIALIZE PARAMETERS
	BHMC(int NP1, int NP2, double R1, double R2, double chi, double eff, double v1, double v2): NP_I(NP1), NP_II(NP2), RADI_I(R1), RADI_II(R2), NATOM(NP_I+NP_II), Chi(chi), eff(eff), v1(v1), v2(v2), PI(3.141592653589793238463), SIGMA(1.0), HEIGHT(10.0*RADI_I), DENS(double(NATOM)/pow(15*RADI_I,2)/HEIGHT), LSIMBOX(pow((double(NATOM)/DENS)/HEIGHT,(1.0/2.0))), EPOT(0), EPOTR(0), PPOT(0), IPOT(0), RX(new double[NATOM]), RY(new double[NATOM]), RZ(new double[NATOM]), FX(new double[NATOM]), FY(new double[NATOM]), FZ(new double[NATOM]){
		cout<<"PARAMETERS INITIALIZATION \n"<<endl;
	}

	// FUNCTION FOR INITIALIZING THE POSITIONS
	void iniconfig()
	{
		//ASSIGN DISTANCE CUTOFF TO PREVENT OVERLAPS
		double SIGMACUT2 = -0.1*2.0*RADI_I*SIGMA+2.0*RADI_I;
			   SIGMACUT2 = SIGMACUT2*SIGMACUT2;

		const double HEIGHT2 = HEIGHT/2.0;
		const double LSIMBOX2 = LSIMBOX/2.0;

		// RESET TO 1 THE SEED FOR THE RANDOM NUMBER GENERATOR
		srand(time(0));
		//srand(1);

		// ASSIGN POSITION OF FIRST ATOM.
		double POSX = ((double)rand()/(RAND_MAX));
		double POSY = ((double)rand()/(RAND_MAX));
		double POSZ = ((double)rand()/(RAND_MAX));

		//TRANSFORM THIS RANDOM NUMBER BETWEEN 0 AND 1 TO A LENGTH BETWEEN
		//-1/2*SIMULATION BOX LENGTH AND +1/2*SIMULATION BOX LENGTH WITH
		RX[0] = (LSIMBOX*POSX - LSIMBOX2);                   // Array starts from 0 (first atom)
		RY[0] = (LSIMBOX*POSY - LSIMBOX2);
		RZ[0] = (HEIGHT*POSZ - HEIGHT2);

		// MAKE SURE TWO ATOMS DO NOT OVERLAP BY MORE THAN 10% OF SIGMA
		for (int I = 1; I<=NATOM-1; I++)
		{
			bool REPEAT = true;
			while (REPEAT)
			{
				REPEAT = false;
				// ASSIGN THE POSITION OF I'TH ATOM WITHIN (-L/2,+L/2)
				POSX = (LSIMBOX*((double)rand()/(RAND_MAX)) - LSIMBOX2);
				POSY = (LSIMBOX*((double)rand()/(RAND_MAX)) - LSIMBOX2);
				POSZ = (HEIGHT*((double)rand()/(RAND_MAX)) - HEIGHT2);

				// CALCULATE DISTANCES FROM PREVIOUSLY INSTERTED ATOMS
				for (int K = 0; K<=I-1; K++)
				{
					double DRX = POSX - RX[K];
					double DRY = POSY - RY[K];
					double DRZ = POSZ - RZ[K];
						   DRX = DRX - LSIMBOX*round(DRX/LSIMBOX);
						   DRY = DRY - LSIMBOX*round(DRY/LSIMBOX);
						   DRZ = DRZ - HEIGHT*round(DRZ/HEIGHT);
					double DELRSQ = pow(DRX,2) + pow(DRY,2) + pow(DRZ,2);

			// TRY A DIFFERENT RANDOM POSITION FOR THE CURRENT ATOM IF IT IS
			// TOO CLOSE TO ANY OF THE PREVIOUSLY INSERTED ATOMS
					if (DELRSQ < SIGMACUT2)
					{
						REPEAT = true;                   //Reset the position of the atom
						break;
					}
				}
			}

			printf("INSTERTION OF MOL %4i SUCCESSFULL \n", I+1);

			RX[I] = POSX;
			RY[I] = POSY;
			RZ[I] = POSZ;
		}
	}

	// FUNCTION FOR CALCULATING THE FORCE ON EACH ATOM
	void forces(double RCX[], double RCY[], double RCZ[])
	{
		// SET FORCES ON ALL ATOMS, POTENTIAL ENERGY, AND PRESSURE TO ZERO
		for (int I=0;I<=NATOM-1;I++)
		{
			FX[I] = 0;
			FY[I] = 0;
			FZ[I] = 0;
		}

		EPOT = 0;
		EPOTR = 0;                  //Real potential without modification
		PPOT = 0;                   //Modified potential
		IPOT = 0;                   //Interfacial energy


		double ML = 3.0;             // Modify the long range potential at position ML*RADI
		double LS = 1.5;
		// LOOP OVER ALL PAIRWISE INTERACTIONS TO OBTAIN POTENTIAL ENERGY,
		// CONFIGURATIONAL PART OF PRESSURE, AND FORCES
		 for (int I = 0; I<=NATOM-2;I++)
		{
			double FXI = FX[I];
			double FYI = FY[I];
			double FZI = FZ[I];
			double RXI = RCX[I];
			double RYI = RCY[I];
			double RZI = RCZ[I];
			for (int J = I+1;J<=NATOM-1;J++)
			{

				double DETA;            //Shifted distance
				double FPM;             //Modify the repulsion part in the potential
				double BPM;             // Modify the long range part in the potential
				double ACOEF;           // Coefficient in potential and force terms

				if (I < NP_I && J < NP_I)
				{
					// Radius of bigger particles
					DETA = RADI_I + RADI_I -1;
					FPM = 2*RADI_I - 0.5;
					BPM = ML*RADI_I;
					ACOEF = 4.0*v2;    // Coefficient in potential and force terms
				}
				else if (I < NP_I && J >= NP_I )
				{
					// interaction between big and small particles
					DETA = RADI_I + RADI_II -1;
					FPM = RADI_I + RADI_II - 0.5;
					BPM = ML*RADI_I -(RADI_I - RADI_II);
					ACOEF = 4.0*v2*pow(RADI_II/RADI_I,0.5);    // Coefficient in potential and force terms

				}
				else{
					//  Radius of bigger particles
					DETA = RADI_II + RADI_II -1;
					FPM = 2*RADI_II - 0.5;
					BPM = ML*RADI_II;
					ACOEF = 4.0*v2*RADI_II/RADI_I;    // Coefficient in potential and force terms

				}

				double RXIJ = RXI-RCX[J];
				double RYIJ = RYI-RCY[J];
				double RZIJ = RZI-RCZ[J];

				//APPLY MINIMUM IMAGE CONVENTION
					   RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);   //fix-wall boundary conditions will be applied to box in x directions
					   RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);
					   //RZIJ = RZIJ - HEIGHT*round(RZIJ/HEIGHT);
				double RIJSQ = pow(RXIJ,2) + pow(RYIJ,2) + pow(RZIJ,2);
				double RIJMS = pow(RIJSQ,0.5);
					   //DTEST=RIJMS;

				// Storing the surface-surface distance to Array RIJ
				//       RIJ[ns]= RIJMS-2.0*RADI_I;

				//IF DISTANCE IS LESS THAN LJ CUTOFF, CALCULATE FORCE ACCORDING TO
				//GRADIENT OF POTENTIAL ENERGY AND THE CONFIGURATIONAL COMPONENT OF PRESSURE BY SUMMING OVER RIJ*FIJ

				double FIJ = 0;
				double FXIJ  = 0;
				double FYIJ  = 0;
				double FZIJ  = 0;
				if (RIJMS < FPM)               // Use linear function to modify the potential when the region of particle overlap > 0.5
				{
					double RIJLF = FPM;
					double R6I = pow(SIGMA,6)/pow(RIJLF-DETA,6);
					double R12I = R6I*R6I;
					double A = (-12.0*R12I+6.0*R6I)/(RIJLF-DETA);
					double B = (-R6I + R12I) - A*RIJLF;

					PPOT = PPOT + ACOEF*(A*RIJMS + B);
					FIJ = -A;
					FXIJ  = RXIJ/RIJMS * FIJ;
					FYIJ  = RYIJ/RIJMS * FIJ;
					FZIJ  = RZIJ/RIJMS * FIJ;

					EPOTR = EPOTR + ACOEF*(A*RIJMS + B);

				}
				else if (RIJMS > BPM){
					double RIJLF = BPM;
					double R6I = pow(SIGMA,6)/pow(RIJLF-DETA,6);
					double R12I = R6I*R6I;
					double A = LS*(-12.0*R12I+6.0*R6I)/(RIJLF-DETA);              //Modify the SLOPE by LS times
					double B = (-R6I + R12I) - A*RIJLF;

					PPOT = PPOT + ACOEF*(A*RIJMS + B);
					FIJ = -A;
					FXIJ  = RXIJ/RIJMS * FIJ;
					FYIJ  = RYIJ/RIJMS * FIJ;
					FZIJ  = RZIJ/RIJMS * FIJ;

					double R6IR = pow(SIGMA,6)/pow(RIJMS-DETA,6);
					double R12IR = R6IR*R6IR;
						   EPOTR = EPOTR + ACOEF*(-R6IR + R12IR);
				}
				else {
					double R6I = pow(SIGMA,6)/pow(RIJMS-DETA,6);
					double R12I = R6I*R6I;
						   PPOT = PPOT + ACOEF*(-R6I + R12I);
					FIJ = -(-12.0*R12I+6.0*R6I)/(RIJMS-DETA);
					FXIJ  = RXIJ/RIJMS * FIJ;
					FYIJ  = RYIJ/RIJMS * FIJ;
					FZIJ  = RZIJ/RIJMS * FIJ;

					EPOTR = EPOTR + ACOEF*(-R6I + R12I);
				}
				FXI   = FXI + ACOEF*FXIJ;
				FYI   = FYI + ACOEF*FYIJ;
				FZI   = FZI + ACOEF*FZIJ;
				FX[J] = FX[J] - ACOEF*FXIJ;
				FY[J] = FY[J] - ACOEF*FYIJ;
				FZ[J] = FZ[J] - ACOEF*FZIJ;
			}
			FX[I] = FXI;
			FY[I] = FYI;
			FZ[I] = FZI;
		}

		double* FINZ= new double[NATOM];
		//  MULTIPLY BACK THE M const double STPP = 1.35;                    // Surface tension at interfacesISSING 4*EPS FACTOR IN THE COMPUTED FORCES
		for (int I = 0; I<=NATOM-1;I++)
		{
			char GRAFT;
			//double GDEN;
			
			//Potential and forces with grafting density of 0.2
			if (I<NP_I)                                               //NPs who like polymer at z<0
			{
				GRAFT = 'A';
				//GDEN = 0.05;
			}
			else if (I >= NP_I)                                         //NPs who like polymer at z>0
			{
				GRAFT = 'B';
				//GDEN = 0.05;
			}
			else
			{
				GRAFT = 'N';
				//GDEN = 0.0;
			}

			// Potential and forces at interfaces (based on free energy theory)
			if (GRAFT == 'A' || GRAFT == 'N')                 // Grafting polymer may like matrix A
			{
				double REff= RADI_I*eff;               //Effective RADI_I
				if (RCZ[I] >= REff)
				{
					//double dint = REff;
					FINZ[I]=0;
					IPOT = IPOT + 4.0*v1;
				}
				else if (abs(RCZ[I]) < REff)
				{
					double dint=RCZ[I];                        // distance from interface unit(nm)
					FINZ[I]= 2.0*dint/REff/REff + 2.0*v1/REff;
					IPOT = IPOT -(1-pow(dint/REff,2))+2.0*(1+dint/REff)*v1;
				}
				else
				{
					FINZ[I]=0;
					IPOT = IPOT + 0;
				}
			}
			else if (GRAFT == 'B')                 // Grafting polymer may like matrix B
			{
				double REff = RADI_II*eff;
				if (RCZ[I] <= - REff)
				{
					//double dint = REff;
					FINZ[I]=0;
					IPOT = IPOT + 4.0*v1*Chi*pow(RADI_II/RADI_I,2.0);
				}
				else if (abs(RCZ[I]) < REff)
				{
					double dint= -RCZ[I];                        // distance from interface unit(nm)
					FINZ[I]= -(2.0*dint/REff/REff + 2.0*v1*Chi/REff)*pow(RADI_II/RADI_I,2.0);
					IPOT = IPOT +(- (1-pow(dint/REff,2))+2.0*(1+dint/REff)*v1*Chi)*pow(RADI_II/RADI_I,2.0);
				}
				else
				{
					FINZ[I]=0;
					IPOT = IPOT + 0;
				}
			}

			FX[I] = FX[I];
			FY[I] = FY[I];
			FZ[I] = (FZ[I] - FINZ[I]);
		}

		EPOT = (PPOT + IPOT);
		EPOTR = (EPOTR + IPOT);
		
		delete[] FINZ;

	}

	//FUNCTION FOR CONJUGATE GRADIENT ALGORITHM
	void CGs(){

		//ARRAYS FOR OLD GRADIENT
		double* ODGRDX = new double[NATOM];
		double* ODGRDY = new double[NATOM];
		double* ODGRDZ = new double[NATOM];

		// GRADIENT SQUARE
		double GRADSQ = 0;

		// SEARCH DIRECTION VECTORS
		double* DIRX = new double[NATOM];
		double* DIRY = new double[NATOM];
		double* DIRZ = new double[NATOM];

		  // FUNCTION FOR CALCUALTING FORCES AND POTENTIAL
		forces(RX, RY, RZ);

		double oldPotn = EPOT;
		double newPotn = EPOT + 1;
		// OBTAIN SERACH DIRECTIONS
		for (int I=0;I<=NATOM-1;I++)
		{
			DIRX[I] = FX[I];
			DIRY[I] = FY[I];
			DIRZ[I] = FZ[I];

			// GRADIENTS OF POTENTIAL
			ODGRDX[I] = - FX[I];
			ODGRDY[I] = - FY[I];
			ODGRDZ[I] = - FZ[I];
			GRADSQ = GRADSQ + ODGRDX[I]*ODGRDX[I] + ODGRDY[I]*ODGRDY[I] +ODGRDZ[I]*ODGRDZ[I];
		}

		cout<<"\n"<<"STARTING CG ITERATION \n"<<endl;

		//TOLERANCE FOR GRADIENT
		const double eps = 1e-8;     // FOR LOCAL MINIMUN
		int iter = 0;
		const int numIter = 15000;    // Max iteration number

		double odgrad = GRADSQ;
		double negrad = GRADSQ;
		while (negrad > eps*eps && iter < numIter && abs((oldPotn - newPotn)/oldPotn)> eps)                                //if tolerance > 1e-7, then x< 15
		{
			// call line search function to find minimum along search direction within the wolffe conditions
			lineSearch(RX, RY, RZ, -odgrad, DIRX, DIRY, DIRZ);

			GRADSQ = 0;
			double ODNEWG = 0;  // production of old and new gradient

			for (int I=0;I<=NATOM-1;I++)
			{
				// GRADIENTS OF POTENTIAL
				ODNEWG = ODNEWG + ODGRDX[I]*FX[I] + ODGRDY[I]*FY[I] + ODGRDZ[I]*FZ[I];
				GRADSQ = GRADSQ + FX[I]*FX[I] + FY[I]*FY[I] + FZ[I]*FZ[I];
			}

			negrad = GRADSQ;
			for (int I=0;I<=NATOM-1;I++)
			{
				// Fletcher-Reeves algorithm
				//double gama = negrad/odgrad;

				// Polak-Ribierecollaborators
				double gama = (negrad + ODNEWG)/odgrad;

				DIRX[I] = FX[I] + gama*DIRX[I];
				DIRY[I] = FY[I] + gama*DIRY[I];
				DIRZ[I] = FZ[I] + gama*DIRZ[I];
			}
			odgrad = negrad;

			for (int I=0;I<=NATOM-1;I++)
			{
				// GRADIENTS OF POTENTIAL
				ODGRDX[I] = - FX[I];
				ODGRDY[I] = - FY[I];
				ODGRDZ[I] = - FZ[I];
			}

			oldPotn = newPotn;
			newPotn = EPOT;

			iter = iter + 1;
			if (iter%1000==0)
				printf("CG ITERATION %6i \n", iter);

		}

		delete[] ODGRDX;
		delete[] ODGRDY;
		delete[] ODGRDZ;
		delete[] DIRX;
		delete[] DIRY;
		delete[] DIRZ;
	}

	// Based on strong Wolfe conditions to constrain the choice of alpha
	void lineSearch(double xc[], double yc[], double zc[], double locgrad, double DIRX[], double DIRY[], double DIRZ[])
	{
		double* lx = new double[NATOM];
		double* ly = new double[NATOM];
		double* lz = new double[NATOM];

		for (int I=0;I<=NATOM-1;I++)
		{
			lx[I]= xc[I];
			ly[I]= yc[I];
			lz[I]= zc[I];
		}

			//parameter for sufficient decrease condition
			double c1 = 0.001;
			// parameter for curvature condition
			double c2 = 0.1;
			if (c1 > c2)
			{
				printf("******** WARNING: c1 > c2 ******** \n");
				printf("c1 = %5.3f, c2 = %5.3f \n", c1, c2);
				//printf("******** STOPPING  SIMULATION ******** \n");
			}

			double alphaMax = 1000.0;
			double alpha    = 0.1;

			double alpha_0  = 0.0;
			double alpha_1  = alpha;

			double of   = EPOT;
			double of_x = of;
			double of_0 = of;
			int    iter = 0;

			double* xts = new double[NATOM];
			double* yts = new double[NATOM];
			double* zts = new double[NATOM];

			double slope0;         // original slope
			double slopec;         // new test slope

			while (1)
			{
				slope0 = locgrad;

				// Calculating the new test position
				for (int I = 0; I<=NATOM-1;I++)
				{
					xts[I] = lx[I] + alpha_1*DIRX[I];
					yts[I] = ly[I] + alpha_1*DIRY[I];
					zts[I] = lz[I] + alpha_1*DIRZ[I];
				}

				// Calculating potential and forces for new test positions
				forces(xts, yts, zts);

				of     = EPOT;
				slopec = 0;
				for (int I = 0; I<=NATOM-1;I++)
				{
					slopec = slopec - FX[I]*DIRX[I] - FY[I]*DIRY[I] - FZ[I]*DIRZ[I];
				}

				// check if current iterate violates sufficient decrease
				if  ((of > of_0 + slope0*c1*alpha_1) || ((of >= of_x ) && (iter > 0)))
				{   // there has to be an acceptable point between alpha_0 and alpha_1
					// (because c1 < c2)
				//    cout<<"large"<<endl;
					alpha = ZOOM(slope0, alpha_0, alpha_1, of_x, of_0, c1, c2, lx, ly, lz, DIRX, DIRY, DIRZ);
					break;
				}
				//current iterate has sufficient decrease, but are we too close?
				if(abs(slopec) <= -c2*slope0)
				{   //strong wolfe full filled, quit
					alpha = alpha_1;
					break;
				}

				//are we behind the minimum?
				if (slopec >= 0)
				{   //there has to be an acceptable point between alpha_0 and alpha_1
					alpha = ZOOM(slope0, alpha_1, alpha_0, of, of_0, c1, c2, lx, ly, lz, DIRX, DIRY, DIRZ);
					break;
				}

				//cout<<"small"<<endl;

				alpha_0 = alpha_1;
				alpha_1 = min(alphaMax,alpha_1*3.0);
				of_x = of;

				if (alpha_0 == alpha_1)
				{
					//printf("******** ERROR: alphaHi = alphaLo ******** \n");
					//alpha = alpha_0;
					break;
				}

				iter = iter + 1;

			}

			// cout<<"\n"<<"FIND A VALUE OF ALPHA"<<endl;

			for (int I = 0; I<=NATOM-1;I++)
			{
				lx[I] = lx[I] + alpha*DIRX[I];
				ly[I] = ly[I] + alpha*DIRY[I];
				lz[I] = lz[I] + alpha*DIRZ[I];
				lx[I] = lx[I] - LSIMBOX*round(lx[I]/LSIMBOX);
				ly[I] = ly[I] - LSIMBOX*round(ly[I]/LSIMBOX);
				lz[I] = lz[I] - HEIGHT*round(lz[I]/HEIGHT);
			}

			//forces(lx, ly, lz);

			slopec = 0;
			for (int I = 0; I<=NATOM-1;I++)
			{
				slopec = slopec - FX[I]*DIRX[I] - FY[I]*DIRY[I] - FZ[I]*DIRZ[I];
			}

			locgrad = -slopec;


		for (int I = 0; I<=NATOM-1;I++)
		{
			RX[I] = lx[I];
			RY[I] = ly[I];
			RZ[I] = lz[I];
		}

		delete[] lx;
		delete[] ly;
		delete[] lz;
		delete[] xts;
		delete[] yts;
		delete[] zts;

	}

	//Zoom called by lineSearch()
	double ZOOM(double slope0, double alphaLo, double alphaHi, double ofLo, double of_0, double c1, double c2, double lx[], double ly[], double lz[],double DIRX[], double DIRY[], double DIRZ[])
	{
		while (1)
		{
			double alpha = (alphaLo+alphaHi)/2.0;

			double of;
			double* xts = new double[NATOM];
                        double* yts = new double[NATOM];
                        double* zts = new double[NATOM];
		     
			// Calculating the new test position
			for (int I = 0; I<=NATOM-1;I++)
			{
				xts[I] = lx[I] + alpha*DIRX[I];
				yts[I] = ly[I] + alpha*DIRY[I];
				zts[I] = lz[I] + alpha*DIRZ[I];
			}

			// Calculating potential and forces for new test positions
			forces(xts, yts, zts);
			of     = EPOT;

			double slopec = 0;
			for (int I = 0; I<=NATOM-1;I++)
			{
				slopec = slopec - FX[I]*DIRX[I] - FY[I]*DIRY[I] - FZ[I]*DIRZ[I];
			}


			if (of > of_0 + c1*alpha*slope0 || of >= ofLo)
			{
				// if we do not observe sufficient decrease in point alpha, we set
				// the maximum of the feasible interval to alpha
				alphaHi = alpha;
				//cout<<" still large"<<endl;
			}
			else
			{

				if (abs(slopec) <= -c2*slope0)
				{
					return alpha;
				}
				else if (slopec*(alphaHi-alphaLo) >= 0 ) // if slope positive and alphaHi > alphaLo
				{
					alphaHi = alphaLo;
					alphaLo = alpha;
					ofLo    = of;
					//cout<<"positive"<<endl;
				}
				else
				{
					alphaLo = alpha;
					ofLo    = of;
					//cout<<"negative"<<endl;
				}
			}

			 if (abs(alphaHi-alphaLo) < 1e-8)
			{
				//printf("******** ERROR: alphaHi = alphaLo ******** \n");
				//printf("alphaHi = %5.3f, alphaLo = %5.3f \n", alphaHi, alphaLo);
				return alpha;
			}
			
			 delete[] xts;
			 delete[] yts;
			 delete[] zts;
		}

	}

	// Random move of all NPs
	void Randmove(double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		 //reset position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		 
		 double stepsizes = 0.4;

		 for (int I = 0; I< NATOM; I++){         
			RX[I] = RCX[I] + (2*((double)rand()/(RAND_MAX))-1)*stepsizes;
			RY[I] = RCY[I] + (2*((double)rand()/(RAND_MAX))-1)*stepsizes;
			RZ[I] = RCZ[I] + (2*((double)rand()/(RAND_MAX))-1)*stepsizes;

			//APPLY MINIMUM IMAGE CONVENTION
			RX[I] = RX[I] - LSIMBOX*round(RX[I]/LSIMBOX);
			RY[I] = RY[I] - LSIMBOX*round(RY[I]/LSIMBOX);
			RZ[I] = RZ[I] - HEIGHT*round(RZ[I]/HEIGHT);
		 }
	}

	// Geometric center displacement operator
	void GeoCenDis(double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		//reset position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		double dij = 1.1;
		double alphamax = 0.52;
		double alphamin = 0.25;

		double* RI = new double[NATOM];                 // distance of atom i to centroids
		double *CT = Centerofcluster(RCX, RCY, RCZ);
		double RMAX = 0;

		for (int J = 0; J < NATOM; J++){
			double RXIJ = RCX[J]-CT[0];
			double RYIJ = RCY[J]-CT[1];
			double RZIJ = RCZ[J]-CT[2];

			//APPLY MINIMUM IMAGE CONVENTION
			RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);
			RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);
			RZIJ = RZIJ - HEIGHT*round(RZIJ/HEIGHT);

			double RIJSQ = pow(RXIJ,2) + pow(RYIJ,2) + pow(RZIJ,2);
			RI[J] = pow(RIJSQ,0.5);

			if (RI[J] > RMAX){
				RMAX = RI[J];
			}
		}

		for (int I = 0; I< NATOM; I++){
			double Theta = ((double)rand()/(RAND_MAX))*2.0*PI;             //Rotation angle theta
			double Phi = ((double)rand()/(RAND_MAX))*PI;             //Rotation angle Phi
			double stepsize = ((alphamax-alphamin)*pow(RI[I]/RMAX,2)+alphamin)*dij;

			RX[I] = RCX[I] + stepsize*cos(Theta)*sin(Phi);
			RY[I] = RCY[I] + stepsize*sin(Theta)*sin(Phi);
			RZ[I] = RCZ[I] + stepsize*cos(Phi);

			//APPLY MINIMUM IMAGE CONVENTION
			RX[I] = RX[I] - LSIMBOX*round(RX[I]/LSIMBOX);
			RY[I] = RY[I] - LSIMBOX*round(RY[I]/LSIMBOX);
			RZ[I] = RZ[I] - HEIGHT*round(RZ[I]/HEIGHT);
		 }

		delete[] CT;
		delete[] RI;

	}

	// Switch operator (switch atom i of species A and atom j of species B)
	bool SwitchOperator(double RCX[], double RCY[], double RCZ[]){
		
		//reset position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		int LI = 0;
		int UI = NP_I;

		for (int J = 0; J < NATOM; J++){
			if(J < NP_I){
				if ( RCZ[J] > RCZ[LI] ){
						LI = J;
				}
			}
			else{
				 if ( RCZ[J] < RCZ[UI] ){
					UI=J;
				}
			}
		}

		if (RCZ[LI]-RCZ[UI]> 1e-8){
			RX[LI]=RCX[UI];
			RY[LI]=RCY[UI];
			RZ[LI]=RCZ[UI];

			RX[UI]=RCX[LI];
			RY[UI]=RCY[LI];
			RZ[UI]=RCZ[LI];

			return true;
		}
		else {
			return false;
		}

	}

	// Surface Angle operator (move the farthest NPs from (CX, CY) in a plane parallel to interface)
	void SurfaceAngle(double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		//Reset Position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		int FP = 0;
		double DMAX =0;

		//The geometrical centre (centroid) of the cluster
		double * CT = Centerofcluster(RCX, RCY, RCZ);

		for (int J = 0; J < NATOM; J++){    
			double RXIJ = RCX[J]-CT[0];
			double RYIJ = RCY[J]-CT[1];

			//APPLY MINIMUM IMAGE CONVENTION
			RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);
			RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);

			double RIJSQ = pow(RXIJ,2) + pow(RYIJ,2);

			if (RIJSQ > DMAX){
				DMAX = RIJSQ;
				FP = J;
			}
		}
		DMAX = pow(DMAX,0.5);

		double RANG = ((double)rand()/(RAND_MAX))*2.0*PI;             //Rotation angle

		RX[FP] = CT[0] + DMAX*cos(RANG);
		RY[FP] = CT[1] + DMAX*sin(RANG);

		delete[] CT;

	}
	
	// Move the small NPs, whose nearest neighbors is less than 2 , in a plane parallel to the interface 
	void NNeighborMove(double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		//Reset Position 
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		//ARRAYS FOR NEAREST NEIGHBORS
		int* NN = new int[NATOM];
		for (int I=0;I<=NATOM-1;I++)
		{
			NN[I] = 0;
		}
		
		//ARRAYS FOR NEAREST NEIGHBORS THAT IS LESS THAN 2
		int* NN2 = new int[NP_II];
		int SIZE = 0;     // Actual size of NN2
		
		//Determine the number of nearest neighbors of small NPs
		for (int I = NP_I; I<=NATOM-2;I++)
		{
			double NNI = NN[I];
			double RXI = RCX[I];
			double RYI = RCY[I];
			double RZI = RCZ[I];
			for (int J = I+1;J<=NATOM-1;J++)
			{
				double RXIJ = RXI-RCX[J];
				double RYIJ = RYI-RCY[J];
				double RZIJ = RZI-RCZ[J];

				//APPLY MINIMUM IMAGE CONVENTION
					   RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);   //fix-wall boundary conditions will be applied to box in x directions
					   RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);
					   //RZIJ = RZIJ - HEIGHT*round(RZIJ/HEIGHT);
				double RIJSQ = pow(RXIJ,2) + pow(RYIJ,2) + pow(RZIJ,2);
				double RIJMS = pow(RIJSQ,0.5);

				if (RIJMS < 2.6*RADI_II)
				{
					NNI++;
					NN[J] = NN[J]+1;
				}
			}
			NN[I] = NNI;
			//cout<<NNI;
		}
		
		for (int I = NP_I; I<=NATOM-1;I++)
		{
			if (NN[I] < 2)
				{
					NN2[SIZE] = I;
					SIZE++;
				}
		}

		if (SIZE != 0)
		{
				double * CT = Centerofcluster(RCX, RCY, RCZ);
				
				//Randomly pick one atom among the atoms that has nearest neighbors less than 2
				int index = rand()%SIZE;
				int FP = NN2[index];
				
				double RXIJ = RCX[FP]-CT[0];
				double RYIJ = RCY[FP]-CT[1];

				//APPLY MINIMUM IMAGE CONVENTION
				RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);
				RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);

				double RIJSQ = pow(RXIJ,2) + pow(RYIJ,2);
				double DMAX = pow(RIJSQ,0.5)*((double)rand()/(RAND_MAX)/2+0.6);   // random distance between (0.5,1) 

				double RANG = ((double)rand()/(RAND_MAX))*2.0*PI;             //Rotation angle

				RX[FP] = CT[0] + DMAX*cos(RANG);
				RY[FP] = CT[1] + DMAX*sin(RANG);
				RZ[FP] = RZ[FP] + ((double)rand()/(RAND_MAX))*RADI_II/4;
		
				//cout<<FP<<endl;
				//cout<<"Center: "<<CT[0]<<" "<<CT[1]<<" "<<CT[2]<<endl;
				//cout<<"Postion: "<<RX[FP]<<" "<<RY[FP]<<" "<<RZ[FP]<<endl;
				
				delete[] CT;			
		}
		
		delete[] NN;
		delete[] NN2;

	}
	

	// Interior move (move the farthest particle to the range close to interface and place particle and place the NP_I in the outer circle)
	void Interiormove (double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		//reset position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		int FP = 0;                   //Index of farthest particle
		double MFDI = 0;              //MAXIMUM Distance from interface
		double SIN = v1*RADI_I;         //Shifted interface
		double SINS = v1*RADI_II;

		//The geometrical centre (centroid) of the cluster
		double * CT = Centerofcluster(RCX, RCY, RCZ);

		double CX = CT[0];
		double CY = CT[1];

		for (int J = 0; J < NATOM; J++){

			if(J < NP_I){
				double dip = abs (RCZ[J] + SIN);         //distance between particle and shifted interface
				if ( dip > MFDI ){
					MFDI = dip;
					FP = J;
				}
			}
			else{
				double dip = abs (RCZ[J] - SINS);
				 if ( dip > MFDI ){
					MFDI = dip;
					FP = J;
				}
			}
		}

		double RMS =0;                   //Root mean square distance with respect to cluster center
		// Radius of gyration
		for(int J = 0; J < NATOM; J++){
			double RXIC = RCX[J]-CX;
			double RYIC = RCY[J]-CY;
				   RMS = RMS + RXIC*RXIC + RYIC*RYIC;
		}
		RMS = RMS / NATOM;
		RMS = 2.0*pow(RMS, 0.5);

		double RANG = ((double)rand()/(RAND_MAX))*2.0*PI;             //Rotation angle

		RX[FP] = CX + RMS*cos(RANG);
		RY[FP] = CY + RMS*sin(RANG);

		if (FP < NP_I ){
			 RZ[FP] = -SIN + (2*((double)rand()/(RAND_MAX))-1)*RADI_I*0.3;
		}
		else{
			RZ[FP] = SINS + (2*((double)rand()/(RAND_MAX))-1)*RADI_I*0.3;
		}

		delete[] CT;

	}

	//Twist operator (Rotating an random angle of NPs above or below the interface)
	void Twistmove(double RCX[], double RCY[], double RCZ[]){
		srand(time(0));
		
		//reset position
		for (int I = 0; I < NATOM; I++){
			RX[I] = RCX[I]; 
			RY[I] = RCY[I];
			RZ[I] = RCZ[I];
		}
		
		double RANG = ((double)rand()/(RAND_MAX))*2.0*PI;             //Rotation angle
		double CTHTA = cos(RANG);
		double STHTA = sin(RANG);

	 //The geometrical centre (centroid) of the cluster
		double * CT = Centerofcluster(RCX, RCY, RCZ);

		double UD = ((double)rand()/(RAND_MAX));        //Probability to rotate NPs at either UP or Down interface
		for(int J=0; J < NATOM; J++){
		   if(RCZ[J] < 0 && UD < 0.5){
				double RXS = RCX[J] - CT[0];
				double RYS = RCY[J] - CT[1];

				double RXP = RXS*CTHTA - RYS*STHTA;
				double RYP = RXS*STHTA + RYS*CTHTA;

				RX[J] = RXP + CT[0];
				RY[J] = RYP + CT[1];
		   }

		   if(RCZ[J] >= 0 && UD >= 0.5){
				double RXS = RCX[J] - CT[0];
				double RYS = RCY[J] - CT[1];

				double RXP = RXS*CTHTA - RYS*STHTA;
				double RYP = RXS*STHTA + RYS*CTHTA;

				RX[J] = RXP + CT[0];
				RY[J] = RYP + CT[1];
		   }
		}

		delete[] CT;

	}

	//Similarity score of two structures, 1 means maximum similarity, 0 minimum similarity
	double USR(double MTO[], double MTN[]){
		double SSF = 0;                                        //Similarity score function
		for (int i=0; i < 12; i++){
			double SQ = pow(MTO[i] - MTN[i], 2);
			SSF = SSF + SQ;
		}
		SSF = pow(SSF, 0.5)/12.0 + 1.0;
		SSF = 1.0/SSF;

		return SSF;
	}

	//Using Ultrafast shape recognition algorithm to determine duplicated structures
	//The three moments with 12 components
	double *MOMENTS(double RCX[], double RCY[], double RCZ[]){
		double *Md = new double[12];   //Matrix for the 12 descriptors

		//The geometrical centre (centroid) of the cluster
		double * CT = Centerofcluster(RCX, RCY, RCZ);
		double CX = CT[0];
		double CY = CT[1];
		double CZ = CT[2];

		//First moment: average atomic distance to the molecular centroid
		double FMC = 0;
		double* RICMS = new double[NATOM];
		int cst = 0;                       //closest atom to ctd (molecular centroid)
		int fct = 0;                       //the farthest atom to ctd
		double MIND = NATOM*RADI_I;
		double MAXD = 0;
		for (int I=0; I < NATOM; I++){
			double RXIC = RCX[I] - CX;
			double RYIC = RCY[I] - CY;
			double RZIC = RCZ[I] - CZ;

			double RICSQ = pow(RXIC,2) + pow(RYIC,2) + pow(RZIC,2);
				   RICMS[I] = pow(RICSQ,0.5);

			if(RICMS[I] < MIND){
				cst = I;
				MIND = RICMS[I];
			}
			if (RICMS[I] > MAXD){
				fct = I;
				MAXD = RICMS[I];
			}

			FMC = FMC + RICMS[I];
		}
		FMC = FMC / NATOM;

		//Second moment: the variance of these atomic distances about FMC
		double SMC = 0;
		for (int I=0; I < NATOM; I++){
			SMC = SMC + pow(RICMS[I]-FMC, 2);
		}
		SMC = SMC / (NATOM - 1);

		//Third moment: skewness of these atomic distances about FMC
		double TMC = 0;
		double SDV = pow(SMC, 0.5);                  //Standard deviation
		for (int I=0; I < NATOM; I++){
			TMC = TMC + pow((RICMS[I]-FMC)/SDV, 3);
		}
		TMC = TMC / NATOM;

		Md[0] = FMC;
		Md[1] = SMC;
		Md[2] = TMC;

		 //First moment: average atomic distance to the CST
		double FMCS = 0;
		double* RICSMS = new double[NATOM];
		for (int I=0; I < NATOM; I++){
			double RXICS = RCX[I] - RCX[cst];
			double RYICS = RCY[I] - RCY[cst];
			double RZICS = RCZ[I] - RCZ[cst];

			double RICSSQ = pow(RXICS,2) + pow(RYICS,2) + pow(RZICS,2);
				   RICSMS[I] = pow(RICSSQ,0.5);

			FMCS = FMCS + RICSMS[I];
		}
		FMCS = FMCS / NATOM;

		//Second moment: the variance of these atomic distances about CST
		double SMCS = 0;
		for (int I=0; I < NATOM; I++){
			SMCS = SMCS + pow(RICSMS[I]-FMCS, 2);
		}
		SMCS = SMCS / (NATOM - 1);

		//Third moment: skewness of these atomic distances about CST
		double TMCS = 0;
		double SDVS = pow(SMCS, 0.5);                  //Standard deviation
		for (int I=0; I < NATOM; I++){
			TMCS = TMCS + pow((RICSMS[I]-FMCS)/SDVS, 3);
		}
		TMCS = TMCS / NATOM;

		Md[3] = FMCS;
		Md[4] = SMCS;
		Md[5] = TMCS;

		 //First moment: average atomic distance to the FCT
		double FMFC = 0;
		double* RIFCMS= new double[NATOM];
		int ftf = 0;
		double FARD = 0;
		for (int I=0; I < NATOM; I++){
			double RXIFC = RCX[I] - RCX[fct];
			double RYIFC = RCY[I] - RCY[fct];
			double RZIFC = RCZ[I] - RCZ[fct];

			double RIFCSQ = pow(RXIFC,2) + pow(RYIFC,2) + pow(RZIFC,2);
				   RIFCMS[I] = pow(RIFCSQ,0.5);

			if(RIFCMS[I] > FARD){
				ftf = I;
				FARD = RIFCMS[I];
			}

			FMFC = FMFC + RIFCMS[I];
		}
		FMFC = FMFC / NATOM;

		//Second moment: the variance of these atomic distances about FCT
		double SMFC= 0;
		for (int I=0; I < NATOM; I++){
			SMFC = SMFC+ pow(RIFCMS[I]-FMFC, 2);
		}
		SMFC = SMFC / (NATOM - 1);

		//Third moment: skewness of these atomic distances about FCT
		double TMFC = 0;
		double SDVF = pow(SMFC, 0.5);                  //Standard deviation
		for (int I=0; I < NATOM; I++){
			TMFC = TMFC + pow((RIFCMS[I]-FMFC)/SDVF, 3);
		}
		TMFC = TMFC / NATOM;

		Md[6] = FMFC;
		Md[7] = SMFC;
		Md[8] = TMFC;

	  //First moment: average atomic distance to the FTF
		double FMFT = 0;
		double* RIFTMS = new double[NATOM];
		for (int I=0; I < NATOM; I++){
			double RXIFT = RCX[I] - RCX[ftf];
			double RYIFT = RCY[I] - RCY[ftf];
			double RZIFT = RCZ[I] - RCZ[ftf];

			double RIFTSQ = pow(RXIFT,2) + pow(RYIFT,2) + pow(RZIFT,2);
				   RIFTMS[I] = pow(RIFTSQ,0.5);

			FMFT= FMFT + RIFTMS[I];
		}
		FMFT = FMFT / NATOM;

		//Second moment: the variance of these atomic distances about FTF
		double SMFT= 0;
		for (int I=0; I < NATOM; I++){
			SMFT = SMFT+ pow(RIFTMS[I]-FMFT, 2);
		}
		SMFT = SMFT / (NATOM - 1);

		//Third moment: skewness of these atomic distances about FTF
		double TMFT = 0;
		double SDVFT = pow(SMFT, 0.5);                  //Standard deviation
		for (int I=0; I < NATOM; I++){
			TMFT = TMFT + pow((RIFTMS[I]-FMFT)/SDVFT, 3);
		}
		TMFT = TMFT / NATOM;

		Md[9] = FMFT;
		Md[10] = SMFT;
		Md[11] = TMFT;

		delete[] CT;
		delete[] RICMS;
		delete[] RICSMS;
		delete[] RIFCMS;
		delete[] RIFTMS;

		return Md;
	}

	//Function to calculate centroids of clusters
	double * Centerofcluster(double RCX[], double RCY[], double RCZ[]){
		double * Ctds = new double[3];                  // array to store coordinate of Centroids

		//The geometrical centre (centroid) of the cluster
		double CX = RCX[0];
		double CY = RCY[0];
		double CZ = RCZ[0];
		for (int J = 1; J < NATOM; J++){
			double RXIJ = RCX[J]-RCX[0];
			double RYIJ = RCY[J]-RCY[0];
			double RZIJ = RCZ[J]-RCZ[0];

			//APPLY MINIMUM IMAGE CONVENTION
					RXIJ = RXIJ - LSIMBOX*round(RXIJ/LSIMBOX);
					RYIJ = RYIJ - LSIMBOX*round(RYIJ/LSIMBOX);
					RZIJ = RZIJ - HEIGHT*round(RZIJ/HEIGHT);

					RCX[J] = RCX[0] + RXIJ;
					RCY[J] = RCY[0] + RYIJ;
					RCZ[J] = RCZ[0] + RZIJ;

					CX = CX + RCX[J];
					CY = CY + RCY[J];
					CZ = CZ + RCZ[J];
		}
		CX = CX / NATOM;
		CY = CY / NATOM;
		CZ = CZ / NATOM;

		Ctds[0] = CX;
		Ctds[1] = CY;
		Ctds[2] = CZ;

		return Ctds;
	}


	~BHMC(){
		//delete[] RX;
		//delete[] RY;
		//delete[] RZ;
		delete[] FX;
		delete[] FY;
		delete[] FZ;
	}

};
