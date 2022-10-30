//Created on 07/2020 by Yilong Zhou

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string>
#include "BHMC.h"

using namespace std;

int main(int argc, char ** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: encrypt key inputFileName\n");
        return EXIT_FAILURE;
    }

    FILE * f = fopen(argv[1], "r");
    if (f == NULL) {
        perror("Could not open file");
        return EXIT_FAILURE;
    }
        
    char * line = NULL;
    size_t sz = 0;
    double params[8] = {0};
    for (int i=0; i<8;i++){
        getline(&line, &sz, f);
        params[i]= atof(line);
    }
    free(line);

    //index string for filename
    //string idx;                     
    //idx.push_back(argv[1][1]);
    //idx.push_back(argv[1][2]);

    ofstream GM_moviefile;                   // POSITIONS IN MOVIE FORMAT EVERY IANIMATE STEPS
    ofstream GM_energyfile;                  // STORES ENERGIES EVERY IE STEPS
    ofstream LM_moviefile;            // STORE LOCAL MINIMUM
    ofstream LM_energyfile;
    //ofstream testfile;
    LM_moviefile.open("LM_MOVIE.xyz");
    LM_energyfile.open("LM_ENERGIES.dat");
    GM_moviefile.open("GM_MOVIE.xyz");
    GM_energyfile.open("GM_ENERGIES.dat");
    //testfile.open("testc.xyz");

	BHMC STRC((int) params[0], (int) params[1], params[2], params[3], params[4], params[5], params[6], params[7]);                             // Object of BHMC class

	const double TEM = 0.5*pow((1-STRC.v1/1.5),0.05)*pow(STRC.v2,1.0);         //INITIAL TEMPERATURE
    double kT = TEM;
     //INITIALIZE THE POSITIONS AND VELOCITIES OF ATOMS
    cout<<"STARTING INITIALIZATION \n"<<endl;

     // FUNCTION FOR INITIALIZING POSITIONS AND VELOCITIES
    STRC.iniconfig();

    // CALL CONJUGATE GRADIENT METHOD FOR LOCAL MINIMUM
    STRC.CGs();

    //ARRAY FOR GLOBAL MINIMUM POSTIONS FOUND SO FAR
    double* GMRX = new double[STRC.NATOM];
    double* GMRY = new double[STRC.NATOM];
    double* GMRZ = new double[STRC.NATOM];

    //ENERGY GLOBAL MINIMUM
    double GMEPOT = 0;
    // ARRAYS FOR PREVIOUS MINIMUM POSITONS
    double* OLDRX = new double[STRC.NATOM];
    double* OLDRY = new double[STRC.NATOM];
    double* OLDRZ = new double[STRC.NATOM];

    // PREVIOUS LOCAL MINIMUM OF POTENTIAL
    double OLDEPOT = 0;


    const int MCSTEPS = 20000;     // MONTE CARLO STEPS

    int ATTEMPT = 0;
    int ACCEPT = 0;
    int LSTEP = ACCEPT;                // Last step

    OLDEPOT = STRC.EPOTR;
    GMEPOT = STRC.EPOTR;
    for (int I=0;I<=STRC.NATOM-1;I++)
    {
        OLDRX[I] = STRC.RX[I];
        OLDRY[I] = STRC.RY[I];
        OLDRZ[I] = STRC.RZ[I];

        //UPDATE GLOBAL MINIMUM
        GMRX[I] = STRC.RX[I];
        GMRY[I] = STRC.RY[I];
        GMRZ[I] = STRC.RZ[I];
    }

    int CONREJ = 0;           // the number of Consecutive rejections
    int SAMESTRUC = 0;        // the number of consecutively converging into one same structure
    int cindex = 0;
    bool SW = false;

    //STORING THE MOMENTS OF STRUCTURES
    double MTO[12]={0};
    double MTN[12]={0};
    double SMS = 0;    //Similarity score

    // STARTING THE MONTE CARLO
    for (int ST = 1; ST <= MCSTEPS; ST++){

        if(ST % 1 == 0){
            cout<<"\n"<<"STARTING MC STEP "<<ST<<endl;
        }

        //TAKE THE CURRENT BEST STRUCTURE AS A NEW STARTING POINT FOR EVERY CERTAIN AMOUNT OF STEPS
		if (ST % 2000 == 0){
			OLDEPOT = GMEPOT;
			for (int I=0;I<=STRC.NATOM-1;I++)
			{
				OLDRX[I] =  GMRX[I];
				OLDRY[I] = 	GMRY[I];			
				OLDRZ[I] = 	GMRZ[I];			
			}
			
			cindex = 5;
			SW = false;
			SAMESTRUC = 0;
			CONREJ = 0;
			
		}
        // if switch operator is true, do switching
        if (SW == true){
            SW = STRC.SwitchOperator(OLDRX, OLDRY, OLDRZ);
            //cout <<"\n"<<" Switch Operator"<<endl;
        }

		if(SW == true){
			SAMESTRUC = 0;
		}
		
		// MOVE PARTICLE POSITIONS
        int ci = cindex%8;
		
        // if not, do other operators
        if(SW == false){
           switch (ci){
               case 0:
                    STRC.GeoCenDis(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" GeoCenDis"<<endl;
                    break;
                case 1:
                    STRC.Twistmove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" Twistmove"<<endl;
                    break;
                case 2:
                    STRC.Randmove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" Randmove"<<endl;
                    break;
                case 3:
                    STRC.Interiormove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" Interiormove"<<endl;
                    break;
                case 4:
                    STRC.SurfaceAngle(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" SurfaceAngle"<<endl;
                    break;
				case 5:
                    STRC.NNeighborMove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" NNeighborMove"<<endl;
                    break; 
                case 6:
                    STRC.Randmove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" Randmove"<<endl;
                    break;
                case 7:
                    STRC.NNeighborMove(OLDRX, OLDRY, OLDRZ);
                    cout <<"\n"<<" NNeighborMove"<<endl;
                    break;
                }
        }
        SW = false;
        //CALL CONJUGATE GRADIENT
        STRC.CGs();

        //CHECK FOR ACCEPTANCE
        double DELTEP = STRC.EPOTR - OLDEPOT;
		//cout<<"delta_E: "<<DELTEP<<endl;
		
        ATTEMPT++;
        CONREJ++;

        if(DELTEP <= 0){
            OLDEPOT = STRC.EPOTR;
            for (int I=0;I<=STRC.NATOM-1;I++)
            {
                OLDRX[I] = STRC.RX[I];
                OLDRY[I] = STRC.RY[I];
                OLDRZ[I] = STRC.RZ[I];
            }
            ACCEPT++;
			SW = true;

            if(OLDEPOT < GMEPOT){
                GMEPOT = OLDEPOT;

                double NGMEPOT = GMEPOT;           // Normalized by PI*RADI_I*RADI_I
                GM_moviefile<<left<<STRC.NATOM<<endl;
                GM_moviefile<<left<<"R1="<<STRC.RADI_I<<", R2="<<STRC.RADI_II<<"; Chi = "<<STRC.v1<<", Epsilon = "<<STRC.v2<<"; Energy = "<<NGMEPOT<<endl;

                //WRITE OUT ENERGIES EVERY IE STEPS
                GM_energyfile<<setw(8)<<left<<ST<<setw(12)<<NGMEPOT<<endl;

                for (int I=0;I<=STRC.NATOM-1;I++)
                {
                    //UPDATE GLOBAL MINIMUM
                    GMRX[I] = STRC.RX[I];
                    GMRY[I] = STRC.RY[I];
                    GMRZ[I] = STRC.RZ[I];

                    if(I < STRC.NP_I){
                        GM_moviefile<<setw(15)<<left<<"C"<<setw(15)<<GMRX[I]<<setw(15)<<GMRY[I]<<setw(15)<<GMRZ[I]<<endl;
                    }
                    else{
                        GM_moviefile<<setw(15)<<left<<"H"<<setw(15)<<GMRX[I]<<setw(15)<<GMRY[I]<<setw(15)<<GMRZ[I]<<endl;
                    }
                }
            }
        }
        else if ( (double)rand()/(RAND_MAX) <= exp(-DELTEP/kT)){
            OLDEPOT = STRC.EPOTR;
            for (int I=0;I<=STRC.NATOM-1;I++)
            {
                OLDRX[I] = STRC.RX[I];
                OLDRY[I] = STRC.RY[I];
                OLDRZ[I] = STRC.RZ[I];
            }
            ACCEPT++;
        }

        double NEPOTR = STRC.EPOTR;        // Normalized by PI*RADI_I*RADI_I
        if (ACCEPT== LSTEP + 1){
            LSTEP = ACCEPT;
            CONREJ = 0;

            LM_energyfile<<setw(8)<<left<<ACCEPT<<setw(12)<<NEPOTR<<endl;
            LM_moviefile<<left<<STRC.NATOM<<endl;
            LM_moviefile<<left<<"R1="<<STRC.RADI_I<<", R2="<<STRC.RADI_II<<"; Chi = "<<STRC.v1<<", Epsilon = "<<STRC.v2<<"; Energy = "<<NEPOTR<<endl;
            for (int I=0; I<=STRC.NATOM-1;I++)
            {
                if(I < STRC.NP_I){
                    LM_moviefile<<setw(15)<<left<<"C"<<setw(15)<<STRC.RX[I]<<setw(15)<<STRC.RY[I]<<setw(15)<<STRC.RZ[I]<<endl;
                }
                else{
                    LM_moviefile<<setw(15)<<left<<"H"<<setw(15)<<STRC.RX[I]<<setw(15)<<STRC.RY[I]<<setw(15)<<STRC.RZ[I]<<endl;
                }
            }

            double *MD = STRC.MOMENTS(STRC.RX, STRC.RY, STRC.RZ);
            for (int ii=0; ii<12; ii++){
                MTN[ii] = MD [ii];
            }
            delete[] MD;
            SMS = STRC.USR(MTO,MTN);
            if(abs(SMS-1.0) < 3e-2){
                SAMESTRUC++;
            }
            else{
                SAMESTRUC=0;
            }

            for (int ii=0; ii<12; ii++){
                MTO[ii] = MTN [ii];
            }

        }

         //if reach a certain number of consecutive rejections, MOVE TO NEXT OPERATOR
        if (CONREJ > 6){
            cindex ++;
            CONREJ = 0;
            SAMESTRUC = 0;
        }

        // If converge to same structure, move to next operator
        if(SAMESTRUC > 5){
            cindex++;
            CONREJ = 0;
            SAMESTRUC = 0;

            /*
            testfile<<left<<STRC.NATOM<<endl;
            testfile<<left<<"R1="<<STRC.RADI_I<<", R2="<<STRC.RADI_II<<"; v1 = "<<STRC.v1<<", v2 = "<<STRC.v2<<"; Energy = "<<NEPOTR<<endl;
            for (int I=0; I<=STRC.NATOM-1;I++)
            {
                if(I< STRC.NP_I){
                    testfile<<setw(15)<<left<<"CH4"<<setw(15)<<STRC.RX[I]<<setw(15)<<STRC.RY[I]<<setw(15)<<STRC.RZ[I]<<endl;
                }
                else{
                    testfile<<setw(15)<<left<<"H2"<<setw(15)<<STRC.RX[I]<<setw(15)<<STRC.RY[I]<<setw(15)<<STRC.RZ[I]<<endl;
                }
            } */
        }
    }

    LM_moviefile.close();
    LM_energyfile.close();
    GM_moviefile.close();
    GM_energyfile.close();
    //testfile.close();

    delete[] GMRX;
    delete[] GMRY;
    delete[] GMRZ;
    delete[] OLDRX;
    delete[] OLDRY;
    delete[] OLDRZ;
    
    cout<<"\nGMEPOT IS "<<GMEPOT<<endl;
    return 0;
}
