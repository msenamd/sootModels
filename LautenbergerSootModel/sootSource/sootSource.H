/*---------------------------------------------------------------------------*\
Application
    sootSource
Description
    A function that returns the soot source terms given the mixture fraction 
    and the temperature
    
    Ref.: C.W. Lautenberger et al. / Fire Safety Journal 40 (2005) 141–176

    NOTE: For clarity, all polynomial source terms are positive
           
Author: Mohamed Ahmed
Date: January 2021
\*---------------------------------------------------------------------------*/

#ifndef SOOTSOURCE_H
#define SOOTSOURCE_H

using namespace std;

#include<vector>
#include "fvCFD.H"


/*---------------------------------------------------------------------------*\
                    Class sootSource Declaration
\*---------------------------------------------------------------------------*/
class sootSource
{

    // Private data

        //- stoichiometric mixture fraction [-]
        double Z_st;

        //- Peak soot fprmation rate [kg/m3.s]
        double f_Zp_sf;

        //- fuel molecular weight [kg/kmol]
        double MW_fuel;

        //- air molecular weight [kg/kmol]
        double MW_air;

        //- fuel adiabatic flame temperature [K]
        double T_flameAd;

    	//- free-stream temperature [K]
    	double T_inf;

    	//- reference air temperature [K]
    	double T_air;

        //- reference air density [kg/m3]
    	double rho_air;

public: 

	// Constructors
	sootSource();

  sootSource(
        double Z_st_,
        double f_Zp_sf_,
        double MW_fuel_,
        double MW_air_,
        double T_flameAd_,
        double T_inf_,
        double T_air_,
        double rho_air_
        );

  //- Destructor
  virtual ~sootSource();

  // Member Functions
  double F_sf(double Z);
  double G_sf(double Tstar);

  double F_so(double Z);
  double G_so(double Tstar);

  double F_rho(double Z);
  double G_rho(double Tstar);

  //calculate the polynomial coefficients
  void calcCoeff();

private:

  //constants of F_sf
  double psi_L_sf , psi_P_sf , psi_H_sf;
  double Z_L_sf , Z_P_sf , Z_H_sf;
  std::vector<double> c_sf; 

  //constants for G_sf
  double T_L , T_P , T_H;
  std::vector<double> k_sf; 

  //constants for G_so
  double psi_L_so , psi_P_so , psi_H_so , f_Zp_so;
  double Z_L_so , Z_P_so , Z_H_so;
  std::vector<double> c_so;

  //Auxilary function to solve 4x4 system
  std::vector<double> solveEquations(
                      double a00, double a01, double a02, double a03, double a04, 
                      double a10, double a11, double a12, double a13, double a14,
                      double a20, double a21, double a22, double a23, double a24,
                      double a30, double a31, double a32, double a33, double a34
                      );

};


#endif