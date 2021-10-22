/*---------------------------------------------------------------------------*\
Application
    betaPDF

Description
 	Algorithm for integrating beta-PDF of Z for any variable
 	Ref.: F. Liu et al. Int. J. of Thermal Sciences 41 (2002) 763–772
           
Author: Mohamed Ahmed
Date: Jan. 2021

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cmath>

#include "fvCFD.H"
#include "sootSource.H"

using namespace std;

#ifndef ZBETAPDF_H
#define ZBETAPDF_H

// Inputs:
// Ztilde (instantaneous value of mean mixture fraction)
// Zvar (instantaneous value of variance of mixture fracion)


//integratio of: F_sf(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_sf(
			double Ztilde,
			double Zvar,
			sootSource SS
			);


//integratio of: F_so(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_so(
			double Ztilde,
			double Zvar,
			sootSource SS
			);


//integratio of: Frho(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_rho(
			double Ztilde,
			double Zvar,
			sootSource SS
			);

#endif