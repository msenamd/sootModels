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

#ifndef TBETAPDF_H
#define TBETAPDF_H

// Inputs:
// TtildeStar (instantaneous value of normalized temperature)
// TvarStar (instantaneous value of variance of normalized temperature)


//integratio of: G_sf(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_sf(
			double TstarTilde,
			double TstarVar,
			sootSource SS
			);


//integratio of: G_so(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_so(
			double TstarTilde,
			double TstarVar,
			sootSource SS
			);

//integratio of: Frho(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_rho(
			double Ztilde,
			double Zvar,
			sootSource SS
			);

#endif