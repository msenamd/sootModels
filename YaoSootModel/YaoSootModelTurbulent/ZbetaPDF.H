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

using namespace std;

#ifndef ZBETAPDF_H
#define ZBETAPDF_H

// Inputs:
// Ztilde (instantaneous value of mean mixture fraction)
// Zvar (instantaneous value of variance of mixture fracion)


//Itegration of Func(Z)*Ptilde(Z)*dZ 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double Zbeta(
			const double& Ztilde,
			const double& Zvar,
			const double& lower,
			const double& upper
			);

//The function Func(Z)
//~~~~~~~~~~~~~~~~~
double Func(const double& Z, const double& lower, const double& upper);

#endif