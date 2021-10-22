/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::LSPSootModel

Description
    > This soot model is based on fuel laminar smoke point height LSP
    > The fuel dependent parameters are the LSP and the peak soot formatoin rate
    > Dependence of peak soot formation rate on pressure is ignored

    This class includes:
    	- direct calculation of mixutre fraction Z
    	- solution of a transport equation of variance of mixture fraction Zvar_SGS
    	- solution of a transport equation of soot mass fraction Ysoot

    The input dictionary reads like in the radiationProperties dictionary:

            sootModel LSPSootModel<gasHThermoPhysics>;

            LSPSootModelCoeffs
            {
                Xtilde_resolution 1e-3;
                XVariance_resolution 1e-4;
                Xtilde_max  1.0;
                Xvariance_max   0.25;

                solveSoot           true;

                turbulence          true;
                
                oxidation           true;

                rho_soot            1740;

                smokePointHeight                0.106; 
                adiabaticFlameTemperature       2375;

                T_inf                   298;
                T_oxidizer              298;   

                rho_oxidizer            1.1839;
                MW_oxidizer             28.9647;
                
                proportionalityConst        0.4;
            }

SourceFiles
    LSPSootModel.C

\*---------------------------------------------------------------------------*/

#ifndef SPHSootModel_H
#define SPHSootModel_H

#include "fvCFD.H"
#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"

#include "fluidThermo.H"
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"

#include "LESModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"

#include "sootSource.H"
#include "ZbetaPDF.H"
#include "TbetaPDF.H"

#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class LSPSootModel Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class LSPSootModel
:
    public sootModel
{

    // Static functions

        //- Check mixture in thermo
        static const reactingMixture<ThermoType>& checkThermo
        (
            const fluidThermo&
        );


    // Private data

         //- Soot mass fraction
        volScalarField Ysoot;

        //- Soot volume fraction
        volScalarField fv;

        //- Density of the two-phase
        volScalarField rho2;

        //- Mean soot source terms
        volScalarField sootFormationRate;
        volScalarField sootOxidationRate;

        //- Explicit limiter of soot oxidation rate
        volScalarField oxidationLimiter;

        //- Soot time derivative and convection terms
        volScalarField sootConvection;
        volScalarField sootTimeDer;
        volScalarField thermophoresis;

        //- Mixture fraction
        volScalarField Z;

        //- Variance of Mixture fraction
        volScalarField Zvar_SGS;

        //- Normalized temperature
        volScalarField Tstar;

        //- Variance of normalized temperature
        volScalarField TstarVar;

        //- Mean density 
        volScalarField rhobar;

        //- Soot model dictionary
        dictionary coeffsDict_;

        // - Solve switch
        Switch solveSoot_;

        // - Solve switch
        Switch turbulence_;

        // -  Oxidation switch
        Switch oxidation_;

        // - Soot density [kg/m3]
        dimensionedScalar rho_soot;

        // - Smoke pint height
        scalar smokePointHeight;
        
        // - Adiabatic flame temperature
        dimensionedScalar T_adiabatic;

        // - Freestream temperature
        dimensionedScalar T_inf;

        // - Oxidizer temperature
        dimensionedScalar T_oxidizer;

        // - oxidizer density
        scalar rho_oxidizer;

        // - oxidizer molecular weight
        scalar MW_oxidizer;      

        // - Temperature variance proportionality const.
        scalar proportionalityConst;

        //- Thermo packages
        const fluidThermo& thermo_;
        const psiReactionThermo& thermo;

        //- Index of fuel and oxygen
        label fuelIndex;
        label O2Index;

        //- Fuel molecular wight
        scalar MW_fuel;

        //- Fuel stoich. ratio
        scalar s;

        // - Mass fraction of oxygen in oxidizer stream
        scalar YO2Inf;  

        //- Fuel mass fraction in fuel supply stream
        scalar YFInf;

        //- Stoic. mixture fraction
        scalar Z_st;
        
        //- Peak soot formation rate
        scalar peakFormationRate;

        //- Soot source class
        sootSource SS;
        
        //- Lookup table of Beta-PDF integration functions
        autoPtr<interpolation2DTable<scalar>> lookup_Fsf;
        autoPtr<interpolation2DTable<scalar>> lookup_Gsf;
        autoPtr<interpolation2DTable<scalar>> lookup_Fso;
        autoPtr<interpolation2DTable<scalar>> lookup_Gso;
        autoPtr<interpolation2DTable<scalar>> lookup_Frho;
        autoPtr<interpolation2DTable<scalar>> lookup_Grho;

public:

    //- Runtime type information
    TypeName("LSPSootModel");


    // Constructors

        //- Construct from components
        LSPSootModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~LSPSootModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();


        // Access

            //- Return Xsoot
            const volScalarField& soot() const
            {
                return Ysoot;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LSPSootModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
