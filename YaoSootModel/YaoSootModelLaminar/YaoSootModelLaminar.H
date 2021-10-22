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
    Foam::radiation::YaoSootModelLaminar

Description
    This soot model is based on fuel laminar smoke point height LSP, the model was
    originally developed by W. Yao et al.,FireSafetyJournal 46(2011)371–387.

    NOTE: This current version ignores turbulence-soot interactions

    This class includes:
    	- direct calculation of mixutre fraction Z
    	- solution of a transport equation of soot mass fraction Ysoot

    The input dictionary reads like in the radiationProperties dictionary:

            sootModel YaoSootModelLaminar<gasHThermoPhysics>;

            YaoSootModelLaminarCoeffs
            {

                solveSoot   true;
                
                rhoSoot     1800;

                //fuel specific parameters
                Af          4e-5; 
                Z_sf        0.15;
                Z_so        0.108;        
            }

SourceFiles
    YaoSootModelLaminar.C

\*---------------------------------------------------------------------------*/

#ifndef YaoSootModelLaminar_H
#define YaoSootModelLaminar_H

#include "fvCFD.H"
#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"

#include "fluidThermo.H"
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class YaoSootModelLaminar Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class YaoSootModelLaminar
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

        //- Thermo packages
        const psiReactionThermo& thermo;

         //- Soot mass fraction
        volScalarField Ysoot;

        //- Soot volume fraction
        volScalarField fv;

        //- Mean soot source terms
        volScalarField sootFormationRate;
        volScalarField sootOxidationRate;  

        //- Soot time derivative and convection terms
        volScalarField sootConvection;
        volScalarField sootTimeDer;
        volScalarField thermophoresis;

        //- Mixture fraction
        volScalarField Z;

        //- Oxygen concentration
        volScalarField O2Concentration;

        //- Soot model dictionary
        dictionary coeffsDict_;

        // - Solve switch
        Switch solveSoot;

        // - Soot density [kg/m3]
        dimensionedScalar rhoSoot;

        //- Formation pre-exponential factor
        scalar Af;

        //- Fuel stoich. ratio
        scalar s;

        //- Index of oxygen and fuel
        label O2Index;      
        label fuelIndex; 

        // - Mass fraction of oxygen in oxidizer stream
        scalar YO2Inf;  

        //- Fuel mass fraction in fuel supply stream
        scalar YFInf;

        //- Stoic. mixture fraction
        scalar Z_st;

        //- Formation and oxidation mixture frac. limits
        scalar Z_sf;
        scalar Z_so;

        //- Formation rate parameters
        scalar gamma;
        scalar Ta;

        //- Oxidation rate parameters
        scalar Asoot;
        scalar Aox;
        scalar EaOx;

public:

    //- Runtime type information
    TypeName("YaoSootModelLaminar");


    // Constructors

        //- Construct from components
        YaoSootModelLaminar
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~YaoSootModelLaminar();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();


        // Access

            //- Return soot mass fraction
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
    #include "YaoSootModelLaminar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
