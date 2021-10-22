/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::radiation::wsggmCassolSoot

Description

    wsggmCassolSoot calculates radiation absorption and emission coefficients
    for a mixture of CO2, H2O and soot (assuming p_w/p_c=2)

    Ref: 
    F. Cassol et al. International Journal of Heat and Mass Transfer 79 (2014) 796–806

    4 bands are used to describe CO2/H2O in addition to 1 bad for a transparent gas.
    Each of the 5 bands include two bands for soot. The total number of 
    bands is 10. The bands are described as following:
    	- Bands 0:4 are the 5 bands of CO2/H2O at FIRST band of soot
    	- Bands 5:9 are the 5 bands of CO2/H2O at SECOND band of soot 

    Bands 0 and 5 are transparent gases.

    The weighting coefficients (temperaure polynomial) of the transparent band
    listed in the input file are manually summed from the tables of Cassol et al.

    All the bands should have the 2 species and have to be entered
    in the same order. 1-mixture of CO2/H2O, 2-soot

    There are six values in the hiTcoeffs coeff. The first value is the pressure
    absorption coefficient for CO2/H2O mixture (or the volumetric fraction absorption 
    coefficient for soot). The other 5 coeff. are used in a temperature polynomial to 
    calculate the weighting parameter.

    The look Up table file should be in the constant directory.

    The range of validity is 400<T<2500

    band dictionary:

    \verbatim

    BetaSoot 4.1; //check value for different fuels

        band0
        {
            species
            {
                H2O_CO2 //transparent
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (0 0.45717 119.8030e-5 -131.3090e-8 55.4640e-11 -85.9380e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (22313.49 0.95552 -1.431e-3 9.871e-7 -3.39e-10 4.555e-14);
                }
            }
        }
        band1
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (0.192 0.05617 78.44e-5 -85.63e-8 42.46e-11 -74.4e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (22313.49 0.95552 -1.431e-3 9.871e-7 -3.39e-10 4.555e-14);
                }
            }
        }
        band2
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (1.719 0.14260 17.95e-5 -1.077e-8 -6.971e-11 17.74e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (22313.49 0.95552 -1.431e-3 9.871e-7 -3.39e-10 4.555e-14);
                }
            }
        }
        band3
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (11.370 0.13620 25.74e-5 -37.11e-8 15.70e-11 -22.67e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (22313.49 0.95552 -1.431e-3 9.871e-7 -3.39e-10 4.555e-14);
                }
            }
        }        
        band4
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (111.016 0.12220 -2.327e-5 -7.492e-8 4.275e-11 -6.608e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (22313.49 0.95552 -1.431e-3 9.871e-7 -3.39e-10 4.555e-14);
                }
            }
        }
        band5
        {
            species
            {
                H2O_CO2 //transparent
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (0 0.45717 119.8030e-5 -131.3090e-8 55.4640e-11 -85.9380e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (466624.8 0.08010 1.290e-3 -7.874e-7 2.322e-10 -3.084e-14);
                }
            }
        }
        band6
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (0.192 0.05617 78.44e-5 -85.63e-8 42.46e-11 -74.4e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (466624.8 0.08010 1.290e-3 -7.874e-7 2.322e-10 -3.084e-14);
                }
            }
        }
        band7
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (1.719 0.14260 17.95e-5 -1.077e-8 -6.971e-11 17.74e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (466624.8 0.08010 1.290e-3 -7.874e-7 2.322e-10 -3.084e-14);
                }
            }
        }
        band8
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (11.370 0.13620 25.74e-5 -37.11e-8 15.70e-11 -22.67e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (466624.8 0.08010 1.290e-3 -7.874e-7 2.322e-10 -3.084e-14);
                }
            }
        }        
        band9
        {
            species
            {
                H2O_CO2
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (111.016 0.12220 -2.327e-5 -7.492e-8 4.275e-11 -6.608e-15);
                }
                soot
                {
                    Tcommon         400.;
                    Tlow            200.;
                    Thigh           2500.;
                    invTemp         false;
                    loTcoeffs (0 0 0 0 0 0) ;
                    hiTcoeffs (466624.8 0.08010 1.290e-3 -7.874e-7 2.322e-10 -3.084e-14);
                }
            }
        } 

    \endverbatim


SourceFiles
    wsggmCassolSoot.C

\*---------------------------------------------------------------------------*/

#ifndef wsggmCassolSoot_H
#define wsggmCassolSoot_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModelBand.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"

#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H" 
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggmAbsorptionEmissionNewCassolSoot Declaration
\*---------------------------------------------------------------------------*/

class wsggmCassolSoot
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of species considered for absorptivity 
        static const int nSpecies_ = 20; //each band has 2 species

        //- Maximum number of bands
        static const int maxBands_ = 10;

        //-  Absorption coefficients
        FixedList<FixedList<absorptionCoeffs, nSpecies_>, maxBands_> coeffs_;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Hash table with species names
        HashTable<label> speciesNames_;

        //- Indices of species in the look-up table
        FixedList<label, nSpecies_> specieIndex_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Bands
        label nBands_;

        //- Pointer list of species being solved involved in the absorption
        UPtrList<volScalarField> Yj_;

        //- Access to soot coefficient from case file
        const scalar BetaSoot_;

public:

    //- Runtime type information
    TypeName("wsggmCassolSoot");


    // Constructors

        //- Construct from components
        wsggmCassolSoot(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggmCassolSoot();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

        // weighting coefficient for wsggm (reimplemented in radiativeIntensityRay)
                tmp<volScalarField> ggCoeffCont(const label bandI = 0) const;

            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

        //- Add contribution of ILambda to the total radiative intensity in
        //  direction i
        tmp<volScalarField> addIntensity
        (
            const label i,
            const volScalarField& ILambda
        ) const;

        void correct
        (
            volScalarField& a_,
            PtrList<volScalarField>& aLambda

        ) const;

        void correctNew
        (
            volScalarField& ggCoeff_,
            PtrList<volScalarField>& ggCoeffLambda

        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
