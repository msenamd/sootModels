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

\*---------------------------------------------------------------------------*/

#include "wsggmCassolSoot.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmCassolSoot, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmCassolSoot,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmCassolSoot::wsggmCassolSoot
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    speciesNames_(0),
    specieIndex_(label(0)),
    thermo_(mesh.lookupObject<fluidThermo>("thermophysicalProperties")),
    Yj_(nSpecies_),
    BetaSoot_(readScalar(coeffsDict_.lookup("BetaSoot")))

{
    label nBand = 0;
    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& dict = iter().dict();

        label nSpec = 0;
        const dictionary& specDicts = dict.subDict("species");
        forAllConstIter(dictionary, specDicts, iter)
        {
            const word& key = iter().keyword();
            if (nBand == 0)
            {
                speciesNames_.insert(key, nSpec);
            }
            else
            {
                if (!speciesNames_.found(key))
                {
                    FatalErrorIn
                    (
                        "Foam::radiation::wideBandAbsorptionEmission(const"
                        "dictionary& dict, const fvMesh& mesh)"
                    )   << "specie: " << key << "is not in all the bands"
                        << nl << exit(FatalError);
                }
            }

            coeffs_[nSpec][nBand].initialise(specDicts.subDict(key));

	    	Info << "WSGG: Species" << " " << nSpec << " " << key << " " <<  "Band" << " " 
	    	<< nBand << " " << "Coeffs" << nSpec << nBand << " = "<< coeffs_[nSpec][nBand].coeffs(1000) << endl;

            nSpec++;
        }
        nBand++;
    }
    nBands_ = nBand;

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmCassolSoot::~wsggmCassolSoot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmCassolSoot::aCont(const label bandI) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    const psiReactionThermo& thermo= mesh_.lookupObject<psiReactionThermo>("thermophysicalProperties");

    //access species mass fractions
    const PtrList<volScalarField>& Y = thermo.composition().Y();

    //access specie thermo data 
    const PtrList<gasHThermoPhysics> & specieThermo =
        dynamic_cast<const reactingMixture<gasHThermoPhysics>&>  (thermo).speciesData();

    // fraction volume for soot term
    const volScalarField& fv =mesh_.lookupObject<volScalarField>("fv");

    // get index of CO2 in mixture
    label indexCO2= dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["CO2"];    

    // get index of H2O in mixture
    label indexH2O= dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["H2O"];

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    volScalarField partialPressureCO2
    (
        IOobject
        (
            "partialPressureCO2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0)
    );

    volScalarField partialPressureH2O
    (
        IOobject
        (
            "partialPressureH2O",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0)
    );

    volScalarField wMean 
    (
        IOobject
        (
            "wMean",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0) // kg/kmol
    );

    // calculation of partial pressure in [atm]

    forAll(Y,specieI)
    {
        wMean+=Y[specieI]/specieThermo[specieI].W();
    }
    wMean=1/wMean;

    partialPressureCO2 = wMean*(p/101325/dimensionedScalar("unity",dimPressure,1.0))*(Y[indexCO2]/specieThermo[indexCO2].W());
    partialPressureH2O = wMean*(p/101325/dimensionedScalar("unity",dimPressure,1.0))*(Y[indexH2O]/specieThermo[indexH2O].W());


    // calculation of a in m-1 
    // species [0] is CO2/H2o mixture, species [1] is soot
    // kappa in tables = b[0] is supplied in m-1.atm-1
    forAll(a, i)
    {
        const absorptionCoeffs::coeffArray& b_wc   = coeffs_[0][bandI].coeffs(T[i]);
        const absorptionCoeffs::coeffArray& b_soot = coeffs_[1][bandI].coeffs(T[i]);

	    a[i] = b_wc[0]*(partialPressureCO2[i] + partialPressureH2O[i]) + b_soot[0]*BetaSoot_*fv[i];

    }

    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmCassolSoot::ggCoeffCont(const label bandI) const
{

    const volScalarField& T = thermo_.T();

    tmp<volScalarField> tggCoeff
    (
        new volScalarField
        (
            IOobject
            (
                "ggCoeff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("ggCoeff", dimless, 0.0)
        )
    );

    tmp<volScalarField> tggCoeffBC
    (
        new volScalarField
        (
            IOobject
            (
                "ggCoeffBC",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("ggCoeffBC", dimless, 0.0)
        )
    );

    scalarField& wsggmWeightingCoeff = tggCoeff.ref().primitiveFieldRef();

    forAll(wsggmWeightingCoeff, i)
    {
        const absorptionCoeffs::coeffArray& b_wc   = coeffs_[0][bandI].coeffs(T[i]);
        const absorptionCoeffs::coeffArray& b_soot = coeffs_[1][bandI].coeffs(T[i]);

        if( (bandI == 0) | (bandI == 5) )
        {
        	wsggmWeightingCoeff[i] = (1.0 - (b_wc[1] + b_wc[2]*T[i]+ b_wc[3]*pow(T[i],2.0) + b_wc[4]*pow(T[i],3.0) + b_wc[5]*pow(T[i],4.0))) *
                			(b_soot[1] + b_soot[2]*T[i]+ b_soot[3]*pow(T[i],2.0) + b_soot[4]*pow(T[i],3.0) + b_soot[5]*pow(T[i],4.0));
        }
        else
        {
        	wsggmWeightingCoeff[i] = (b_wc[1] + b_wc[2]*T[i]+ b_wc[3]*pow(T[i],2.0) + b_wc[4]*pow(T[i],3.0) + b_wc[5]*pow(T[i],4.0)) *
                			(b_soot[1] + b_soot[2]*T[i]+ b_soot[3]*pow(T[i],2.0) + b_soot[4]*pow(T[i],3.0) + b_soot[5]*pow(T[i],4.0));
        }
    }

// BOUNDARY FIELD CALCULATION
    forAll(mesh().boundary(), bid)
    {
       scalarField Tw = thermo_.T().boundaryField()[bid];
      
       scalarField& wsggmWeightingCoeffBC = tggCoeffBC.ref().boundaryFieldRef()[bid];
    
	    forAll(wsggmWeightingCoeffBC, i)
	    {
            const absorptionCoeffs::coeffArray& b_wc   = coeffs_[0][bandI].coeffs(Tw[i]);
            const absorptionCoeffs::coeffArray& b_soot = coeffs_[1][bandI].coeffs(Tw[i]);

	        if( (bandI == 0) | (bandI == 5) )
	        {
	        	wsggmWeightingCoeffBC[i] = (1.0 - (b_wc[1] + b_wc[2]*Tw[i]+ b_wc[3]*pow(Tw[i],2.0) + b_wc[4]*pow(Tw[i],3.0) + b_wc[5]*pow(Tw[i],4.0))) *
	                			(b_soot[1] + b_soot[2]*Tw[i]+ b_soot[3]*pow(Tw[i],2.0) + b_soot[4]*pow(Tw[i],3.0) + b_soot[5]*pow(Tw[i],4.0));
	        }
	        else
	        {
	        	wsggmWeightingCoeffBC[i] = (b_wc[1] + b_wc[2]*Tw[i]+ b_wc[3]*pow(Tw[i],2.0) + b_wc[4]*pow(Tw[i],3.0) + b_wc[5]*pow(Tw[i],4.0)) *
	                			(b_soot[1] + b_soot[2]*Tw[i]+ b_soot[3]*pow(Tw[i],2.0) + b_soot[4]*pow(Tw[i],3.0) + b_soot[5]*pow(Tw[i],4.0));
	        }

	    }
    }

    return tggCoeff + tggCoeffBC;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmCassolSoot::eCont(const label bandI) const
{
    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmCassolSoot::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    return E;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmCassolSoot::addIntensity
(
    const label i,
    const volScalarField& ILambda
) const
{

    return ILambda;
}


void Foam::radiation::wsggmCassolSoot::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda

) const
{
    a = dimensionedScalar("zero", dimless/dimLength, 0.0);

    for (label j=0; j<nBands_; j++)
    {
        Info<< "Calculating absorption in band: " << j << endl;
        aLambda[j].primitiveFieldRef() = this->a(j);

        a.primitiveFieldRef() =
            aLambda[j].primitiveField();

    }

}

void Foam::radiation::wsggmCassolSoot::correctNew // modified to include ggCoeffLambda -> ggCoeff (13/10/2014)
(
    volScalarField& ggCoeff,
    PtrList<volScalarField>& ggCoeffLambda

) const
{
    ggCoeff = dimensionedScalar("zero", dimless, 0.0);


    for (label j=0; j<nBands_; j++)
    {
        Info<< "Calculating weighting coefficient in band: " << j << endl;
        ggCoeffLambda[j] = this->ggCoeff(j);

        ggCoeff = ggCoeffLambda[j];

    }

}

// ************************************************************************* //
