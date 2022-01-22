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
    ANY WARuANTY; without even the implied waRuanty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "YaoSootModelTurbulent.H"
#include "betaPDF.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::YaoSootModelTurbulent<ThermoType>::YaoSootModelTurbulent
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),

    thermo(mesh.lookupObject<psiReactionThermo>("thermophysicalProperties")),

    Ysoot
    (
        IOobject
        (
            "Ysoot",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    fv
    (
        IOobject
        (
            "fv",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("soot", dimless, scalar(0.0))
    ),
    sootFormationRate
    (
        IOobject
        (
            "sootFormationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sootFormationRate", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
    ),    
    sootOxidationRate
    (
        IOobject
        (
            "sootOxidationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sootOxidationRate", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
    ),        
    sootConvection
    (
        IOobject
        (
            "sootConvection",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sootConvection", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
    ),      
    sootTimeDer
    (
        IOobject
        (
            "sootTimeDer",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sootTimeDer", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
    ),   
    diffusion
    (
        IOobject
        (
            "diffusion",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("diffusion", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
    ), 
    Z
    (
        IOobject
        (
            "Z",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Z", dimless, scalar(0.0))
    ),
    Zvar
    (
        IOobject
        (
            "Zvar",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_st
    (
        IOobject
        (
            "T_st",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T_st", dimensionSet(0,0,0,1,0,0,0), scalar(0.0))
    ),
    T_check
    (
        IOobject
        (
            "T_check",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T_check", dimensionSet(0,0,0,1,0,0,0), scalar(0.0))
    ),       
    rhoBar
    (
        IOobject
        (
            "rhoBar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoBar", dimensionSet(1,-3,0,0,0,0,0), scalar(0.0))
    ),           
    O2Concentration
    (
        IOobject
        (
            "O2Concentration",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("O2Concentration", dimensionSet(0,-3,0,0,1,0,0), scalar(0.0))
    ),


    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),

    solveSoot(coeffsDict_.lookup("solveSoot")),
    SGSFilter(coeffsDict_.lookup("SGSFilter")),
    	          
    rhoSoot
    (
        dimensionedScalar("rhoSoot", 
        dimensionSet(1,-3,0,0,0,0,0),
        readScalar(coeffsDict_.lookup("rhoSoot")))
    ),

    Af( readScalar(coeffsDict_.lookup("Af")) ),

    s (dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&> (thermo).s().value() ),
    O2Index (dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["O2"]),
    fuelIndex (dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()[thermo.lookup("fuel")]),

    YO2Inf(readScalar(thermo.lookup("YO2Inf"))),
    YFInf(readScalar(thermo.lookup("YFInf"))),

    Z_st( (YO2Inf/s)/(YFInf+YO2Inf/s) ),
    
    Z_sf( 2.5*Z_st ),
    Z_so( 1.25*Z_st ),
    
    gamma(2.25),
    Ta(2000.0),

    Asoot(160e3),
    Aox(120.0),
    EaOx(163540.0),

    MW_O2("MW_O2", dimensionSet(1,0,0,0,-1,0,0), thermo.composition().W(O2Index)*1e-3),    //(kg/mol)
    MW_Fuel("MW_Fuel", dimensionSet(1,0,0,0,-1,0,0), thermo.composition().W(fuelIndex)*1e-3), //(kg/mol)
    Ru("Ru", dimensionSet(1,2,-2,-1,-1,0,0), scalar(8.3145)),

    rho_ref( readScalar(coeffsDict_.lookup("rho_ref")) ),
    MW_ref( readScalar(coeffsDict_.lookup("MW_ref")) ),
    T_ref( readScalar(coeffsDict_.lookup("T_ref")) ),

    dX( readScalar(coeffsDict_.lookup("Xtilde_resolution")) ),  
    dXVar( readScalar(coeffsDict_.lookup("Xvar_resolution")) ),
    X_max( readScalar(coeffsDict_.lookup("Xtilde_max")) ),
    XVar_max( readScalar(coeffsDict_.lookup("Xvar_max")) ), 

    lookup_Tst()

{

    Info << "Z_st =  " << Z_st << endl;
    Info << "fuel molecular weight (kg/mol) = " << MW_Fuel.value() << endl;
    Info << "oxidizer molecular weight (kg/mol) = " << MW_O2.value() << endl;

    Info << "Soot formation/oxidation mix. frac. limits are "
         << Z_sf << " and " << Z_so << endl;

    Info << "Soot model Af = " << Af << endl;     


    //- Generating look-up tableLookup table of Beta-PDF integration functions
/*
    if (SGSFilter)
    {
        Info << "Generating lookup tables of Beta-PDF integrals" << endl;

        generateLookup("Tst");

        lookup_Tst.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Tst_data.dat"));
    }
*/
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::YaoSootModelTurbulent<ThermoType>::~YaoSootModelTurbulent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiation::YaoSootModelTurbulent<ThermoType>::correct()
{
    if (solveSoot)
    {  
        //access flow properties and species data
        const volScalarField& YO2 = thermo.composition().Y()[O2Index];
        const volScalarField& YFuel = thermo.composition().Y()[fuelIndex];
        const volScalarField& T = thermo.T();

        const volScalarField& rho = mesh().objectRegistry::template lookupObject<volScalarField>("rho");        

        const surfaceScalarField& phi = mesh().objectRegistry::template lookupObject<surfaceScalarField>("phi");

        const compressible::LESModel& lesModel = mesh().lookupObject<compressible::LESModel>
                                                (
                                                 turbulenceModel::propertiesName
                                                );
        
        //updating concentration of O2
        O2Concentration == YO2*rho/MW_O2;

        // Updating mixture fraction
        Z == (s*YFuel-YO2+YO2Inf)/(s*YFInf+YO2Inf);
        Z.max(0.0);
        Z.min(1.0);

        if (SGSFilter)
        {    
            // Solving transport equation of mixture fraction variance
            dimensionedScalar  k_small("k_small", dimensionSet(0,2,-2,0,0,0,0),1e-12);

            solve
            (
                    fvm::ddt(rho, Zvar)
                +   fvm::div(phi, Zvar)
                ==  fvm::laplacian(lesModel.alphaEff(), Zvar)
                +   2.0 * lesModel.alphat() * (fvc::grad(Z) & fvc::grad(Z))
                -   fvm::Sp(2.0 * rho * lesModel.epsilon()/max(lesModel.k(), k_small) , Zvar)
            );

            Zvar.max(0.0);
            Zvar.min(0.25);   
            Info << "Zvar min/max       = "   << min(Zvar).value() << " , " << max(Zvar).value() << endl;

            // Calculate formation and oxidation rates
            Info <<"updating soot formation/oxidation rates (Turbulent)" << endl;

            scalar intPDF_Ztriangle = 0.0;  //integration of f(Z)*P(Z)*dZ

            forAll(Ysoot, cellI)
            {
                intPDF_Ztriangle = integratePDF("Ztriangle", T_st[cellI], Z[cellI], Zvar[cellI]);

                T_st[cellI]    = (T[cellI] - T_ref+ T_ref*intPDF_Ztriangle) / max(intPDF_Ztriangle, 1e-6);

                T_check[cellI] = T_ref + T_st[cellI]*intPDF_Ztriangle - T_ref*intPDF_Ztriangle;

                rhoBar[cellI]  = 1.0 / max(integratePDF("density", T_st[cellI], Z[cellI], Zvar[cellI]), 1e-6);  
                                            
                sootFormationRate[cellI] = rhoBar[cellI] * integratePDF("formationRate", T_st[cellI], Z[cellI], Zvar[cellI]);

                sootOxidationRate[cellI] = rhoBar[cellI] * integratePDF("oxidationRate", T_st[cellI], Z[cellI], Zvar[cellI])
                                         * rho[cellI] * Ysoot[cellI] * Asoot;
            }

            Info << "T_st min/max       = "   << min(T_st).value() << " , " << max(T_st).value() << endl;
            Info << "T_check min/max    = "   << min(T_check).value() << " , " << max(T_check).value() << endl;

        }
        else
        {
            Info <<"updating soot formation/oxidation rates (Turbulent, with NO SGS model)" << endl;

            forAll(Ysoot, cellI)
            {
                sootFormationRate[cellI] = 0.0;
                sootOxidationRate[cellI] = 0.0;

                if ((Z[cellI] >= Z_so) && (Z[cellI] <= Z_sf))
                {

                    sootFormationRate[cellI] = Af * Foam::pow(rho[cellI], 2.0)
                                                * YFInf*(Z[cellI]-Z_st)/(1.0-Z_st)
                                                * Foam::pow(T[cellI], gamma)
                                                * Foam::exp(-Ta/T[cellI]);

                }

                if ((Z[cellI] >= 0.0) && (Z[cellI] <= Z_sf))
                {
                    sootOxidationRate[cellI] = rho[cellI] * Ysoot[cellI] * Asoot 
                                                * Aox 
                                                * O2Concentration[cellI]
                                                * Foam::pow(T[cellI], 0.5)
                                                * Foam::exp(-EaOx/Ru.value()/T[cellI]);

                }
            }
        }

        //Safety: limiting soot oxidation
        volScalarField sootOxidationLimiter = rho*Ysoot/mesh().time().deltaT()
                                             - fvc::div(phi, Ysoot)
                                             + fvc::laplacian(lesModel.alphat(), Ysoot)
                                             + sootFormationRate;

        forAll (sootOxidationRate, cellI)
        {
            sootOxidationRate[cellI] = max(0.0, min(sootOxidationRate[cellI], sootOxidationLimiter[cellI]));
        }         

        Info << "soot formation rate min/max = " << min(sootFormationRate).value() 
             << " , " << max(sootFormationRate).value() << endl;

        Info << "soot oxidation rate min/max = " << min(sootOxidationRate).value() 
             << " , " << max(sootOxidationRate).value() << endl;

        // Solve soot mass conservation equation
        fvScalarMatrix SootEqn
            (
                    fvm::ddt(rho, Ysoot)
                +   fvm::div(phi, Ysoot)
                ==  
                    fvm::laplacian(lesModel.alphat(), Ysoot)
                +   sootFormationRate
                -   sootOxidationRate
            );

        SootEqn.solve();      
        Ysoot.max(0.0);
        Ysoot.min(1.0);

        Info << "Ysoot min/max       = "   << min(Ysoot).value() << " , " << max(Ysoot).value() << endl;

        // Updating the density of the two-phase and soot vol. frac.
        fv == rho * Ysoot / rhoSoot; 

        Info << "soot vol fraction max = " << max(fv).value() << endl;

        //for diagnostic purposes only
        sootTimeDer     = fvc::ddt(rho, Ysoot);
        sootConvection  = fvc::div(phi, Ysoot);
        diffusion       = fvc::laplacian(lesModel.alphat(), Ysoot);        
   }

}

template<class ThermoType>
double Foam::radiation::YaoSootModelTurbulent<ThermoType>::sourceFunc(
                                                const word& sourceName,
                                                const double& T_st,
                                                const double& Z
                                                )
{

    double Ztriangle    = 0.0;
    double rho_approx   = 0.0;
    double T_approx     = 0.0;
    double MW_approx    = 0.0;
    double YO2_approx   = 0.0;

    //approximate temperature
    if (Z >=0 && Z <= Z_st)
    {
        Ztriangle = Z/Z_st;
    }
    else
    {
        Ztriangle = (1.0-Z)/(1-Z_st);
    } 
    T_approx = T_ref + (T_st-T_ref) * Ztriangle;

    //approximate YO2
    double Z_L = 0.875 * Z_st;
    double B = Z_st-Z_L;
    double A = YO2Inf * (1-Z_L/Z_st)*Foam::exp(Z_L/B);

    if (Z >=0 && Z <= Z_L)
    {
        YO2_approx = YO2Inf * (1.0 - Z/Z_st);
    }
    else
    {
        YO2_approx = A*Foam::exp(-Z/B);        
    } 

    //approximate density
    MW_approx   = 1.0/ ( (1.0-Z)/MW_ref + Z/MW_Fuel.value() );
    rho_approx  = rho_ref * MW_approx/MW_ref * T_ref/T_approx;

    // return sootFormationRate/rho
    if (sourceName == "formationRate")
    {
        if (Z>=Z_so && Z<=Z_sf)
        {
             return Af * Foam::pow(rho_approx, 2.0)
                    * YFInf*(Z-Z_st)/(1.0-Z_st)
                    * Foam::pow(T_approx, gamma)
                    * Foam::exp(-Ta/T_approx)
                    / rho_approx;           
        }
        else
        {
            return 0;
        }        
    }

    // return sootOxidationRate/rho / Ysoot
    else if (sourceName == "oxidationRate")
    {
        if (Z>=0.0 && Z<=Z_sf)
        {
             return Aox 
                    * rho_approx * YO2_approx/MW_O2.value()
                    * Foam::pow(T_approx, 0.5)
                    * Foam::exp(-EaOx/Ru.value()/T_approx)
                    / rho_approx;           
        }
        else
        {
            return 0;
        }        
    } 

    //return 1/rho
    else if (sourceName == "density")
    {
        return 1.0/rho_approx;
    }

    // return Z triangle profile
    else if (sourceName == "Ztriangle")
    {
        return Ztriangle;

    }

    else
    {
        FatalErrorInFunction
        << "Attempt to use a function in PDF integration that is not listed"
        << abort(FatalError);
        return 0;        
    }        
                    
}

/*
template<class ThermoType>
void Foam::radiation::YaoSootModelTurbulent<ThermoType>::generateLookup(
                                                const word& sourceName
                                                )
{
        int m = round(X_max/dX) + 1 ;
        int n = round(XVar_max/dXVar) + 1;

        List<scalar> X, XVar;
        X.resize(m,0.0);
        XVar.resize(n,0.0);

        X[0]    = 0.0;
        XVar[0] = 0.0;
        for(int i=1 ; i<m ; i++)
        {
            X[i] = X[i-1] + dX;
        }
        for(int i=1 ; i<n ; i++)
        {
            XVar[i] = XVar[i-1] + dXVar;
        }

        scalar Func;

        List< Tuple2<scalar,scalar> > Func_column;

        List< Tuple2<scalar, List< Tuple2<scalar,scalar> > > > Func_data; 


        // creating lookup table of Func(sourceName)    
        Func_column.setSize(n); 
        Func_data.setSize(m);

        ofstream Func_CSV(mesh().time().path()/"constant"/sourceName+"_data.csv");
        OFstream Func_os(mesh().time().path()/"constant"/sourceName+"_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {
                Func = integratePDF(sourceName, X[i] , XVar[j]);
                Func_column[j] = Tuple2<scalar,scalar>(XVar[j] , Func);

                Func_CSV << Func << ",";
            }
            Func_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Func_column);
            Func_CSV << "\n";
        }
        Func_os << Func_data << endl;

        Func_column.resize(0); 
        Func_data.resize(0);

        Info << "   Done writing lookup table of " << sourceName << endl;        

}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
