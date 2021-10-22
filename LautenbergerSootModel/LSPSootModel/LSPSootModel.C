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

\*---------------------------------------------------------------------------*/

#include "LSPSootModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::LSPSootModel<ThermoType>::LSPSootModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),
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
    rho2
    (
        IOobject
        (
            "rho2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho2", dimMass/dimVol, scalar(0.0))
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
    oxidationLimiter
    (
        IOobject
        (
            "oxidationLimiter",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("oxidationLimiter", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
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
    thermophoresis
    (
        IOobject
        (
            "thermophoresis",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("thermophoresis", dimensionSet(1,-3,-1,0,0,0,0), scalar(0.0))
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
    Zvar_SGS
    (
        IOobject
        (
            "Zvar_SGS",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Tstar
    (
        IOobject
        (
            "Tstar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tstar", dimless, scalar(0.0))
    ),
    TstarVar
    (
        IOobject
        (
            "TstarVar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TstarVar", dimless, scalar(0.0))
    ),
    rhobar
    (
        IOobject
        (
            "rhobar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhobar", dimensionSet(1,-3,0,0,0,0,0), scalar(0.0))
    ),

    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    
    solveSoot_(coeffsDict_.lookup("solveSoot")),

	turbulence_(coeffsDict_.lookup("turbulence")),

    oxidation_(coeffsDict_.lookup("oxidation")),
    	
	rho_soot
    (
        dimensionedScalar("rho_soot", 
        dimensionSet(1,-3,0,0,0,0,0),
        readScalar(coeffsDict_.lookup("rho_soot")))
    ),
        
    smokePointHeight(readScalar(coeffsDict_.lookup("smokePointHeight"))),
        
    T_adiabatic
    (
        dimensionedScalar("T_adiabatic", 
        dimensionSet(0,0,0,1,0,0,0),
        readScalar(coeffsDict_.lookup("adiabaticFlameTemperature")))
    ),
    
    T_inf
    (
        dimensionedScalar("T_inf", 
        dimensionSet(0,0,0,1,0,0,0),
        readScalar(coeffsDict_.lookup("T_inf")))
    ),
    
    T_oxidizer
    (
        dimensionedScalar("T_oxidizer", 
        dimensionSet(0,0,0,1,0,0,0),
        readScalar(coeffsDict_.lookup("T_oxidizer")))
    ),
    
    rho_oxidizer(readScalar(coeffsDict_.lookup("rho_oxidizer"))),
    
    MW_oxidizer(readScalar(coeffsDict_.lookup("MW_oxidizer"))),
   	    
    proportionalityConst(readScalar(coeffsDict_.lookup("proportionalityConst"))),
    
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    
    thermo(mesh.lookupObject<psiReactionThermo>("thermophysicalProperties")),
    
    fuelIndex (dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()[thermo.lookup("fuel")]),
    
    O2Index (dynamic_cast<const reactingMixture<gasHThermoPhysics>&> (thermo).species()["O2"]),
    
    MW_fuel (dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&> (thermo).speciesData()[fuelIndex].W() ),
    
    s (dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&> (thermo).s().value() ),
    
    YO2Inf(readScalar(thermo.lookup("YO2Inf"))),

    YFInf(readScalar(thermo.lookup("YFInf"))),

    Z_st ( (YO2Inf/s)/(YFInf+YO2Inf/s) ),

    peakFormationRate (1.1 * (0.106/smokePointHeight) * (28.0/MW_fuel) * YFInf),

    SS(    			
		    	Z_st,
		    	peakFormationRate,
		    	MW_fuel,
		    	MW_oxidizer,
		    	T_adiabatic.value(),
		    	T_inf.value(),
		    	T_oxidizer.value(),
		    	rho_oxidizer
    	),

    lookup_Fsf(),
    lookup_Gsf(),
    lookup_Fso(),
    lookup_Gso(),
    lookup_Frho(),
    lookup_Grho()
{

    Info << "fuel molecular weight = " << MW_fuel << endl;
    Info << "oxidizer molecular weight = " << MW_oxidizer << endl;
	Info << "stoichiometric mixture fraction: Z_st = " << Z_st << endl;

    Info << "peak soot formation rate = " << peakFormationRate << endl;

    //-calculating soot polynomial coefficients
	Info << "Calculating SPH model polynomial coefficients" << endl;

    SS.calcCoeff();


	//- Generating look-up tableLookup table of Beta-PDF integration functions
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(turbulence_)
    {
        Info << "Generating lookup tables of Beta-PDF integrals" << endl;

	    scalar dX 		= readScalar(coeffsDict_.lookup("Xtilde_resolution")); 	//lookup table resolution
	    scalar dXVar 	= readScalar(coeffsDict_.lookup("XVariance_resolution"));

	    scalar X_max 	= readScalar(coeffsDict_.lookup("Xtilde_max"));		//upper bounds
	    scalar XVar_max = readScalar(coeffsDict_.lookup("Xvariance_max")); 

	    int m = round(X_max/dX) + 1 ;
	    int n = round(XVar_max/dXVar) + 1;

        List<scalar> X, XVar;
        X.resize(m,0.0);
        XVar.resize(n,0.0);

        X[0] 	= 0.0;
        XVar[0] = 0.0;
        for(int i=1 ; i<m ; i++)
        {
            X[i] = X[i-1] + dX;
        }
        for(int i=1 ; i<n ; i++)
        {
            XVar[i] = XVar[i-1] + dXVar;
        }

        scalar Fsf , Gsf, Fso, Gso, Frho, Grho;

        List< Tuple2<scalar,scalar> > Fsf_column, Gsf_column, Fso_column, Gso_column,
        							  Frho_column, Grho_column;

        List< Tuple2<scalar, List< Tuple2<scalar,scalar> > > > Fsf_data, Gsf_data, Fso_data, Gso_data,
        													   Frho_data, Grho_data;      


    // creating lookup table of F_sf    
        Fsf_column.setSize(n); 
        Fsf_data.setSize(m);

        ofstream Fsf_CSV(mesh.time().path()/"constant"/"Fsf_data.csv");
        OFstream Fsf_os(mesh.time().path()/"constant"/"Fsf_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {
                Fsf = FZbeta_sf(X[i] , XVar[j] , SS);
                Fsf_column[j] = Tuple2<scalar,scalar>(XVar[j] , Fsf);

                Fsf_CSV << Fsf << ",";
            }
            Fsf_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Fsf_column);
            Fsf_CSV << "\n";
        }
        Fsf_os << Fsf_data << endl;
        lookup_Fsf.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Fsf_data.dat"));
        Fsf_column.resize(0); 
        Fsf_data.resize(0);
        Info << "	Done Generating of F_sf" << endl;


	// creating lookup table of G_sf 
        Gsf_column.setSize(n);
        Gsf_data.setSize(m);

        ofstream Gsf_CSV(mesh.time().path()/"constant"/"Gsf_data.csv");
        OFstream Gsf_os(mesh.time().path()/"constant"/"Gsf_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {
                Gsf = GTbeta_sf(X[i] , XVar[j] , SS);
                Gsf_column[j] = Tuple2<scalar,scalar>(XVar[j] , Gsf);

                Gsf_CSV << Gsf<< ",";
            }
            Gsf_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Gsf_column);
            Gsf_CSV << "\n";
        }
        Gsf_os << Gsf_data << endl;
        lookup_Gsf.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Gsf_data.dat"));
        Gsf_column.resize(0);
        Gsf_data.resize(0);
        Info << "	Done Generating of G_sf" << endl;


	// creating lookup table of F_so
        Fso_column.setSize(n);
        Fso_data.setSize(m);

        ofstream Fso_CSV(mesh.time().path()/"constant"/"Fso_data.csv");
        OFstream Fso_os(mesh.time().path()/"constant"/"Fso_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {
                Fso = FZbeta_so(X[i] , XVar[j] , SS);
                Fso_column[j] = Tuple2<scalar,scalar>(XVar[j] , Fso);

                Fso_CSV << Fso << ",";
            }
            Fso_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Fso_column);
            Fso_CSV << "\n";
        }
        Fso_os << Fso_data << endl;
        lookup_Fso.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Fso_data.dat"));
        Fso_column.resize(0);
        Fso_data.resize(0);        
        Info << "	Done Generating of F_so" << endl;        


	// creating lookup table of G_so 
        Gso_column.setSize(n);
        Gso_data.setSize(m);

        ofstream Gso_CSV(mesh.time().path()/"constant"/"Gso_data.csv");
        OFstream Gso_os(mesh.time().path()/"constant"/"Gso_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {

                Gso = GTbeta_so(X[i] , XVar[j] , SS);
                Gso_column[j] = Tuple2<scalar,scalar>(XVar[j] , Gso);

                Gso_CSV << Gso << ",";
            }
            Gso_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Gso_column);
            Gso_CSV << "\n";
        }
        Gso_os << Gso_data << endl;
        lookup_Gso.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Gso_data.dat"));
        Gso_column.resize(0);
        Gso_data.resize(0);
        Info << "	Done Generating of G_so" << endl;


	// creating lookup table of F_rho
        Frho_column.setSize(n);
        Frho_data.setSize(m);

        ofstream Frho_CSV(mesh.time().path()/"constant"/"Frho_data.csv");
        OFstream Frho_os(mesh.time().path()/"constant"/"Frho_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {
                Frho = FZbeta_rho(X[i] , XVar[j] , SS);
                Frho_column[j] = Tuple2<scalar,scalar>(XVar[j] , Frho);

                Frho_CSV << Frho << ",";
            }
            Frho_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Frho_column);
            Frho_CSV << "\n";
        }
        Frho_os << Frho_data << endl;
        lookup_Frho.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Frho_data.dat"));
        Frho_column.resize(0);
        Frho_data.resize(0);        
        Info << "	Done Generating of F_rho" << endl;        


	// creating lookup table of G_rho
        Grho_column.setSize(n);
        Grho_data.setSize(m);

        ofstream Grho_CSV(mesh.time().path()/"constant"/"Grho_data.csv");
        OFstream Grho_os(mesh.time().path()/"constant"/"Grho_data.dat");

        for(int i=0 ; i<m ; i++)
        {
            for(int j=0 ; j<n ; j++)
            {

                Grho = GTbeta_rho(X[i] , XVar[j] , SS);
                Grho_column[j] = Tuple2<scalar,scalar>(XVar[j] , Grho);

                Grho_CSV << Grho << ",";
            }
            Grho_data[i] = Tuple2<scalar, List< Tuple2<scalar,scalar>>> (X[i] , Grho_column);
            Grho_CSV << "\n";
        }
        Grho_os << Grho_data << endl;
        lookup_Grho.reset(new interpolation2DTable<scalar> (mesh.time().path()/"constant"/"Grho_data.dat"));
        Grho_column.resize(0);
        Grho_data.resize(0);
        Info << "	Done Generating of G_rho" << endl;

    }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::LSPSootModel<ThermoType>::~LSPSootModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiation::LSPSootModel<ThermoType>::correct()
{

    //access flow properties and species data
    const volScalarField& rho = mesh().objectRegistry::template lookupObject<volScalarField>("rho");
    const volScalarField& T = thermo.T();
    const surfaceScalarField& phi = mesh().objectRegistry::template lookupObject<surfaceScalarField>("phi");
    const volVectorField& U = mesh().objectRegistry::template lookupObject<volVectorField>("U");

    // Updating mixture fraction
    Z = mesh().objectRegistry::template lookupObject<volScalarField>("ft");

    // Calculating normalized temperature
    Tstar = (T - T_inf) / (T_adiabatic-T_inf); 
    Tstar.max(0.0);
    Tstar.min(1.0);

    // Updating the density of the two-phase soot+gas
    rho2 = rho/(1.0-Ysoot);

    surfaceScalarField phi2("phi2", linearInterpolate(rho2*U) & mesh().Sf());

    if (solveSoot_)
    {    

        const turbulenceModel& turbulence 
            = mesh().objectRegistry::template lookupObject<turbulenceModel>("turbulenceProperties");
     
        const scalar Prt = 0.5; 

        //update soot oxidation source limiter
        volScalarField Ysoot_01 = Ysoot;
        Ysoot_01.min(1.0);
        Ysoot_01.max(0.0);
        oxidationLimiter = rho2*Ysoot_01/mesh().time().deltaT()
                         - fvc::div(phi2, Ysoot_01) 
                         + fvc::laplacian(rho2*turbulence.nut()/Prt, Ysoot_01)
                         + fvc::div(0.54*thermo.mu()/T * fvc::grad(T) * Ysoot_01)
                         + sootFormationRate;


        if(turbulence_)
        {
       		dimensionedScalar  k_small("k_small", dimensionSet(0,2,-2,0,0,0,0),1e-12);

	        //Transport equation of variance of mixture fraction
	        solve
	            (
	                	fvm::ddt(rho, Zvar_SGS)
	                +	fvm::div(phi, Zvar_SGS)
	                == 	fvm::laplacian( thermo.alpha() + rho*turbulence.nut()/Prt , Zvar_SGS)
	                + 	2.0 * rho*turbulence.nut()/Prt * (fvc::grad(Z) & fvc::grad(Z))
	                - 	fvm::Sp(2.0 * rho * turbulence.epsilon()/max(turbulence.k(),k_small) , Zvar_SGS)
	            );

	        Zvar_SGS.max(0.0);
	        Zvar_SGS.min(0.25);   

	        Info << "Zvar min       = "   << min(Zvar_SGS).value() << endl;
	        Info << "Zvar max       = "   << max(Zvar_SGS).value() << endl;

	        //calculating normalized temperature variance
	    	TstarVar = pow(proportionalityConst * Tstar/max(Z,1e-4) , 2.0) * Zvar_SGS;
	        TstarVar.max(0.0);
	        TstarVar.min(0.25);

	        Info << "TstarVar min    = "   << min(TstarVar).value() << endl;
	        Info << "TstarVar max    = "   << max(TstarVar).value() << endl;
        	
	        // Updating soot source term from PDF integration

	        Info <<"calculating soot source terms from Beta-PDF interpolation" << endl;

	        forAll(Ysoot, cellI)
	        {
	        	rhobar[cellI] = 1.0 / 
	        					(
	        						lookup_Frho()(Z[cellI] , Zvar_SGS[cellI]) *
	        						lookup_Grho()(Tstar[cellI] , TstarVar[cellI]) 
	        					);

	            sootFormationRate[cellI] = rhobar[cellI] * 
			                                (
			                                    lookup_Fsf()(Z[cellI] , Zvar_SGS[cellI]) 
			                                    * 
			                                    lookup_Gsf()(Tstar[cellI] , TstarVar[cellI])
			                                 );         		
	 
                sootOxidationRate[cellI] = rhobar[cellI] * 
                                            (
                                                lookup_Fso()(Z[cellI] , Zvar_SGS[cellI]) 
                                                * 
                                                lookup_Gso()(Tstar[cellI] , TstarVar[cellI])
                                            );                    

                sootOxidationRate[cellI] = max(0.0, min(sootOxidationRate[cellI], oxidationLimiter[cellI]));
	        } 
    	}
    	else
    	{
	        // Updating soot source term from lamianr polynomials

	        Info <<"calculating soot source terms (laminar)" << endl;

	        forAll(Ysoot, cellI)
	        {
	        	rhobar[cellI] = 1.0 / 
	        					(
	        						SS.F_rho(Z[cellI]) *
	        						SS.G_rho(Tstar[cellI]) 
	        					);

	            sootFormationRate[cellI] = rhobar[cellI] * 
			                                (
			                                    SS.F_sf(Z[cellI]) 
			                                    * 
			                                    SS.G_sf(Tstar[cellI])
			                                 );         		

                sootOxidationRate[cellI] = rhobar[cellI] * 
                                            (
                                                SS.F_so(Z[cellI]) 
                                                * 
                                                SS.G_so(Tstar[cellI])
                                            );  
                
                sootOxidationRate[cellI] = max(0.0, min(sootOxidationRate[cellI], oxidationLimiter[cellI]));
	        }
    	}

        if(!oxidation_)
        {
            sootOxidationRate *= scalar(0.0);
        }

        Info << "soot formation rate max = " << max(sootFormationRate).value() << endl;
        Info << "soot oxidation rate max = " << max(sootOxidationRate).value() << endl;

        
        // Solve soot mass conservation equation
        fvScalarMatrix SootEqn
            (
                    fvm::ddt(rho2, Ysoot)
                +   fvm::div(phi2, Ysoot)
                ==  
                    fvm::laplacian(rho2*turbulence.nut()/Prt, Ysoot)
                +   fvc::div(0.54*thermo.mu()/T * fvc::grad(T) * Ysoot)
                +   sootFormationRate
                -	sootOxidationRate
            );

        SootEqn.solve(); 

        Info << "soot mass fraction min = " << min(Ysoot).value() << endl;
        Info << "soot mass fraction max = " << max(Ysoot).value() << endl;

        //for diagnostic purposes only
        sootTimeDer     = fvc::ddt(rho2, Ysoot);
        sootConvection  = fvc::div(phi2, Ysoot);
        thermophoresis  = fvc::div(0.54*thermo.mu()/T * fvc::grad(T) * Ysoot);

        fv = rho2 * Ysoot / rho_soot; 
        fv.max(0.0);

        Info << "soot vol fraction max = " << max(fv).value() << endl;

    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
