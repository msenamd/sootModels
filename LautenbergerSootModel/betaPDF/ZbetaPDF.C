
#include "ZbetaPDF.H"


//integratio of: F_sf(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_sf(
			double Ztilde,
			double Zvar,
			sootSource SS
			)
{
	if(Zvar < 1e-6)
	{
		return SS.F_sf(Ztilde);
	}

	//beta function parameters
	double a = Ztilde * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);
	double b = (1.0-Ztilde) * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);

	double Zmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>1500)
	{
		a = 1500;
		b = (a-1-Zmax*(a-2))/Zmax;
	}
	if(b>1500)
	{
		b = 1500;
		a = (1+Zmax*(b-2))/(1-Zmax);
	}

	if(a>0 && b>0)
	{
		// integration variables
		std::vector<double> v;
		v.assign(12,0.0);

		v[0]=epsilon;
		v[1]=1e-5;
		v[2]=1e-4;
		v[3]=1e-3;
		v[4]=0.05;
		v[5]=0.15;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double Z1,Z2;
		std::vector<double> ZZ , fBeta_num , fBeta_denom;
		int N=20;
		int M=200;

		//Numerator
		double n1,n2,n3;
		n1 = SS.F_sf(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.F_sf(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between Z1 & Z2
		// ZZ has values between Z1 & Z2
		// fBeta_num is Y(Z)*Beta(Z)
		// fBeta_denom is Beta(Z)		
		for(int j=0 ; j<11 ; j++)
		{
			Z1=v[j];
			Z2=v[j+1];

			if((Z1==0.05) | (Z1==0.15))
			{
				ZZ.assign(M+1,0.0);// use M
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
				fBeta_num[i]   	= SS.F_sf(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);
				
				for(int i=1 ; i<M+1 ; i++)
				{
					ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
					fBeta_num[i]   	= SS.F_sf(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				ZZ.assign(N+1,0.0);// use N
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);				


				int i=0;
				ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
				fBeta_num[i]  	= SS.F_sf(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				for(int i=1 ; i<N+1 ; i++)
				{
					ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
					fBeta_num[i]  	= SS.F_sf(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		if ((n1+n2+n3) == 0.0 && (d1+d2+d3) ==0.0)
		{
			return 0;
		}
		else
		{
			return (n1+n2+n3)/(d1+d2+d3);
		}
		
	}
	else
	{
		return 0;
	}
}
  

//integratio of: F_so(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_so(
			double Ztilde,
			double Zvar,
			sootSource SS
			)
{

	if(Zvar < 1e-6)
	{
		return SS.F_so(Ztilde);
	}

	//beta function parameters
	double a = Ztilde * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);
	double b = (1.0-Ztilde) * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);

	double Zmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>2000)
	{
		a = 2000;
		b = (a-1-Zmax*(a-2))/Zmax;
	}
	if(b>2000)
	{
		b = 2000;
		a = (1+Zmax*(b-2))/(1-Zmax);
	}
	
	if(a>0 && b>0)
	{
		// integration variables
		std::vector<double> v;
		v.assign(12,0.0);

		v[0]=epsilon;
		v[1]=1e-5;
		v[2]=1e-4;
		v[3]=1e-3;
		v[4]=0.02;
		v[5]=0.07;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double Z1,Z2;
		std::vector<double> ZZ , fBeta_num , fBeta_denom;
		int N=20;
		int M=200;

		//Numerator
		double n1,n2,n3;
		n1 = SS.F_so(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.F_so(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between Z1 & Z2
		// ZZ has values between Z1 & Z2
		// fBeta_num is Y(Z)*Beta(Z)
		// fBeta_denom is Beta(Z)		
		for(int j=0 ; j<11 ; j++)
		{
			Z1=v[j];
			Z2=v[j+1];

			if((Z1==0.02) | (Z1==0.07))
			{
				ZZ.assign(M+1,0.0);// use M
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
				fBeta_num[i]   	= SS.F_so(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);
				for(int i=1 ; i<M+1 ; i++)
				{
					ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
					fBeta_num[i]   	= SS.F_so(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				ZZ.assign(N+1,0.0);// use N
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);

				int i=0;
				ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
				fBeta_num[i]  	= SS.F_so(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				for(int i=1 ; i<N+1 ; i++)
				{
					ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
					fBeta_num[i]  	= SS.F_so(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		if ((n1+n2+n3) == 0.0 && (d1+d2+d3) ==0.0)
		{
			return 0;
		}
		else
		{
			return (n1+n2+n3)/(d1+d2+d3);
		}
	}
	else
	{
		return 0;
	}
}



//integratio of: F_rho(Z) Ptilde(Z) dZ
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FZbeta_rho(
			double Ztilde,
			double Zvar,
			sootSource SS
			)
{

	if(Zvar < 1e-6)
	{
		return SS.F_rho(Ztilde);
	}
	
	//beta function parameters
	double a = Ztilde * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);
	double b = (1.0-Ztilde) * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);

	double Zmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>500)
	{
		a = 500;
		b = (a-1-Zmax*(a-2))/Zmax;
	}
	if(b>500)
	{
		b = 500;
		a = (1+Zmax*(b-2))/(1-Zmax);
	}	

	if(a>0 && b>0)
	{
		// integration variables
		std::vector<double> v;
		v.assign(12,0.0);

		v[0]=epsilon;
		v[1]=1e-5;
		v[2]=1e-4;
		v[3]=1e-3;
		v[4]=0.02;
		v[5]=0.07;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double Z1,Z2;
		std::vector<double> ZZ , fBeta_num , fBeta_denom;
		int N=20;
		int M=200;

		//Numerator
		double n1,n2,n3;
		n1 = SS.F_rho(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.F_rho(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between Z1 & Z2
		// ZZ has values between Z1 & Z2
		// fBeta_num is Y(Z)*Beta(Z)
		// fBeta_denom is Beta(Z)		
		for(int j=0 ; j<11 ; j++)
		{
			Z1=v[j];
			Z2=v[j+1];

			if((Z1==0.02) | (Z1==0.07))
			{
				ZZ.assign(M+1,0.0);// use M
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
				fBeta_num[i]   	= SS.F_rho(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);
				for(int i=1 ; i<M+1 ; i++)
				{
					ZZ[i]   		= Z1 + i* (Z2-Z1) / M; 
					fBeta_num[i]   	= SS.F_rho(ZZ[i]) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				ZZ.assign(N+1,0.0);// use N
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);

				int i=0;
				ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
				fBeta_num[i]  	= SS.F_rho(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
				for(int i=1 ; i<N+1 ; i++)
				{
					ZZ[i]    	  	= Z1 + i* (Z2-Z1) / N; 
					fBeta_num[i]  	= SS.F_rho(ZZ[i]) * std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);
					fBeta_denom[i] 	= std::pow(ZZ[i],a-1) * std::pow(1-ZZ[i],b-1);

					//integration step
					n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		return (n1+n2+n3)/(d1+d2+d3);

	}
	else
	{
		return SS.F_rho(Ztilde);
	}
}


// ************************************************************************* //
