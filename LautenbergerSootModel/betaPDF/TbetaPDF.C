
#include "TbetaPDF.H"


//integratio of: G_sf(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_sf(
			double TstarTilde,
			double TstarVar,
			sootSource SS
			)
{

	if(TstarVar < 1e-6)
	{
		return SS.G_sf(TstarTilde);
	}

	//beta function parameters
	double a = TstarTilde * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);
	double b = (1.0-TstarTilde) * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);

	double Tmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>500)
	{
		a = 500;
		b = (a-1-Tmax*(a-2))/Tmax;
	}
	if(b>500)
	{
		b = 500;
		a = (1+Tmax*(b-2))/(1-Tmax);
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
		v[4]=1e-2;
		v[5]=0.1;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double T1,T2;
		std::vector<double> TT , fBeta_num , fBeta_denom;
		int N=20;
		int M=50;

		//Numerator
		double n1,n2,n3;
		n1 = SS.G_sf(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.G_sf(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between T1 & T2
		// TT has values between T1 & T2
		// fBeta_num is Y(T)*Beta(T)
		// fBeta_debin is Beta(T)
		for(int j=0 ; j<11 ; j++)
		{
			T1=v[j];
			T2=v[j+1];

			if(T1==0.1)
			{
				TT.assign(M+1,0.0);// use M=50
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				TT[i]   	   = T1 + i* (T2-T1) / M; 
				fBeta_num[i]   = SS.G_sf(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<M+1 ; i++)
				{
					TT[i]   	   = T1 + i* (T2-T1) / M; 
					fBeta_num[i]   = SS.G_sf(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				TT.assign(N+1,0.0);// use N=20
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);

				int i=0;
				TT[i]          = T1 + i* (T2-T1) / N; 
				fBeta_num[i]   = SS.G_sf(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<N+1 ; i++)
				{
					TT[i]          = T1 + i* (T2-T1) / N; 
					fBeta_num[i]   = SS.G_sf(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		return (n1+n2+n3)/(d1+d2+d3);
	}
	else
	{
		return 0;
	}

}
  

//integratio of: G_so(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_so(
			double TstarTilde,
			double TstarVar,
			sootSource SS
			)
{

	if(TstarVar < 1e-6)
	{
		return SS.G_so(TstarTilde);
	}

	//beta function parameters
	double a = TstarTilde * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);
	double b = (1.0-TstarTilde) * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);
	
	double Tmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>500)
	{
		a = 500;
		b = (a-1-Tmax*(a-2))/Tmax;
	}
	if(b>500)
	{
		b = 500;
		a = (1+Tmax*(b-2))/(1-Tmax);
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
		v[4]=1e-2;
		v[5]=0.1;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double T1,T2;
		std::vector<double> TT , fBeta_num , fBeta_denom;
		int N=20;
		int M=50;

		//Numerator
		double n1,n2,n3;
		n1 = SS.G_so(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.G_so(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between T1 & T2
		// TT has values between T1 & T2
		// fBeta_num is Y(T)*Beta(T)
		// fBeta_debin is Beta(T)
		for(int j=0 ; j<11 ; j++)
		{
			T1=v[j];
			T2=v[j+1];

			if(T1==0.1)
			{
				TT.assign(M+1,0.0);// use M=50
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				TT[i]   	   = T1 + i* (T2-T1) / M; 
				fBeta_num[i]   = SS.G_so(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<M+1 ; i++)
				{
					TT[i]   	   = T1 + i* (T2-T1) / M; 
					fBeta_num[i]   = SS.G_so(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				TT.assign(N+1,0.0);// use N=20
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);

				int i=0;
				TT[i]          = T1 + i* (T2-T1) / N; 
				fBeta_num[i]   = SS.G_so(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<N+1 ; i++)
				{
					TT[i]          = T1 + i* (T2-T1) / N; 
					fBeta_num[i]   = SS.G_so(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		return (n1+n2+n3)/(d1+d2+d3);
	}
	else
	{
		return 0;
	}

}



//integratio of: G_rho(T*) Ptilde(T*) dT*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double GTbeta_rho(
			double TstarTilde,
			double TstarVar,
			sootSource SS
			)
{

	if(TstarVar < 1e-6)
	{
		return SS.G_rho(TstarTilde);
	}

	//beta function parameters
	double a = TstarTilde * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);
	double b = (1.0-TstarTilde) * (TstarTilde*(1.0-TstarTilde)/max(TstarVar,1e-6) - 1.0);
	
	double Tmax = 1.0/(1.0+(b-1.0)/(a-1.0));
	double epsilon = 1E-6;

	//capping of a&b
	if(a>500)
	{
		a = 500;
		b = (a-1-Tmax*(a-2))/Tmax;
	}
	if(b>500)
	{
		b = 500;
		a = (1+Tmax*(b-2))/(1-Tmax);
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
		v[4]=1e-2;
		v[5]=0.1;
		v[6]=0.9;
		v[7]=1.0-1e-2;
		v[8]=1.0-1e-3;
		v[9]=1.0-1e-4;
		v[10]=1.0-1e-5;
		v[11]=1.0-epsilon;

		double T1,T2;
		std::vector<double> TT , fBeta_num , fBeta_denom;
		int N=20;
		int M=50;

		//Numerator
		double n1,n2,n3;
		n1 = SS.G_rho(0.0) * std::pow(epsilon,a)/a;
		n2 = 0.0;
		n3 = SS.G_rho(1.0) * std::pow(epsilon,b)/b;

		//Denominator
		double d1,d2,d3;
		d1 = std::pow(epsilon,a)/a;
		d2 = 0.0;
		d3 = std::pow(epsilon,b)/b;

		// setting up the bounded integral data between T1 & T2
		// TT has values between T1 & T2
		// fBeta_num is Y(T)*Beta(T)
		// fBeta_debin is Beta(T)
		for(int j=0 ; j<11 ; j++)
		{
			T1=v[j];
			T2=v[j+1];

			if(T1==0.1)
			{
				TT.assign(M+1,0.0);// use M=50
				fBeta_num.assign(M+1,0.0);
				fBeta_denom.assign(M+1,0.0);

				int i=0;
				TT[i]   	   = T1 + i* (T2-T1) / M; 
				fBeta_num[i]   = SS.G_rho(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<M+1 ; i++)
				{
					TT[i]   	   = T1 + i* (T2-T1) / M; 
					fBeta_num[i]   = SS.G_rho(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
			else
			{
				TT.assign(N+1,0.0);// use N=20
				fBeta_num.assign(N+1,0.0);
				fBeta_denom.assign(N+1,0.0);

				int i=0;
				TT[i]          = T1 + i* (T2-T1) / N; 
				fBeta_num[i]   = SS.G_rho(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
				for(int i=1 ; i<N+1 ; i++)
				{
					TT[i]          = T1 + i* (T2-T1) / N; 
					fBeta_num[i]   = SS.G_rho(TT[i]) * std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);
					fBeta_denom[i] = std::pow(TT[i],a-1.0) * std::pow(1.0-TT[i],b-1.0);

					//integration step
					n2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
					d2 += 0.5 * (TT[i]-TT[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
				}
			}
		}

		return (n1+n2+n3)/(d1+d2+d3);
	}
	else
	{
		return SS.G_rho(TstarTilde);
	}

}

// ************************************************************************* //
