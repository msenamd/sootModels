#include <vector>

template<class ThermoType>
double Foam::radiation::YaoSootModelTurbulent<ThermoType>::integratePDF(
                const word& sourceName,
                const double& T_st,
                const double& etaTilde,
                const double& etaVar               
			)
{

double epsilon = 1e-6;

if(etaVar > epsilon)
{
	//beta function parameters
	double a = etaTilde * (etaTilde*(1.0-etaTilde)/etaVar - 1.0);
	double b = (1.0-etaTilde) * (etaTilde*(1.0-etaTilde)/etaVar - 1.0);

	//safety statements to prevent singularity or overflow
	double 	etaMax, anew, bnew;

	if(a>500)
	{
		etaMax = 1.0 / (1.0 + (b-1.0)/(a-1.0));
		anew = 500;
		bnew = (anew-1-etaMax*(anew-2))/etaMax;

		a = anew;
		b = bnew;		
	}
	if(b>500)
	{
		etaMax = 1.0 / (1.0 + (b-1.0)/(a-1.0));
		bnew = 500;
		anew = (1+etaMax*(bnew-2))/(1-etaMax);

		a = anew;
		b = bnew;
	}
	if(a<=1e-4 && b>1.0)
	{
		return sourceFunc(sourceName, T_st, 0.0);
	}
	if(b<=1e-4 && a>1.0)
	{
		return sourceFunc(sourceName, T_st, 1.0);
	}

	a = max (a, 1e-9);
	b = max (b, 1e-9);


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

	double eta1,eta2;
	std::vector<double> eta , fBeta_num , fBeta_denom;
	int N=20;
	int M=50;
	int size;

	//Numerator
	double n1,n2,n3;
	n1 = sourceFunc(sourceName, T_st, 0.0) * std::pow(epsilon,a)/a;
	n2 = 0.0;
	n3 = sourceFunc(sourceName, T_st, 1.0) * std::pow(epsilon,b)/b;

	//Denominator
	double d1,d2,d3;
	d1 = std::pow(epsilon,a)/a;
	d2 = 0.0;
	d3 = std::pow(epsilon,b)/b;

	// setting up the bounded integral data between eta1 & eta2
	// eta has values between eta1 & eta2
	// fBeta_num is sourceFunc(eta)*Beta(eta)
	// fBeta_denom is Beta(eta)		
	for(int j=0 ; j<11 ; j++)
	{
		eta1=v[j];
		eta2=v[j+1];

		if(eta1==0.1)
		{
			// use M
			size = M;
		}
		else
		{
			// use N
			size = N;
		}	

		eta.assign(size+1, 0.0);
		fBeta_num.assign(size+1, 0.0);
		fBeta_denom.assign(size+1, 0.0);

		int i=0;
		eta[i]   		= eta1 + i* (eta2-eta1) / size; 
		fBeta_num[i]   	= sourceFunc(sourceName, T_st, eta[i]) * std::pow(eta[i],a-1.0) * std::pow(1.0-eta[i],b-1.0);					
		fBeta_denom[i] 	= std::pow(eta[i],a-1.0) * std::pow(1.0-eta[i],b-1.0);
		
		for(int i=1 ; i<size+1 ; i++)
		{
			eta[i]   		= eta1 + i* (eta2-eta1) / size; 
			fBeta_num[i]   	= sourceFunc(sourceName, T_st, eta[i]) * std::pow(eta[i],a-1.0) * std::pow(1.0-eta[i],b-1.0);					
			fBeta_denom[i] 	= std::pow(eta[i],a-1.0) * std::pow(1.0-eta[i],b-1.0);

			//integration step
			n2 += 0.5 * (eta[i]-eta[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
			d2 += 0.5 * (eta[i]-eta[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
		}
	}

	return (n1+n2+n3)/max(d1+d2+d3, 1e-6);
}

else
{
	return sourceFunc(sourceName, T_st, etaTilde);
}

}


