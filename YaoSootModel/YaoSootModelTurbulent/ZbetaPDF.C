#include "ZbetaPDF.H"


double Func(const double& Z, const double& lower, const double& upper)
{
	if (Z>=lower && Z <=upper)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}


double Zbeta(
			const double& Ztilde,
			const double& Zvar,
			const double& lower,
			const double& upper
			)
{

	//beta function parameters
	double a = Ztilde * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);
	double b = (1.0-Ztilde) * (Ztilde*(1.0-Ztilde)/max(Zvar,1e-6) - 1.0);

	//bounds of a&b
	a = max(a,1e-6);
	b = max(b,1e-6);

	double Zmax = 1.0/(1.0+(b-1.0)/(a-1.0));
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

	// integration variables
	double epsilon = 1E-6;
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
	int M=50;
	int size;

	//Numerator
	double n1,n2,n3;
	n1 = Func(0.0, lower, upper) * std::pow(epsilon,a)/a;
	n2 = 0.0;
	n3 = Func(1.0, lower, upper) * std::pow(epsilon,b)/b;

	//Denominator
	double d1,d2,d3;
	d1 = std::pow(epsilon,a)/a;
	d2 = 0.0;
	d3 = std::pow(epsilon,b)/b;

	// setting up the bounded integral data between Z1 & Z2
	// ZZ has values between Z1 & Z2
	// fBeta_num is Func(Z)*Beta(Z)
	// fBeta_denom is Beta(Z)		
	for(int j=0 ; j<11 ; j++)
	{
		Z1=v[j];
		Z2=v[j+1];

		if((Z1==0.05) | (Z1==0.15))
		{
			// use M
			size = M;
		}
		else
		{
			// use N
			size = N;
		}	

		ZZ.assign(size+1, 0.0);
		fBeta_num.assign(size+1, 0.0);
		fBeta_denom.assign(size+1, 0.0);

		int i=0;
		ZZ[i]   		= Z1 + i* (Z2-Z1) / size; 
		fBeta_num[i]   	= Func(ZZ[i], lower, upper) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
		fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);
		
		for(int i=1 ; i<size+1 ; i++)
		{
			ZZ[i]   		= Z1 + i* (Z2-Z1) / size; 
			fBeta_num[i]   	= Func(ZZ[i], lower, upper) * std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);					
			fBeta_denom[i] 	= std::pow(ZZ[i],a-1.0) * std::pow(1.0-ZZ[i],b-1.0);

			//integration step
			n2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_num[i]   + fBeta_num[i-1]);
			d2 += 0.5 * (ZZ[i]-ZZ[i-1]) * (fBeta_denom[i] + fBeta_denom[i-1]);
		}
	}

	return (n1+n2+n3)/(d1+d2+d3);

}





