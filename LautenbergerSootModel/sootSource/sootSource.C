#include "sootSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
sootSource::sootSource()
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
sootSource::~sootSource()
{
}

sootSource::sootSource
(
                double Z_st_,
                double f_Zp_sf_,
                double MW_fuel_,
                double MW_air_,
                double T_flameAd_,
                double T_inf_,
                double T_air_,
                double rho_air_
):                
    Z_st(Z_st_),
    f_Zp_sf(f_Zp_sf_),
    MW_fuel(MW_fuel_),
    MW_air(MW_air_),
    T_flameAd(T_flameAd_),
    T_inf(T_inf_),
    T_air(T_air_),
    rho_air(rho_air_)
{

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

// A function to solve system of 4x4 linear equations
std::vector<double> sootSource::solveEquations(
                      double a00, double a01, double a02, double a03, double a04, 
                      double a10, double a11, double a12, double a13, double a14,
                      double a20, double a21, double a22, double a23, double a24,
                      double a30, double a31, double a32, double a33, double a34
                      )
{

    int i,j,k;

    int n=4;

    double mat[n][n+1];
    
    std::vector<double> res;
    res.assign(n,0.0);
   
    mat[0][0] = a00;
    mat[0][1] = a01;
    mat[0][2] = a02;
    mat[0][3] = a03;
    mat[0][4] = a04;

    mat[1][0] = a10;
    mat[1][1] = a11;
    mat[1][2] = a12;
    mat[1][3] = a13;
    mat[1][4] = a14; 

    mat[2][0] = a20;
    mat[2][1] = a21;
    mat[2][2] = a22;
    mat[2][3] = a23;
    mat[2][4] = a24;

    mat[3][0] = a30;
    mat[3][1] = a31;
    mat[3][2] = a32;
    mat[3][3] = a33;
    mat[3][4] = a34;


    for(i=0;i<n;i++) 
    {                   
        for(j=i+1;j<n;j++)
        {
            if(abs(mat[i][i]) < abs(mat[j][i]))
            {
                for(k=0;k<n+1;k++)
                {
        /* swapping mat[i][k] and mat[j][k] */
        mat[i][k]=mat[i][k]+mat[j][k];
                    mat[j][k]=mat[i][k]-mat[j][k];
                    mat[i][k]=mat[i][k]-mat[j][k];
                }
            }
      }
    }
   
     /* performing Gaussian elimination */
    for(i=0;i<n-1;i++)
    {
        for(j=i+1;j<n;j++)
        {
            double f=mat[j][i]/mat[i][i];
            for(k=0;k<n+1;k++)
            {
              mat[j][k]=mat[j][k]-f*mat[i][k];
            }
        }
    }
    /* Backward substitution for discovering values of unknowns */
    for(i=n-1;i>=0;i--)          
    {                     
        res[i]=mat[i][n];
                    
        for(j=i+1;j<n;j++)
        {
          if(i!=j)
            {
              res[i]=res[i]-mat[i][j]*res[j];
            }          
        }
    res[i]=res[i]/mat[i][i];  
    }

    return res;
}


// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

//calculate the polynomial coefficients
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void sootSource::calcCoeff()
{

    int n=4;

    //calculating constants of F_sf polynomial
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    psi_L_sf = 1.05;
    psi_P_sf = 1.77;
    psi_H_sf = 2.15;

    Z_L_sf = psi_L_sf * Z_st;
    Z_H_sf = psi_H_sf * Z_st;
    Z_P_sf = psi_P_sf * Z_st;


    c_sf.assign(n,0.0);

    c_sf = solveEquations(
                            1.0, Z_L_sf, std::pow(Z_L_sf,2.0), std::pow(Z_L_sf,3.0), 0.0, 
                            1.0, Z_H_sf, std::pow(Z_H_sf,2.0), std::pow(Z_H_sf,3.0), 0.0,
                            0.0, 1.0, 2*Z_P_sf, 3*std::pow(Z_P_sf,2.0), 0.0,
                            1.0, Z_P_sf, std::pow(Z_P_sf,2.0), std::pow(Z_P_sf,3.0), f_Zp_sf
                      );

    //calculating constants of G_sf polynomial
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    T_L = 1375;
    T_P = 1625;
    T_H = 1825;

    k_sf.assign(n,0.0);   
  
    k_sf = solveEquations(
                            1.0, T_L, std::pow(T_L,2.0), std::pow(T_L,3.0), 0.0, 
                            1.0, T_H, std::pow(T_H,2.0), std::pow(T_H,3.0), 0.0,
                            0.0, 1.0, 2*T_P, 3*std::pow(T_P,2.0), 0.0,
                            1.0, T_P, std::pow(T_P,2.0), std::pow(T_P,3.0), 1.0
                      );  


    //calculating constants of F_so polynomial
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    psi_L_so = 0.56;
    psi_P_so = 0.84;
    psi_H_so = 1.05;
    f_Zp_so  = 0.85;

    Z_L_so = psi_L_so * Z_st;
    Z_H_so = psi_H_so * Z_st;
    Z_P_so = psi_P_so * Z_st;

    c_so.assign(n,0.0);

    c_so = solveEquations(
                            1.0, Z_L_so, std::pow(Z_L_so,2.0), std::pow(Z_L_so,3.0), 0.0, 
                            1.0, Z_H_so, std::pow(Z_H_so,2.0), std::pow(Z_H_so,3.0), 0.0,
                            0.0, 1.0, 2*Z_P_so, 3*std::pow(Z_P_so,2.0), 0.0,
                            1.0, Z_P_so, std::pow(Z_P_so,2.0), std::pow(Z_P_so,3.0), f_Zp_so
                      );
}

// Soot formation function F_sf
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sootSource::F_sf(double Z)
{

    double f_sf = 0.0;

      if(Z>=Z_L_sf && Z<=Z_H_sf) //note that function is zero otherwise
      {
        f_sf = c_sf[0] + c_sf[1] * Z + c_sf[2] * std::pow(Z,2.0) + c_sf[3] * std::pow(Z,3.0); 
      }  
    
    double MW = MW_air * (1.0-Z) + MW_fuel * Z;

    return f_sf/(rho_air*(MW/MW_air));
}


// Soot formation function G_sf
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sootSource::G_sf(double Tstar)
{
    //returning the non-normalized temperature from T*
    double T = T_inf + Tstar * (T_flameAd - T_inf);
        
    double g_sf = 0.0;

      if(T>=T_L && T<=T_H) //note that function is zero otherwise
      {
        g_sf = k_sf[0] + k_sf[1] * T + k_sf[2] * std::pow(T,2.0) + k_sf[3] * std::pow(T,3.0); 
      }  

    return g_sf/(T_air/T);
}


// Soot oxidation function F_so
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sootSource::F_so(double Z)
{
     
    double f_so = 0.0;

      if(Z>=Z_L_so && Z<=Z_H_so) //note that function is zero otherwise
      {
        f_so = c_so[0] + c_so[1] * Z + c_so[2] * std::pow(Z,2.0) + c_so[3] * std::pow(Z,3.0); 
      }  
           
    double MW = MW_air * (1.0-Z) + MW_fuel * Z;

    return f_so/(rho_air*(MW/MW_air));
}


// Soot oxidation function g_so
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sootSource::G_so(double Tstar)
{
    //returning the non-normalized temperature from T*
    double T = T_inf + Tstar * (T_flameAd - T_inf);

    double g_so = 0.0;   


    g_so = max(0.0, 0.006*(T-1375));

    return g_so/(T_air/T);
}


// Inverse density functions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double sootSource::F_rho(double Z)
{
               
    double MW = MW_air * (1.0-Z) + MW_fuel * Z;

    return 1.0/(rho_air*(MW/MW_air));
}

double sootSource::G_rho(double Tstar)
{
    //returning the non-normalized temperature from T*
    double T = T_inf + Tstar * (T_flameAd - T_inf);

    return 1.0/(T_air/T);
}