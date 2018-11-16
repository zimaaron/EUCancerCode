// ///////////////////////////////////////////////////
// AOZ, work with LM and JW
// Oct 2018
// Template file for EU mort and incidence cancer model

// ///////////////////////////////////////////////////

// ///////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. Anything in the density namespace (ie SEPARABLE) returns the NEGATIVE log likelihood and is thus added to jnll accumulator
//    also, other density function such as dnorm and dbinom return POSITIVE log likelihood and are thus subtracted away.
// 3. 
// 4. 
// 5. refs 
// ///////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Robust Inverse Logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // cooresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}


// AR funtion from neal m
template<class Type>
SparseMatrix<Type> ar_Q(int N, Type rho, Type sigma) {
  SparseMatrix<Type> Q(N,N);
  Q.insert(0,0) = (1.) / pow(sigma, 2.);
  for (size_t n = 1; n < N; n++) {
    Q.insert(n,n) = (1. + pow(rho, 2.)) / pow(sigma, 2.);
    Q.insert(n-1,n) = (-1. * rho) / pow(sigma, 2.);
    Q.insert(n,n-1) = (-1. * rho) / pow(sigma, 2.);
  }
  Q.coeffRef(N-1,N-1) = (1.) / pow(sigma, 2.);
  return Q;
}

// Function for preparing a LCAR RE structure matrix given a neighborhood graph and rho parameter
// taken from Laura DL's small-area work
template<class Type> 
SparseMatrix<Type> lcar_strmat(SparseMatrix<Type> graph, Type rho) {
  SparseMatrix<Type> K = rho * graph; 
  for (size_t i = 0; i < K.rows(); i++)
    K.coeffRef(i,i) += (1 - rho);
  return K; 
}

// Function for constructing LCAR precision matrix
template<class Type> 
SparseMatrix<Type> lcar_prec(SparseMatrix<Type> I, SparseMatrix<Type> K, Type lambda, Type tau) {
  SparseMatrix<Type> Q = tau * ((1 - lambda) * I + lambda * K);
  return Q; 
}



// objective function (ie the likelihood function for the model), returns the evaluated negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{

  
  ///////////////////////////////////////////////////////////////////////////////
  // INPUTS
  ///////////////////////////////////////////////////////////////////////////////

  // TODO // options to calculate determinant only in outer loop
  // DATA_INTEGER(flag); // flag=0 => only prior

  // -- data for type 1 and 2 countries -- //
  DATA_INTEGER(NY);   // number of type 1 and 2 regions
  DATA_IVECTOR(popY); // population for the type 1 and 2 incidence cases
  DATA_VECTOR(Y);    // the InciYdence Cases, type 1 and 2
  DATA_VECTOR(Z);    // the Deaths, type 1 and 2
  DATA_IVECTOR(Yc);   // indices for type 1 and 2 countries
  DATA_IVECTOR(Ya);   // indices for type 1 and 2 countries
  DATA_IVECTOR(Yt);   // years for type 1 and 2 countries

  // -- data for type 2 and 3 countries -- //
  DATA_INTEGER(Ndeaths);   // number of type 2 and 3 countries
  DATA_IVECTOR(popDeaths); //population for the type 2 and 3 countries
  DATA_VECTOR(Deaths);    // Deatht outcome for type 2 and 3
  DATA_IVECTOR(Deaths_c);  //index for countries
  DATA_IVECTOR(Deaths_a);  //index for countries
  DATA_IVECTOR(Deaths_t);  //index for countries

  // -- general -- //
  DATA_INTEGER(Countries); // number of countries
  DATA_INTEGER(Ages);      // number of ages
  //  DATA_INTEGER(Years);    // number of years // TODO
  DATA_SPARSE_MATRIX(K);   // The structure DATA_MATRIX
  DATA_SPARSE_MATRIX(I);   // the Identity matrix
  // DATA_INTEGER(Nrates);    // total number of country-year-age
  DATA_IVECTOR(index_c);   // country index across country-year-age grid
  DATA_IVECTOR(index_a);   // age index across country-year-age grid
  DATA_IVECTOR(index_t);   // time index across country-year-age grid

  
  // -- Options -- //
  DATA_VECTOR(options);      // boolean vector of options to be used to select different models/modelling options:
  //                         // 0: Include priors. All are default settings right now
  //                         // 1: If 0, ADREPORT is on. Used for testing for now
  //                         // 2:

  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETERS
  ///////////////////////////////////////////////////////////////////////////////

  // -- spatial effects -- //
  PARAMETER_VECTOR(bI);   // LCAR Incidence
  PARAMETER_VECTOR(bMI);  // LCAR MI
  PARAMETER(lambdaI);     // proportion spatial // TODO add limits (0,1)
  PARAMETER(lambdaMI);    // proportion spatial // TODO add limits (0,1)
  PARAMETER(log_tau_bi);  // precision of bI    
  PARAMETER(log_tau_bmi); // precision of bMI   

  // -- linear time -- //
  PARAMETER(betI);           // intercept
  PARAMETER(betMI);          // intercept
  PARAMETER_VECTOR(betaI);   // random slope
  PARAMETER_VECTOR(betaMI);  // random slope
  PARAMETER(log_tau_betai);  // precision of bI  
  PARAMETER(log_tau_betami); // precision of bMI 

  // -- prior for space-age effect -- //
  PARAMETER_MATRIX(deltaI);   // LCAR Incidence
  PARAMETER_MATRIX(deltaMI);  // LCAR MI
  PARAMETER(log_tau_delti);  // precision of bI  
  PARAMETER(log_tau_deltmi); // precision of bMI 

  // -- age effects -- //
  PARAMETER_VECTOR(gammaI);  // RW2 Incidence
  PARAMETER_VECTOR(gammaMI); // RW2 MI
  PARAMETER(log_tau_gami);   // precision of bI  
  PARAMETER(log_tau_gammi);  // precision of bMI 

  ///////////////////////////////////////////////////////////////////////////////
  // TRANSFORMATIONS - DERIVED PARAMETERS
  ///////////////////////////////////////////////////////////////////////////////

  // -- lcar precision matrices -- //
  SparseMatrix<Type> TauI  = lcar_prec(I, K, lambdaI,  exp(log_tau_bi));
  SparseMatrix<Type> TauMI = lcar_prec(I, K, lambdaMI, exp(log_tau_bmi));

  // -- log and logit rate vecs -- //    
  vector<Type> log_muY(NY);  // log of poisson rate for Y
  vector<Type> logit_rY(NY); // logit of binomial prob/rate for Y

  vector<Type> log_muZ(Ndeaths);  // log of poisson rate for Y
  vector<Type> logit_rZ(Ndeaths); // logit of binomial prob/rate for Y
  vector<Type> lamb(Ndeaths);     // poisson lambda rate

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // initialize our joint negative log-likelihood
  Type jnll = 0;

  // IF NEEDED, PRINT THINGS FOR DEBUGGING
  // e.g.
  // printf("Epsilon_stz size: %ld \n", Epsilon_stz.size());
  // print parallel info
  // max_parallel_regions = omp_get_max_threads();
  // printf("This is thread %ld\n", max_parallel_regions);
  // max_parallel_regions = omp_get_max_threads();

  ///////////////////////////////////////////////////////////////////////////////
  // PRIORS
  ///////////////////////////////////////////////////////////////////////////////
  if(options[0]==1){

    // -- hyperpriors -- // 
    jnll -= dgamma(exp(log_tau_bi),      Type(1.0), Type(0.0001), true); 
    jnll -= dgamma(exp(log_tau_bmi),     Type(1.0), Type(0.0001), true); 
    jnll -= dgamma(exp(log_tau_betai),   Type(1.0), Type(0.0001), true); 
    jnll -= dgamma(exp(log_tau_betami),  Type(1.0), Type(0.0001), true); 
     
    jnll -= dgamma(exp(log_tau_delti),  Type(1.0), Type(0.0001), true); 
    jnll -= dgamma(exp(log_tau_deltmi), Type(1.0), Type(0.0001), true); 
      
    jnll -= dgamma(exp(log_tau_gami),    Type(1.0), Type(0.0001), true); 
    jnll -= dgamma(exp(log_tau_gammi),   Type(1.0), Type(0.0001), true); 
      
    for(int s = 0; s < Countries; s++){
      for(int t = 0; t < Ages; t++){
	//TODO unlist these into vector?
	jnll -= dnorm( deltaI(s, t), Type(0.0), 1/sqrt(exp(log_tau_delti)),  true);
	jnll -= dnorm(deltaMI(s, t), Type(0.0), 1/sqrt(exp(log_tau_deltmi)), true);
      }
    }

    // -- precision matrix for CAR -- //
    // WISHART // TODO

    // -- means for spline coefficients -- //
    jnll -= dnorm(betI,  Type(0.0), Type(sqrt(100000.0)), true);
    jnll -= dnorm(betMI, Type(0.0), Type(sqrt(100000.0)), true);

    jnll -= sum(dnorm(betaI,  betI,  1/sqrt(exp(log_tau_betai)),  true));
    jnll -= sum(dnorm(betaMI, betMI, 1/sqrt(exp(log_tau_betami)), true));
  }


  ///////////////////////////////////////////////////////////////////////////////
  // RANDOM EFFECTS
  ///////////////////////////////////////////////////////////////////////////////

  // -- CAR model using precision version of MVN -- //
  jnll += GMRF(TauI)(bI);   // TODO right now, normalization flag defaults to TRUE. 
  jnll += GMRF(TauMI)(bMI); // TODO could speed it up by setting to FALSE and normalizing in R like w/ SPDE

  for(int k = 0; k < 2; k++){
    jnll -= dnorm(gammaI(k),  Type(0.0), 1/sqrt(exp(log_tau_gami))  * 1000.0, true); // TODO why the 1000?
    jnll -= dnorm(gammaMI(k), Type(0.0), 1/sqrt(exp(log_tau_gammi)) * 1000.0, true);
  }

  for(int k = 2; k < Ages; k++){
    jnll -= dnorm(gammaI(k),  2 * gammaI(k - 1)  - gammaI(k - 2),  1/sqrt(exp(log_tau_gami)),  true);
    jnll -= dnorm(gammaMI(k), 2 * gammaMI(k - 1) - gammaMI(k - 2), 1/sqrt(exp(log_tau_gammi)), true);
  }

  ///////////////////////////////////////////////////////////////////////////////
  // DATA LIKELIHOOD
  ///////////////////////////////////////////////////////////////////////////////

  // -- type 1 and 2: incidence and mortality -- //
  for(int k = 0; k < NY; k++){

    //
    log_muY[k] = log(popY[k])
      + bI[Yc[k]]
      + gammaI[Ya[k]]
      + deltaI(Yc[k], Ya[k])
      + betaI[Yc[k]] * Yt[k];

    // 
    logit_rY[k] = bMI[Yc[k]]
      + gammaMI[Ya[k]]
      + deltaMI(Yc[k], Ya[k])
      + betaMI[Yc[k]] * Yt[k];
  }

  jnll -= sum(dpois(Y, exp(log_muY), true));

  jnll -= sum(dbinom_robust(Z, Y, logit_rY, true));

  
  // -- Type 2 and 3: Mortality model -- //
  for(int i = 0; i < Ndeaths; i++){

    //
    log_muZ[i] = log(popDeaths[i])
      + bI[Deaths_c[i]]
      + gammaI[Deaths_a[i]]
      + deltaI(Deaths_c[i], Deaths_a[i])
      + betaI[Deaths_c[i]] * Deaths_t[i];

    //
    logit_rZ[i] = bMI[Deaths_c[i]]
      + gammaMI[Deaths_a[i]]
      + deltaMI(Deaths_c[i], Deaths_a[i])
      + betaMI[Deaths_c[i]] * Deaths_t[i];

    //
    lamb[i] = exp(log_muZ[i]) * 1 / (1 + exp(-logit_rZ[i]));
  }

  jnll -= sum(dpois(Deaths, lamb, true));


  // // Report estimates
  // if(options[1] == 0){
  //   ADREPORT(alpha_j);
  //   ADREPORT(Epsilon_stz);
  // }

  
  return jnll;
}
