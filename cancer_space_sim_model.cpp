// ///////////////////////////////////////////////////
// AOZ, work with LM and JW
// Oct 2018
// single year, single age, space only model for use with simulated data
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
  Q.coeffRef(N-1, N-1) = (1.) / pow(sigma, 2.);
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

  // -- data for type 1, 2 and 3 countries -- //
  DATA_INTEGER(N);       // number of observations
  DATA_VECTOR(Ctype);    // vector of which country is which type of observations
  DATA_VECTOR(pop_nat);  // national population
  DATA_VECTOR(pop_reg);  // regional population
  DATA_VECTOR(r_cases);  // 
  DATA_VECTOR(r_deaths); // 
  DATA_VECTOR(l_cases);  // 
  DATA_VECTOR(l_deaths); // 

   // -- Options -- //
  DATA_VECTOR(options);      // boolean vector of options to be used to select different models/modelling options:
  //                         // 0: Include priors. All are default settings right now
  //                         // 1: If 0, ADREPORT is on. Used for testing for now
  
  // -- spatial effect matrices -- //
  DATA_SPARSE_MATRIX(K);   // The structure DATA_MATRIX
  DATA_SPARSE_MATRIX(I);   // the Identity matrix
   
  
  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETERS
  ///////////////////////////////////////////////////////////////////////////////

  // -- intercept -- //
  PARAMETER(aI);  // intercept
  PARAMETER(aMI); // intercept

  // -- country random effects -- //
  PARAMETER_VECTOR(vI);   // random intercept
  PARAMETER_VECTOR(vMI);  // random intercept
  PARAMETER(log_tau_vI);  // precision of bI  
  PARAMETER(log_tau_vMI); // precision of bMI

  // -- spatial effects -- //
  PARAMETER_VECTOR(uI);   // LCAR spatial effect for incidence
  PARAMETER_VECTOR(uMI);  // LCAR spatial effect for MI
  PARAMETER(lambdaI);     // proportion spatial // TODO add limits (0,1)
  PARAMETER(lambdaMI);    // proportion spatial // TODO add limits (0,1)
  PARAMETER(log_tau_uI);  // precision of uI    
  PARAMETER(log_tau_uMI); // precision of uMI
 
  ///////////////////////////////////////////////////////////////////////////////
  // TRANSFORMATIONS - DERIVED PARAMETERS
  ///////////////////////////////////////////////////////////////////////////////

  // -- rates and whatnot -- //    
  vector<Type> pc(N); // incidence rate
  vector<Type> rc(N); // mortality prob | cancer
  vector<Type> qc(N); // inc * mort

  // -- country random effect objects -- // TODO! get correct model from Rue 2005
  // -- sum-to-0 constraints -- //
  Type vI_sum  = sum(vI);
  Type vMI_sum = sum(vMI);
  
  // -- space relevant objects -- //
  // -- sum-to-0 constraints -- // TODO! get correct model from Rue 2005
  Type uI_sum  = sum(uI);
  Type uMI_sum = sum(uMI);
    
  // -- lcar precision matrices -- //
  SparseMatrix<Type> TauI  = lcar_prec(I, K, lambdaI,  exp(log_tau_uI));
  SparseMatrix<Type> TauMI = lcar_prec(I, K, lambdaMI, exp(log_tau_uMI));

  
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // initialize our joint negative log-likelihood
  vector<Type> jnll(3);
  // jnll(0): hyperpriors
  // jnll(1): priors
  // jnll(2): data liklihood

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

    // ----------------- //
    // -- hyperpriors -- //
    // ----------------- //

    // - country random effects - //
    if(options[1] == 1){
      jnll(0) -= dgamma(exp(log_tau_vI),  Type(1.0), Type(0.0001), true); 
      jnll(0) -= dgamma(exp(log_tau_vMI), Type(1.0), Type(0.0001), true);
    }

    // - spatial random effects - //
    if(options[2] == 1){
      jnll -= dgamma(exp(log_tau_uI),  Type(1.0), Type(0.0001), true); 
      jnll -= dgamma(exp(log_tau_uMI), Type(1.0), Type(0.0001), true);

      // priors for lambda
      // FIXME maybe better to use transformed lambda on (-Inf, Inf)?
      jnll -= dbeta(lambdaI, Type(1.0), Type(1.0), true);
      jnll -= dbeta(lambdaI, Type(1.0), Type(1.0), true);
    }

    // ------------ //
    // -- priors -- //
    // ------------ //
 
    // -- priors for intercepts -- //
    jnll(1) -= dnorm(aI,  Type(0.0), Type(sqrt(100000.0)), true);
    jnll(1) -= dnorm(aMI, Type(0.0), Type(sqrt(100000.0)), true);

     // - country random effects - //
    if(options[1] == 1){
      jnll(1) -= sum(dnorm(vI,  Type(0.0), 1/sqrt(exp(log_tau_vI)),  true));
      jnll(1) -= sum(dnorm(vMI, Type(0.0), 1/sqrt(exp(log_tau_vMI)), true));

      // -- constraints -- //
      jnll(1) -= dnorm(vI_sum,  Type(0.0), Type(.0000001), true);
      jnll(1) -= dnorm(vMI_sum, Type(0.0), Type(.0000001), true);
    }

    // -- spatial random effects -- //
    if(options[2] == 1){
      jnll += GMRF(TauI)(uI);
      jnll += GMRF(TauMI)(uMI);
      
      // -- constraints -- //
      jnll -= dnorm(uI_sum,  Type(0.0), Type(.0000001), true);
      jnll -= dnorm(uMI_sum, Type(0.0), Type(.0000001), true);

      // -- correlated precision matrix for CAR -- //
      // WISHART // TODO?
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  // DATA LIKELIHOOD
  ///////////////////////////////////////////////////////////////////////////////

  // -- generate incidence, mort, and (inc * mort) vectors -- //
  for(int i = 0; i < N; i++){

    // cancer incidence baseline
    pc[i] = aI;

    // mortality prob | cancer (i.e. mi-ratio) baseline
    rc[i] = aMI;

    // if country REs
    if(options[1] == 1){
      pc[i] += vI[i];
      rc[i] += vMI[i];
    }

    // if spatial effect
    if(options[2] == 1){
      pc[i] += uI[i];
      rc[i] += uMI[i];
    }

    // convert with link functions
    pc[i] = exp(pc[i]);
    rc[i] = exp(rc[i]) / (1 + exp(rc[i]));

    // incidence * mort
    qc[i] = pc[i] * rc[i]; 
 
  }

  // -- data contribution to jnll(2) -- //
  for(int i = 0; i < N; i++){

    // -- type 1 countries -- //
    if(Ctype[i] == 1){

      // incidence
      jnll(2) -= dpois(r_cases[i], pop_nat[i] * pc[i], true);

      // mortality
      jnll(2) -= dbinom(r_deaths[i], r_cases[i], rc[i], true);
      
    }
    
    // -- type 2 countries -- //
    if(Ctype[i] == 2){

      // local/registry incidence
      jnll(2) -= dpois(l_cases[i], pop_reg[i] * pc[i], true);

      // local/registry mortality
      jnll(2) -= dbinom(l_deaths[i], l_cases[i], rc[i], true);

      // remainder mortality
      jnll(2) -= dpois(r_deaths[i], (pop_nat[i] - pop_reg[i]) * qc[i], true);
    }
    
    
    // -- type 3 countries -- //
    if(Ctype[i] == 3){

      // national mortality
      jnll(2) -= dpois(r_deaths[i], pop_nat[i] * qc[i], true);
    }

    // -- type 4 countries have no data -- //
    
  }
  
  
  // // Report estimates
  if(options[1] == 1){
    ADREPORT(pc); // just an example
  }

  printf("\n\nHyper contrib: %f \n", asDouble(jnll(0)));
  printf("Prior contrib: %f \n", asDouble(jnll(1)));
  printf("Data  contrib: %f \n\n", asDouble(jnll(2)));

  Type jnllSum  = sum(jnll);
  
  // return the jnll
  return jnllSum;
}
