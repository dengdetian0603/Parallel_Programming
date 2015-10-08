/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Alassoprob->ny KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
// Optimization problem
//      \arg \min_x ||Ax-y|| + \lambda |x|_1
//

#include "common.h"
#include <random>

// Problem definition
shotgun_data * lassoprob;


// Major optimization: always keep updated vector 'Ax'
 
void initialize_feature(int feat_idx) {
    sparse_array& col = lassoprob->A_cols[feat_idx];
    feature& feat = lassoprob->feature_consts[feat_idx];
    lassoprob->x[feat_idx] = 0.0;

    // Precompute covariance of a feature
    feat.covar = 0;

    for(int i=0; i<col.length(); i++) {
        feat.covar += sqr(col.values[i]);
    }
    feat.covar *= 2;
    
    // Precompute (Ay)_i
    feat.Ay_i = 0.0;
    for(int i=0; i<col.length(); i++) {
        feat.Ay_i += col.values[i]*lassoprob->y[col.idxs[i]];
    }
    feat.Ay_i *= 2;

}

void initialize() {
    lassoprob->feature_consts.reserve(lassoprob->nx);
    lassoprob->x.resize(lassoprob->nx);
    lassoprob->Ax.resize(lassoprob->ny);
    #pragma omp for
    for(int i=0; i<lassoprob->nx; i++) {
        initialize_feature(i);
    }
}

valuetype_t soft_thresholdO(valuetype_t _lambda, valuetype_t shootDiff) {
    return (shootDiff > _lambda)* (_lambda - shootDiff) + 
	               (shootDiff < -_lambda) * (-_lambda- shootDiff) ;
}

valuetype_t soft_threshold(valuetype_t _lambda, valuetype_t shootDiff) {
  if (shootDiff > _lambda) return _lambda - shootDiff;
  if (shootDiff < -_lambda) return -_lambda - shootDiff ;
  else return 0;
}


double GetUniform()
{
    static std::default_random_engine re;
    static std::uniform_real_distribution<double> Dist(0,1);
    return Dist(re);
}

// John D. Cook, http://stackoverflow.com/a/311716/15485
void SampleWithoutReplacement
(
    int populationSize,    // size of set sampling from
    int sampleSize,        // size of each sample
    std::vector<int> & samples  // output, zero-offset indicies to selected items
)
{
    // Use Knuth's variable names
    int& n = sampleSize;
    int& N = populationSize;

    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    while (m < n)
    {
        u = GetUniform(); // call a uniform(0,1) random number generator

        if ( (N - t)*u >= n - m )
        {
            t++;
        }
        else
        {
            samples[m] = t;
            t++; m++;
        }
    }
}


double shoot(int x_i, valuetype_t lambda) {
    feature& feat = lassoprob->feature_consts[x_i];
    valuetype_t oldvalue = lassoprob->x[x_i];
    
    // Compute dotproduct A'_i*(Ax)
    valuetype_t AtAxj = 0.0;
    sparse_array& col = lassoprob->A_cols[x_i];
    int len=col.length();
    for(int i=0; i<len; i++) {
        AtAxj += col.values[i] * lassoprob->Ax[col.idxs[i]];// the second term explodes
    }
    
    valuetype_t S_j = 2 * AtAxj - feat.covar * oldvalue - feat.Ay_i;
    valuetype_t newvalue = soft_threshold(lambda,S_j)/feat.covar;
    valuetype_t Delta = (newvalue - oldvalue);  // delta x of x_i th coordinate
    
    // Update ax
    if (Delta != 0.0) 
    {
        // for(int i=0; i<len; i++) 
        // {
	       // lassoprob->Ax.add(col.idxs[i], col.values[i] * Delta); 
        //    // lassoprob->Ax is accessed by different thread at the same time
        //    // update the collective coordinate descent step
        // }
        lassoprob->x[x_i] = newvalue;
    }

    if (x_i%lassoprob->nx==0)
    {
        mexPrintf("Ax[10] = %f\n",lassoprob->Ax[col.idxs[10]]);
        //mexPrintf("AtAxj = %f\n", AtAxj);
        mexPrintf("delta %d = %f\n", x_i, Delta);
        
    }
    //return std::abs(Delta);
    return Delta;
}

double shoot1(int x_i, valuetype_t Lambda) {
    feature& feat = lassoprob->feature_consts[x_i];
    valuetype_t oldvalue = lassoprob->x[x_i];
    
    // Compute dotproduct A'_i*(Ax)
    valuetype_t AtAxj = 0.0;
    sparse_array& col = lassoprob->A_cols[x_i];
    int len=col.length();
    for(int i=0; i<len; i++) {
        AtAxj += col.values[i] * lassoprob->Ax[col.idxs[i]];
    }
    
    valuetype_t S_j = 2 * AtAxj - feat.covar * oldvalue - feat.Ay_i;
    valuetype_t newvalue = soft_threshold(Lambda,S_j)/feat.covar;
    valuetype_t Delta = (newvalue - oldvalue);  // try debug
    
    // Update ax
    if (Delta != 0.0) {
        for(int i=0; i<len; i++) {
           lassoprob->Ax.add(col.idxs[i], col.values[i] * Delta);
        }
        lassoprob->x[x_i] = newvalue;
    }
    return std::abs(Delta);
}

double shoot2(int x_i, valuetype_t Lambda) {
    feature& feat = lassoprob->feature_consts[x_i];
    valuetype_t oldvalue = lassoprob->x[x_i];
    
    // Compute dotproduct A'_i*(Ax)
    valuetype_t AtAxj = 0.0;
    sparse_array& col = lassoprob->A_cols[x_i];
    int len=col.length();
    for(int i=0; i<len; i++) {
        AtAxj += col.values[i] * lassoprob->Ax[col.idxs[i]];
    }
    
    valuetype_t S_j = 2 * AtAxj - feat.covar * oldvalue - feat.Ay_i;
    valuetype_t newvalue = soft_threshold(Lambda,S_j)/feat.covar;
    valuetype_t Delta = (newvalue - oldvalue);  // try debug
    
    // Update ax
    if (Delta != 0.0) {
        #pragma omp parallel for
        for(int i=0; i<len; i++) {
           lassoprob->Ax.add(col.idxs[i], col.values[i] * Delta);
        }
        lassoprob->x[x_i] = newvalue;
    }

    return std::abs(Delta);
}


// Find such lambda that if used for regularization,
// optimum would have all weights zero.
valuetype_t compute_max_lambda() {
    valuetype_t maxlambda = 0.0;
    for(int i=0; i<lassoprob->nx; i++) {
        maxlambda = std::max(maxlambda, std::abs(lassoprob->feature_consts[i].Ay_i));
    }
    return maxlambda;
}

valuetype_t get_term_threshold(int k, int K, double delta_threshold) {
  // Stricter termination threshold for last step in the optimization.
  return (k == 0 ? delta_threshold  : (delta_threshold + k*(delta_threshold*50)/K));
}

valuetype_t compute_objective(valuetype_t _lambda, std::vector<valuetype_t>& x, valuetype_t * l1x = NULL, valuetype_t * l2err = NULL) {
    double least_sqr = 0.0;
    #pragma omp parallel for reduction(+: least_sqr) 
    for(int i=0; i<lassoprob->ny; i++) 
    {
        least_sqr += (lassoprob->Ax[i]-lassoprob->y[i])*(lassoprob->Ax[i]-lassoprob->y[i]);					
    }
     // Penalty check for feature 0
    double penalty = 0.0;
    #pragma omp parallel for reduction(+: penalty) 
    for(int i=0; i<lassoprob->nx; i++) 
    {
        penalty += std::abs(lassoprob->x[i]);
    }
    //
    if (l1x != NULL) *l1x = penalty;
    if (l2err != NULL) *l2err = least_sqr;
    return penalty * _lambda + least_sqr;
}

void main_optimization_loop(double lambda, int regpathlength, double threshold, int maxiter, int verbose, int numthreads) {
    // Empirically found heuristic for deciding how malassoprob->ny steps
    // to take on the regularization path
    int regularization_path_length = 0;
    //mexPrintf("regularization_path_length = %d\n",regularization_path_length);

    //valuetype_t lambda_max = compute_max_lambda();
    valuetype_t lambda_min = lambda;
    //valuetype_t alpha = pow(lambda_max/lambda_min, 1.0/(1.0*regularization_path_length));
    int regularization_path_step = regularization_path_length;

    double delta_threshold = threshold;
    //long long int num_of_shoots = 0;
    //int counter = 0;
    int iterations = 0;
    double *delta = new double[lassoprob->nx];

    #pragma omp parallel for  
    for(int i=0; i<lassoprob->nx; i++) 
    {
        delta[i] = 1.0;
    }

    double l1x = 0, l2err = 0;
    valuetype_t obj0 = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
    mexPrintf("Initial Objective: %g L1 : %g L2err: %g\n", obj0, l1x, l2err);

    do {
        ++iterations;
        if (iterations >= maxiter && maxiter > 0) {
            mexPrintf("Exceeded max iterations: %d", maxiter);
            valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
            mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
            return;
        }
        double maxChange=0.0;
        //lambda = lambda_min * pow(alpha, regularization_path_step);
        
        
        std::vector<int> subCols(numthreads);
        SampleWithoutReplacement(lassoprob->nx, numthreads, subCols);

        // Parallel loop
        #pragma omp parallel for  
        for(int i=0; i<numthreads; i++) 
        {
            delta[subCols[i]] = shoot(subCols[i], lambda);
            //mexPrintf("thread ID: %d, coordinate: %d \n",omp_get_thread_num(),i);
            //mexPrintf("Coordinate: %d, delta : %f \n",subCols[i],delta[subCols[i]]);
        }
        #pragma omp barrier

        #pragma omp parallel for 
        for(int i=0; i < lassoprob->ny; i++) 
        {
            for(int j=0; j< numthreads; j++) 
            {
                sparse_array& Acol = lassoprob->A_cols[subCols[j]];      
                lassoprob->Ax.add(Acol.idxs[i], Acol.values[i] * delta[subCols[j]]); 
                // update the collective coordinate descent step
            }
        }

        // TODO: use OpenMP reductions
        maxChange = 0.0;
        for(int i=0; i<lassoprob->nx; i++)
        {
                maxChange = (maxChange < abs(delta[i]) ? abs(delta[i]) : maxChange);
        }

        if (iterations%1000==0){
            mexPrintf("iteration: %d \n", iterations);
            mexPrintf("maxChange = %f \n",maxChange);
        }

        //num_of_shoots += lassoprob->nx;
        //counter++;
        // Convergence check.
        bool converged1 = (maxChange <= get_term_threshold(regularization_path_step,regularization_path_length,delta_threshold));
        if (converged1) {     //|| counter>std::min(100, (100-regularization_path_step)*2)) {
            //counter = 0;
            //mexPrintf("Checking convergence cyclicaly...\n");
            //for(int i=0; i<lassoprob->nx; i++) 
            //{
            //   delta[i] = shoot2(i, lambda);
            //}

            // cyclic check
            // maxChange = 0.0;
            // for(int i=0; i<lassoprob->nx; i++)
            // {
            //         maxChange = (maxChange < delta[i] ? delta[i] : maxChange);
            // }
            // if (maxChange <= delta_threshold){
                valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
                mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
                mexPrintf("maxChange = %f \n",maxChange);
                regularization_path_step--; 
            //}
        }
        if (verbose){
          valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
          mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
        }
    } while (regularization_path_step >= 0);
    if (!verbose){
        double l1x = 0, l2err = 0;
        valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
        mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
    }
    delete[] delta;
    //mexPrintf("Num of shoots = %lld\n", num_of_shoots);
    mexPrintf("Num of interrations = %d\n", iterations);
}

void main_optimization_loop1(double lambda, int regpathlength, double threshold, int maxiter, int verbose, int numthreads) {
    // Empirically found heuristic for deciding how malassoprob->ny steps
    // to take on the regularization path
    int regularization_path_length = 0;
    //mexPrintf("regularization_path_length = %d\n",regularization_path_length);

    //valuetype_t lambda_max = compute_max_lambda();
    valuetype_t lambda_min = lambda;
    //valuetype_t alpha = pow(lambda_max/lambda_min, 1.0/(1.0*regularization_path_length));
    int regularization_path_step = regularization_path_length;

    double delta_threshold = threshold;
    //long long int num_of_shoots = 0;
    //int counter = 0;
    int iterations = 0;
    double *delta = new double[lassoprob->nx];

    double l1x = 0, l2err = 0;
    valuetype_t obj0 = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
    mexPrintf("Initial Objective: %g L1 : %g L2err: %g\n", obj0, l1x, l2err);

    do {
        ++iterations;
        if (iterations >= maxiter && maxiter > 0) {
            mexPrintf("Exceeded max iterations: %d", maxiter);
            valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
            mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
            return;
        }
        double maxChange=0.0;
        //lambda = lambda_min * pow(alpha, regularization_path_step);

        // Parallel loop
        #pragma omp parallel for  
        for(int i=0; i<lassoprob->nx; i++) 
        {
            delta[i] = shoot1(i, lambda);
        }

        // TODO: use OpenMP reductions
        maxChange = 0.0;
        for(int i=0; i<lassoprob->nx; i++)
        {
                maxChange = (maxChange < delta[i] ? delta[i] : maxChange);
        }

        if (iterations%1000==0){
            mexPrintf("iteration: %d \n", iterations);
            mexPrintf("maxChange = %f \n",maxChange);
        }

        //num_of_shoots += lassoprob->nx;
        //counter++;
        // Convergence check.
        bool converged1 = (maxChange <= get_term_threshold(regularization_path_step,regularization_path_length,delta_threshold));
        if (converged1) {     //|| counter>std::min(100, (100-regularization_path_step)*2)) {
            //counter = 0;
            valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
            mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
            regularization_path_step--; 
        }
        if (verbose){
          valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
          mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
        }
    } while (regularization_path_step >= 0);
    if (!verbose){
        double l1x = 0, l2err = 0;
        valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
        mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
    }
    delete[] delta;
    mexPrintf("Num of interrations = %d\n", iterations);
}

void main_optimization_loop2(double lambda, int regpathlength, double threshold, int maxiter, int verbose, int numthreads) {
    // Empirically found heuristic for deciding how malassoprob->ny steps
    // to take on the regularization path
    int regularization_path_length = 0;
    //mexPrintf("regularization_path_length = %d\n",regularization_path_length);

    //valuetype_t lambda_max = compute_max_lambda();
    valuetype_t lambda_min = lambda;
    //valuetype_t alpha = pow(lambda_max/lambda_min, 1.0/(1.0*regularization_path_length));
    int regularization_path_step = regularization_path_length;

    double delta_threshold = threshold;
    //long long int num_of_shoots = 0;
    //int counter = 0;
    int iterations = 0;
    double *delta = new double[lassoprob->nx];

    double l1x = 0, l2err = 0;
    valuetype_t obj0 = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
    mexPrintf("Initial Objective: %g L1 : %g L2err: %g\n", obj0, l1x, l2err);

    do {
        ++iterations;
        if (iterations >= maxiter && maxiter > 0) {
            mexPrintf("Exceeded max iterations: %d", maxiter);
            valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
            mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
            return;
        }
        double maxChange=0.0;
        //lambda = lambda_min * pow(alpha, regularization_path_step);

        // Parallel loop 
        for(int i=0; i<lassoprob->nx; i++) 
        {
            delta[i] = shoot2(i, lambda);
        }

        // TODO: use OpenMP reductions
        maxChange = 0.0;
        for(int i=0; i<lassoprob->nx; i++)
        {
                maxChange = (maxChange < delta[i] ? delta[i] : maxChange);
        }

        if (iterations%1000==0){
            mexPrintf("iteration: %d \n", iterations);
            mexPrintf("maxChange = %f \n",maxChange);
        }

        //num_of_shoots += lassoprob->nx;
        //counter++;
        // Convergence check.
        bool converged1 = (maxChange <= get_term_threshold(regularization_path_step,regularization_path_length,delta_threshold));
        if (converged1) {     //|| counter>std::min(100, (100-regularization_path_step)*2)) {
            //counter = 0;
            valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
            mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
            mexPrintf("maxChange = %f \n",maxChange);
            regularization_path_step--; 
        }
        if (verbose){
          valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
          mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
        }
    } while (regularization_path_step >= 0);
    if (!verbose){
        double l1x = 0, l2err = 0;
        valuetype_t obj = compute_objective(lambda, lassoprob->x, &l1x, &l2err);
        mexPrintf("Objective: %g L1 : %g L2err: %g\n", obj, l1x, l2err);
    }
    delete[] delta;
    //mexPrintf("Num of shoots = %lld\n", num_of_shoots);
    mexPrintf("Num of interrations = %d\n", iterations);
}





double solveLasso(shotgun_data  * probdef, double lambda, int regpathlength, double threshold, int maxiter, int verbose, int numthreads) {
    lassoprob = probdef;
    initialize();
    if (regpathlength == 0){ // pSCD
        main_optimization_loop(lambda, regpathlength, threshold, maxiter, verbose, numthreads); 
    } else if (regpathlength == 1){ // sCCD
        main_optimization_loop1(lambda, regpathlength, threshold, maxiter, verbose, numthreads); 
    } else if (regpathlength == 2){ // uCCD
        main_optimization_loop2(lambda, regpathlength, threshold, maxiter, verbose, numthreads); 
    } else{
        throw std::runtime_error("Unknown algorithm solver");
    }
    return 0;
}


