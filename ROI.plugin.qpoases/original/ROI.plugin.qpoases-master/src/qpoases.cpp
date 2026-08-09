
#include "qpoases.h"

/* HessianType 
    #define HST_ZERO             0
    #define HST_IDENTITY         1
    #define HST_POSDEF           2
    #define HST_POSDEF_NULLSPACE 3
    #define HST_SEMIDEF          4
    #define HST_INDEF            5
    #define HST_UNKNOWN          6
*/


// Solving Simply Bounded QPs
// [[Rcpp::export]]
SEXP qproblemb(int number_of_varibales, int hessian_type, int alloc_dense_matrix) {

    qpOASES::BooleanType alloc_dense_mat = (alloc_dense_matrix < 1) ? qpOASES::BT_FALSE : qpOASES::BT_TRUE;
    qpOASES::HessianType ht = static_cast<qpOASES::HessianType>(hessian_type);
    Rcpp::XPtr< qpOASES::QProblemB > p(new qpOASES::QProblemB(number_of_varibales, ht, alloc_dense_mat), true);
    return p;
}


// [[Rcpp::export]]
SEXP qproblem(int number_of_varibales, int number_of_constraints, int hessian_type) {

    qpOASES::HessianType ht = static_cast<qpOASES::HessianType>(hessian_type);
    XPtr< qpOASES::QProblem > p(new qpOASES::QProblem(number_of_varibales, number_of_constraints, ht), true);
    return p;
}


// [[Rcpp::export]]
SEXP sqproblem(int number_of_varibales, int number_of_constraints, int hessian_type, int alloc_dense_matrix) {

    qpOASES::BooleanType alloc_dense_mat = (alloc_dense_matrix < 1) ? qpOASES::BT_FALSE : qpOASES::BT_TRUE;
    qpOASES::HessianType ht = static_cast<qpOASES::HessianType>(hessian_type);
    Rcpp::XPtr< qpOASES::SQProblem > p(new qpOASES::SQProblem(number_of_varibales, number_of_constraints, ht, alloc_dense_mat), true);
    return p;
}


// [[Rcpp::export]]
SEXP set_options(SEXP r_model, Rcpp::List control) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    qpOASES::Options options;

    options.printLevel = static_cast<qpOASES::PrintLevel>((int) control["printLevel"]);

    options.enableRamping = static_cast<qpOASES::BooleanType>((int) control["enableRamping"]);
    options.enableFarBounds = static_cast<qpOASES::BooleanType>((int) control["enableFarBounds"]);
    options.enableFlippingBounds = static_cast<qpOASES::BooleanType>((int) control["enableFlippingBounds"]);
    options.enableRegularisation = static_cast<qpOASES::BooleanType>((int) control["enableRegularisation"]);
    options.enableFullLITests = static_cast<qpOASES::BooleanType>((int) control["enableFullLITests"]);
    options.enableNZCTests = static_cast<qpOASES::BooleanType>((int) control["enableNZCTests"]);

    options.enableDriftCorrection = (int) control["enableDriftCorrection"];
    options.enableCholeskyRefactorisation = (int) control["enableCholeskyRefactorisation"];
    options.enableEqualities = static_cast<qpOASES::BooleanType>((int) control["enableEqualities"]);

    options.terminationTolerance = (double) control["terminationTolerance"];
    options.boundTolerance = (double) control["boundTolerance"];

    options.boundRelaxation = (double) control["boundRelaxation"];
    options.epsNum = (double) control["epsNum"];
    options.epsDen = (double) control["epsDen"];
    options.maxPrimalJump = (double) control["maxPrimalJump"];
    options.maxDualJump = (double) control["maxDualJump"];

    options.initialRamping = (double) control["initialRamping"];
    options.finalRamping = (double) control["finalRamping"];
    options.initialFarBounds = (double) control["initialFarBounds"];
    options.growFarBounds = (double) control["growFarBounds"];
    options.rcondSMin = (double) control["rcondSMin"];

    options.epsFlipping = (double) control["epsFlipping"];
    options.epsRegularisation = (double) control["epsRegularisation"];
    options.epsIterRef = (double) control["epsIterRef"];
    options.epsLITests = (double) control["epsLITests"];
    options.epsNZCTests = (double) control["epsNZCTests"];

    options.numRegularisationSteps = (int) control["numRegularisationSteps"];
    options.numRefinementSteps = (int) control["numRefinementSteps"];
    options.dropBoundPriority = (int) control["dropBoundPriority"];
    options.dropEqConPriority = (int) control["dropEqConPriority"];
    options.dropIneqConPriority = (int) control["dropIneqConPriority"];

    options.enableInertiaCorrection = static_cast<qpOASES::BooleanType>((int) control["enableInertiaCorrection"]);
    options.enableDropInfeasibles = static_cast<qpOASES::BooleanType>((int) control["enableDropInfeasibles"]);

    options.initialStatusBounds = static_cast<qpOASES::SubjectToStatus>((int) control["initialStatusBounds"]);

    model->setOptions( options );
    return R_NilValue;
}


/***
    H                           const real_t* const         Hessian matrix (a shallow copy is made).
                                                            If Hessian matrix is trivial, a NULL pointer can be passed.
    g                           const real_t* const         Gradient vector.
    A                           const real_t* const         Constraint matrix (a shallow copy is made).
    lb                          const real_t* const         Lower bound vector (on variables).
                                                            If no lower bounds exist, a NULL pointer can be passed.
    ub                          const real_t* const         Upper bound vector (on variables).
                                                            If no upper bounds exist, a NULL pointer can be passed.
    lbA                         const real_t* const         Lower constraints' bound vector.
                                                            If no lower constraints' bounds exist, a NULL pointer can be passed.
    ubA                         const real_t* const         Upper constraints' bound vector.
                                                            If no lower constraints' bounds exist, a NULL pointer can be passed.
    nWSR                        int_t&                      Input: Maximum number of working set recalculations when using initial homotopy.
                                                            Output: Number of performed working set recalculations.
    cputime = 0                 real_t* const               Input: Maximum CPU time allowed for QP initialisation.
                                                            Output: CPU time spent for QP initialisation (if pointer passed).
    xOpt = 0                    const real_t* const         Optimal primal solution vector.
                                                            (If a null pointer is passed, the old primal solution is kept!)
    yOpt = 0                    const real_t* const         Optimal dual solution vector.
                                                            (If a null pointer is passed, the old dual solution is kept!)
    guessedBounds = 0           const Bounds* const         Optimal working set of bounds for solution (xOpt,yOpt).
                                                            (If a null pointer is passed, all bounds are assumed inactive!)
    guessedConstraints = 0      const Constraints* const    Optimal working set of constraints for solution (xOpt,yOpt).
                                                            (If a null pointer is passed, all constraints are assumed inactive!)
    _R = 0                      const real_t* const         Pre-computed (upper triangular) Cholesky factor of Hessian matrix.
                                                            The Cholesky factor must be stored in a real_t array of size nV*nV
                                                            in row-major format. Note: Only used if xOpt/yOpt and gB are NULL!
                                                            (If a null pointer is passed, Cholesky decomposition is computed internally!)
***/

// [[Rcpp::export]]
SEXP init_qproblem(SEXP r_model, NumericVector r_H, NumericVector r_g, NumericVector r_A, 
                   NumericVector r_lb, NumericVector r_ub, NumericVector r_lbA, 
                   NumericVector r_ubA, int r_nWSRm, double r_cputime) {
                   // NumericVector r_xopt, NumericVector r_yopt, NumericVector r_guessed_bounds,
                   // NumericVector ) {

    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    using namespace qpOASES;

    const real_t* H   = &*r_H.begin();
    const real_t* g   = &*r_g.begin();
    const real_t* A   = &*r_A.begin();
    const real_t* lb  = &*r_lb.begin();
    const real_t* ub  = &*r_ub.begin();
    const real_t* lbA = &*r_lbA.begin();
    const real_t* ubA = &*r_ubA.begin();
    real_t cpu_time = r_cputime;

    int_t nwsrm = r_nWSRm;
    returnValue status_enum = model->init(&H[0], &g[0], &A[0], &lb[0], &ub[0], 
                                          &lbA[0], &ubA[0], nwsrm, &cpu_time);
    int status_int = static_cast<int>(status_enum);

    return Rcpp::wrap(status_int);
}

// [[Rcpp::export]]
SEXP init_qproblemb(SEXP r_model, NumericVector r_H, NumericVector r_g, 
                    NumericVector r_lb, NumericVector r_ub, int r_nWSRm, double r_cputime) {
                    // NumericVector r_xopt, NumericVector r_yopt, NumericVector r_guessed_bounds,
                    // NumericVector ) {

    Rcpp::XPtr< qpOASES::QProblemB >model(r_model);
    using namespace qpOASES;

    const real_t* H   = &*r_H.begin();
    const real_t* g   = &*r_g.begin();
    const real_t* lb  = &*r_lb.begin();
    const real_t* ub  = &*r_ub.begin();

    int_t nwsrm = r_nWSRm;
    returnValue status_enum = model->init(&H[0], &g[0], &lb[0], &ub[0], nwsrm);
    int status_int = static_cast<int>(status_enum);

    return Rcpp::wrap(status_int);
}

// [[Rcpp::export]]
SEXP init_sqproblem(SEXP r_model, NumericVector r_H, NumericVector r_g, NumericVector r_A, 
                    NumericVector r_lb, NumericVector r_ub, NumericVector r_lbA, 
                    NumericVector r_ubA, int r_nWSRm, double r_cputime) {
                    // NumericVector r_xopt, NumericVector r_yopt, NumericVector r_guessed_bounds,
                    // NumericVector ) {

    Rcpp::XPtr< qpOASES::SQProblem >model(r_model);
    using namespace qpOASES;

    const real_t* H   = &*r_H.begin();
    const real_t* g   = &*r_g.begin();
    const real_t* A   = &*r_A.begin();
    const real_t* lb  = &*r_lb.begin();
    const real_t* ub  = &*r_ub.begin();
    const real_t* lbA = &*r_lbA.begin();
    const real_t* ubA = &*r_ubA.begin();

    int_t nwsrm = r_nWSRm;
    returnValue status_enum = model->init(&H[0], &g[0], &A[0], &lb[0], &ub[0], &lbA[0], &ubA[0], nwsrm);
    int status_int = static_cast<int>(status_enum);

    return Rcpp::wrap(status_int);
}

// [[Rcpp::export]]
SEXP hotstart_qproblem(SEXP r_model, NumericVector r_g, NumericVector r_lb, NumericVector r_ub, 
                       NumericVector r_lbA, NumericVector r_ubA, int r_nWSR) {

    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    using namespace qpOASES;

    const real_t* g   = &*r_g.begin();
    const real_t* lb  = &*r_lb.begin();
    const real_t* ub  = &*r_ub.begin();
    const real_t* lbA = &*r_lbA.begin();
    const real_t* ubA = &*r_ubA.begin();

    int_t nwsr = r_nWSR;
    model->hotstart(&g[0], &lb[0], &ub[0], &lbA[0], &ubA[0], nwsr);

    return R_NilValue;
}

// [[Rcpp::export]]
SEXP hotstart_qproblemb(SEXP r_model, NumericVector r_g, NumericVector r_lb, NumericVector r_ub, 
                  int r_nWSR) {

    Rcpp::XPtr< qpOASES::QProblemB >model(r_model);
    using namespace qpOASES;

    const real_t* g   = &*r_g.begin();
    const real_t* lb  = &*r_lb.begin();
    const real_t* ub  = &*r_ub.begin();

    int_t nwsr = r_nWSR;
    model->hotstart(&g[0], &lb[0], &ub[0], nwsr);

    return R_NilValue;
}



// [[Rcpp::export]]
SEXP print_options(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    model->printOptions();
    return R_NilValue;   
}

// [[Rcpp::export]]
double get_objval(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    double objval = model->getObjVal();
    return objval;   
}


// [[Rcpp::export]]
int get_number_of_variables(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNV();
}


// [[Rcpp::export]]
int get_number_of_free_variables(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNFR();
}


// [[Rcpp::export]]
int get_number_of_fixed_variables(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNFX();
}


// [[Rcpp::export]]
int get_number_of_constraints(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNC();
}


// [[Rcpp::export]]
int get_number_of_equality_constraints(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNEC();
}


// [[Rcpp::export]]
int get_number_of_active_constraints(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNAC();
}


// [[Rcpp::export]]
int get_number_of_inactive_constraints(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->getNIAC();
}


// [[Rcpp::export]]
int is_initialised(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->isInitialised();
}


// [[Rcpp::export]]
int is_solved(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->isSolved();
}


// [[Rcpp::export]]
int is_infeasible(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->isInfeasible();
}


// [[Rcpp::export]]
int is_unbounded(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);
    return (int) model->isUnbounded();
}


// [[Rcpp::export]]
SEXP get_primal_solution(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);

    int n = model->getNV();
    std::vector<qpOASES::real_t> vec(n);
    model->getPrimalSolution( &vec[0] );

    // qpOASES::real_t xopt[n];
    // model->getPrimalSolution( xopt );
    // NumericVector vec(n);
    // for ( int i = 0 ; i < n; i++ ) {
    //     vec[i] = xopt[i];
    // }
    return Rcpp::wrap(vec);
}


// [[Rcpp::export]]
SEXP get_dual_solution(SEXP r_model) {
    Rcpp::XPtr< qpOASES::QProblem >model(r_model);

    int n = model->getNV() + model->getNC();
    std::vector<qpOASES::real_t> vec(n);
    model->getDualSolution( &vec[0] );

    // qpOASES::real_t yopt[n];
    // model->getDualSolution( yopt );
    // NumericVector vec(n);
    // for ( int i = 0 ; i < n; i++ ) {
    //     vec[i] = yopt[i];
    // }
    return Rcpp::wrap(vec);
}


/*
 *  Read Test Problems from the QP Benchmark Collection
 */

// readOQPdimensions
// [[Rcpp::export]]
SEXP read_oqp_dimensions(std::string r_path) {
    using namespace qpOASES;
    int_t nQP, nV, nC, nEC;
    const char* path = r_path.c_str();
    readOqpDimensions(path, nQP, nV, nC, nEC);
    return List::create(Named("number_of_qps") = nQP, 
        Named("number_of_varibales") = nV, Named("number_of_constraints") = nC, 
        Named("number_of_equality_constraints") = nEC);
}

// readOQPdata
SEXP read_oqp_data(std::string r_path) {
    using namespace qpOASES;
    int_t nQP, nV, nC, nEC;
    const char* path = r_path.c_str();
    // TODO: !
    return R_NilValue;
}


// returnValue readOqpData(    const char* path,   /**< Full path of the data files (without trailing slash!). */
//                             int_t& nQP,         /**< Output: Number of QPs. */
//                             int_t& nV,          /**< Output: Number of variables. */
//                             int_t& nC,          /**< Output: Number of constraints. */
//                             int_t& nEC,         /**< Output: Number of equality constraints. */
//                             real_t** H,         /**< Output: Hessian matrix. */
//                             real_t** g,         /**< Output: Sequence of gradient vectors. */
//                             real_t** A,         /**< Output: Constraint matrix. */
//                             real_t** lb,        /**< Output: Sequence of lower bound vectors (on variables). */
//                             real_t** ub,        /**< Output: Sequence of upper bound vectors (on variables). */
//                             real_t** lbA,       /**< Output: Sequence of lower constraints' bound vectors. */
//                             real_t** ubA,       /**< Output: Sequence of upper constraints' bound vectors. */
//                             real_t** xOpt,      /**< Output: Sequence of primal solution vectors
//                                                  *           (not read if a null pointer is passed). */
//                             real_t** yOpt,      /**< Output: Sequence of dual solution vectors
//                                                  *           (not read if a null pointer is passed). */
//                             real_t** objOpt     /**< Output: Sequence of optimal objective function values
//                                                  *           (not read if a null pointer is passed). */
//                             );


// solveOQPbenchmark
// runOQPbenchmark



/*

Page 19

SQProblem
  - H and A varry in the differnt problems therefore they can be provided in
    the hotstart method.

QProblemB
  - This problems have no 

*/






