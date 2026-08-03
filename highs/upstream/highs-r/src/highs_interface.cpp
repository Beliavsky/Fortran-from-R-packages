
#include "highs_interface.h"


using namespace Rcpp;


/*
class HighsR: public Highs {
    public:

    void set_ncol(int32_t ncol) {
        model_.lp_.num_col_ = ncol;
    }

};
*/

static void R_message_handler(HighsLogType type, const char* message, void* log_callback_data) {
    // Only print if log_callback_data is valid and output_flag is true
    if (log_callback_data) {
        Highs* highs = static_cast<Highs*>(log_callback_data);
        const HighsOptions& options = highs->getOptions();
        if (options.log_options.output_flag &&
            *options.log_options.output_flag) {
            Rcpp::Rcout << message << std::endl;
        }
    }
}


std::vector<HighsVarType> to_vartype(std::vector<int32_t> type) {
    std::vector<HighsVarType> htype;
    std::vector<HighsVarType> variable_types = {
        HighsVarType::kContinuous, HighsVarType::kInteger, HighsVarType::kSemiContinuous,
        HighsVarType::kSemiInteger, HighsVarType::kImplicitInteger};
    for (std::size_t i = 0; i < type.size(); ++i) {
        htype.push_back(variable_types[type[i]]);
    }
    return htype;
}


// [[Rcpp::export]]
SEXP new_model() {
    Rcpp::XPtr<HighsModel> highs_model(new HighsModel(), true);
    return highs_model;
}

// [[Rcpp::export]]
SEXP model_set_ncol(SEXP mpt, int32_t ncol) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_col_ = ncol;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_nrow(SEXP mpt, int32_t nrow) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_row_ = nrow;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_sense(SEXP mpt, bool maximum) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (maximum) {
        model->lp_.sense_ = ObjSense::kMaximize;
    } else {
        model->lp_.sense_ = ObjSense::kMinimize;
    }
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_offset(SEXP mpt, double_t offset) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.offset_ = offset;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_objective(SEXP mpt, std::vector<double> objective) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_cost_ = objective;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lower(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_upper(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_constraint_matrix_(SEXP mpt, std::string format,
                                  std::vector<int32_t> start,
                                  std::vector<int32_t> index,
                                  std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (format == "colwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    } else if (format == "rowwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    } else if (format == "rowwise_partitioned") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwisePartitioned;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->lp_.a_matrix_.start_ = start;
    model->lp_.a_matrix_.index_ = index;
    model->lp_.a_matrix_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lhs(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_rhs(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_hessian_(SEXP mpt, std::string format, int32_t dim,
                        std::vector<int32_t> start,
                        std::vector<int32_t> index,
                        std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->hessian_.dim_ = dim;
    if (format == "triangular") {
        model->hessian_.format_ = HessianFormat::kTriangular;
    } else if (format == "square") {
        model->hessian_.format_ = HessianFormat::kSquare;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->hessian_.start_ = start;
    model->hessian_.index_ = index;
    model->hessian_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_vartype(SEXP mpt, std::vector<int32_t> type) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (model->lp_.integrality_.size() < type.size()) {
        model->lp_.integrality_.resize(type.size());
    }
    std::vector<HighsVarType> variable_types = {
        HighsVarType::kContinuous, HighsVarType::kInteger, HighsVarType::kSemiContinuous,
        HighsVarType::kSemiInteger, HighsVarType::kImplicitInteger};
    for (std::size_t i = 0; i < type.size(); ++i) {
        model->lp_.integrality_[i] = variable_types[type[i]];
    }
    return R_NilValue;
}

// [[Rcpp::export]]
int32_t model_get_nvars(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_col_);
}

// [[Rcpp::export]]
int32_t model_get_ncons(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_row_);
}


RCPP_MODULE(RcppHighs) {
    class_<Highs>("Highs")

    .constructor()

    
    .method("getObjectiveValue", &Highs::getObjectiveValue)
    .method("getNumCol", &Highs::getNumCol)
    .method("getNumRow", &Highs::getNumRow)
    .method("getNumNz", &Highs::getNumNz)
    .method("getHessianNumNz", &Highs::getHessianNumNz)
    // .method("getObjectiveSense", &Highs::getObjectiveSense)
    // .method("getObjectiveOffset", &Highs::getObjectiveOffset)

    ;
}


// [[Rcpp::export]]
SEXP new_solver(SEXP mpt) {
    Rcpp::XPtr<Highs> highs(new Highs(), true);
    highs->setLogCallback(R_message_handler, highs.get());

    if (Rf_isNull(mpt)) {
        return Rcpp::List::create(
            Named("status") = 0,
            Named("solver") = highs
        );
    }
    Rcpp::XPtr<HighsModel>model(mpt);
    HighsStatus return_status = highs->passModel(*model.get());
    if (return_status == HighsStatus::kError) {
        return Rcpp::List::create(
            Named("status") = static_cast<int32_t>(return_status),
            Named("solver") = R_NilValue
        );
    }
    return Rcpp::List::create(
        Named("status") = static_cast<int32_t>(return_status),
        Named("solver") = highs
    );
}


// [[Rcpp::export]]
SEXP highs_pass_model(SEXP hi,
                   const int32_t num_col, const int32_t num_row,
                   const int32_t num_nz, const int32_t a_format,
                   const int32_t sense, const double_t offset,
                   NumericVector col_cost, NumericVector col_lower,
                   NumericVector col_upper, NumericVector row_lower,
                   NumericVector row_upper, IntegerVector a_start,
                   IntegerVector a_index, NumericVector a_value,
                   IntegerVector integrality) {
    Rcpp::XPtr<Highs>highs(hi);
    highs->passModel(num_col,
                     num_row,
                     num_nz,
                     a_format,
                     sense,
                     offset,
                     &(col_cost[0]),
                     &(col_lower[0]),
                     &(col_upper[0]),
                     &(row_lower[0]),
                     &(row_upper[0]),
                     &(a_start[0]),
                     &(a_index[0]),
                     &(a_value[0]),
                     &(integrality[0]));
    return R_NilValue;
}


// [[Rcpp::export]]
SEXP solver_pass_hessian() {
    return R_NilValue;
}


// [[Rcpp::export]]
SEXP solver_pass_constraints () {
    return R_NilValue;
}


// [[Rcpp::export]]
int32_t solver_get_sense(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    ObjSense sense;
    HighsStatus status = highs->getObjectiveSense(sense);
    if (status == HighsStatus::kOk) {
        return sense == ObjSense::kMaximize;
    } else {
        Rcpp::stop("could not obtain the sense of the objective function.");
    }
}

// [[Rcpp::export]]
int32_t solver_set_sense(SEXP hi, bool maximum) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    if (maximum) {
        status = highs->changeObjectiveSense(ObjSense::kMaximize);
    } else {
        status = highs->changeObjectiveSense(ObjSense::kMinimize);
    }
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_offset(SEXP hi, const double ext_offset) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeObjectiveOffset(ext_offset);
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_integrality(SEXP hi, std::vector<int32_t> index, std::vector<int32_t> type) {
    Rcpp::XPtr<Highs>highs(hi);
    const int32_t num_set_entries = index.size();
    const HighsInt* set = static_cast<HighsInt*>(index.data());
    std::vector<HighsVarType> var_types = to_vartype(type);
    const HighsVarType* integrality = var_types.data();
    HighsStatus status = highs->changeColsIntegrality(num_set_entries, set, integrality);
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_objective(SEXP hi, std::vector<int32_t> index, std::vector<double_t> obj) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeColsCost(index.size(), &(index[0]), &(obj[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_set_variable_bounds(SEXP hi, std::vector<int32_t> index, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeColsBounds(index.size(), &(index[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_set_constraint_bounds(SEXP hi, std::vector<int32_t> index, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeRowsBounds(index.size(), &(index[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
SEXP solver_set_coeff(SEXP hi, std::vector<int32_t> row, std::vector<int32_t> col, std::vector<double_t> val) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    for (std::size_t i = 0; i < row.size(); ++i) {
        status = highs->changeCoeff(row[i], col[i], val[i]);
        if (status != HighsStatus::kOk) {
            Rcpp::stop("could not change the coefficient.");
        }
    }
    return R_NilValue;
}


// [[Rcpp::export]]
int32_t solver_get_objective_sense(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    ObjSense sense;
    HighsStatus status = highs->getObjectiveSense(sense);
    if (status == HighsStatus::kOk) {
        return static_cast<int32_t>(sense);
    } else {
        Rcpp::stop("could not obtain the sense of the objective function.");
    }
}


// [[Rcpp::export]]
double_t solver_get_objective_offset(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    double_t offset;
    HighsStatus status = highs->getObjectiveOffset(offset);
    if (status == HighsStatus::kOk) {
        return offset;
    } else {
        Rcpp::stop("could not obtain the offset of the objective function.");
    }
}


// [[Rcpp::export]]
SEXP solver_get_cols(SEXP hi, std::vector<int32_t> col_indices) {
    Rcpp::XPtr<Highs>highs(hi);
    if (col_indices.size() == 0) {
        return List::create();
    }

    HighsInt* indices_ptr = static_cast<HighsInt*>(col_indices.data());

    std::vector<double> cost(col_indices.size());
    std::vector<double> lower(col_indices.size());
    std::vector<double> upper(col_indices.size());

    double* cost_ptr = static_cast<double*>(cost.data());
    double* lower_ptr = static_cast<double*>(lower.data());
    double* upper_ptr = static_cast<double*>(upper.data());
    HighsInt get_num_col;
    HighsInt get_num_nz;
    HighsStatus status = highs->getCols(col_indices.size(), indices_ptr, get_num_col, cost_ptr,
                                        lower_ptr, upper_ptr, get_num_nz, nullptr, nullptr, nullptr);
    if (status != HighsStatus::kOk) {
        Rcpp::stop("could not obtain the columns.");
    }
    Rcpp::List ret = Rcpp::List::create(
        Named("status") = static_cast<int32_t>(status),
        Named("num_col") = get_num_col,
        Named("num_nz") = get_num_nz,
        Named("cost") = cost,
        Named("lower") = lower,
        Named("upper") = upper
    );
    return ret;
}


// [[Rcpp::export]]
int32_t solver_get_col_integrality(SEXP hi, int32_t col_index) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsVarType integrality;
    HighsStatus status = highs->getColIntegrality(col_index, integrality);
    if (status != HighsStatus::kOk) {
        Rcpp::stop("could not obtain the integrality of the column.");
    }
    return static_cast<int32_t>(integrality);
}


// [[Rcpp::export]]
SEXP solver_get_rows(SEXP hi, std::vector<int32_t> row_indices) {
    Rcpp::XPtr<Highs>highs(hi);
    if (row_indices.size() == 0) {
        return List::create();
    }

    HighsInt* indices_ptr = static_cast<HighsInt*>(row_indices.data());
    std::vector<double> lower(row_indices.size());
    std::vector<double> upper(row_indices.size());
    double* lower_ptr = static_cast<double*>(lower.data());
    double* upper_ptr = static_cast<double*>(upper.data());
    HighsInt get_num_row;
    HighsInt get_num_nz;
    HighsStatus status = highs->getRows(row_indices.size(), indices_ptr, get_num_row, lower_ptr, upper_ptr, get_num_nz, nullptr, nullptr, nullptr);
    if (status != HighsStatus::kOk) {
        Rcpp::stop("could not obtain the rows.");
    }
    Rcpp::List ret = Rcpp::List::create(
        Named("status") = static_cast<int32_t>(status),
        Named("num_row") = get_num_row,
        Named("num_nz") = get_num_nz,
        Named("lower") = lower,
        Named("upper") = upper
    );
    return ret;
}


// [[Rcpp::export]]
SEXP solver_get_coeff(SEXP hi,  int32_t row_index, int32_t col_index) {
    Rcpp::XPtr<Highs>highs(hi);
    double val;
    HighsStatus status = highs->getCoeff(row_index, col_index, val);
    if (status != HighsStatus::kOk) {
        Rcpp::stop("could not obtain the coefficient.");
    }
    return Rcpp::wrap(val);
}


// [[Rcpp::export]]
int32_t solver_delete_cols(SEXP hi, std::vector<int32_t> col_indices) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->deleteCols(col_indices.size(), &(col_indices[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_delete_rows(SEXP hi, std::vector<int32_t> row_indices) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->deleteRows(row_indices.size(), &(row_indices[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_scale_col(SEXP hi, int32_t col_index, double_t scale) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->scaleCol(col_index, scale);
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_scale_row(SEXP hi, int32_t row_index, double_t scale) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->scaleRow(row_index, scale);
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_postsolve(SEXP hi, SEXP solution, SEXP basis) {
    Rcpp::stop("not implemented.");
}


// [[Rcpp::export]]
int32_t solver_mip_postsolve(SEXP hi, SEXP solution) {
    Rcpp::stop("not implemented.");
}


// TODO
// HighsStatus postsolve(const HighsSolution& solution, const HighsBasis& basis);




// [[Rcpp::export]]
int32_t solver_add_vars(SEXP hi, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->addVars(lower.size(), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}

/*
model_set_lower
model_set_upper
model_set_constraint_matrix_
model_set_lhs
model_set_rhs
model_set_hessian_
model_set_vartype
*/
//
// changeObjectiveSenseInterface
// 
// addRowsInterface
// addColsInterface
// setSolution
// getHotStart

// [[Rcpp::export]]
int32_t solver_set_option(SEXP hi, std::string key, SEXP value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    if (Rf_isLogical(value)) {
        bool logval = Rcpp::as<bool>(value);
        status = highs->setOptionValue(key, logval);
    } else if (Rf_isInteger(value)) {
        HighsInt intval = Rcpp::as<int32_t>(value);
        status = highs->setOptionValue(key, intval);
    } else if (Rf_isNumeric(value)) {
        double numval = Rcpp::as<double>(value);
        status = highs->setOptionValue(key, numval);
    } else if (Rf_isString(value)) {
        std::string strval = Rcpp::as<std::string>(value);
        status = highs->setOptionValue(key, strval);
    } else {
        Rcpp::stop("unkown type of value.");
    }
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_clear(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clear();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_model(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearModel();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_solver(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearSolver();
    return static_cast<int32_t>(return_status);
}



// [[Rcpp::export]]
int32_t solver_run(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->run();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
SEXP solver_get_model(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel mpt = highs->getModel();
    Rcpp::XPtr<HighsModel> highs_model(&mpt, true);
    return highs_model;
}


// [[Rcpp::export]]
int32_t solver_get_num_col(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getNumCol();
}


// [[Rcpp::export]]
int32_t solver_get_num_row(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getNumRow();
}



//
// int32_t solver_presolve(SEXP hi) {
//     Rcpp::XPtr<Highs>highs(hi);
//     HighsStatus return_status = highs->presolve();
//     return static_cast<int32_t>(return_status);
// }


// TODO
// const HighsLp& getPresolvedLp()


// TODO
// const HighsModel& getPresolvedModel()


// [[Rcpp::export]]
int32_t solver_write_model(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeModel(filename);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_write_basis(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeBasis(filename);
    return static_cast<int32_t>(return_status);
}


// style int in {-1, 0, 1, 2, 3, 4}
// [[Rcpp::export]]
int32_t solver_write_solution(SEXP hi, const std::string filename, const int32_t style) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeSolution(filename, style);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
std::string solver_status_message(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return highs->modelStatusToString(model_status);
}


// [[Rcpp::export]]
int32_t solver_status(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return static_cast<int32_t>(model_status);
}


// [[Rcpp::export]]
double_t  solver_infinity() {
    Highs highs;
    double_t dbl_inf = highs.getInfinity();
    return dbl_inf;
}


// [[Rcpp::export]]
SEXP reset_global_scheduler(bool blocking) {
    Highs highs;
    highs.resetGlobalScheduler(blocking);
    return R_NilValue;
}


// [[Rcpp::export]]
Rcpp::List solver_info(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsInfo& info = highs->getInfo();
    List z = List::create(
        Named("valid") = info.valid,
        Named("mip_node_count") = info.mip_node_count,
        Named("simplex_iteration_count") = info.simplex_iteration_count,
        Named("ipm_iteration_count") = info.ipm_iteration_count,
        Named("qp_iteration_count") = info.qp_iteration_count,
        Named("crossover_iteration_count") = info.crossover_iteration_count,
        Named("primal_solution_status") = highs->solutionStatusToString(info.primal_solution_status),
        Named("dual_solution_status") = highs->solutionStatusToString(info.dual_solution_status),
        Named("basis_validity") = info.basis_validity,
        Named("objective_function_value") = info.objective_function_value,
        Named("mip_dual_bound") = info.mip_dual_bound,
        Named("mip_gap") = info.mip_gap,
        Named("num_primal_infeasibilities") = info.num_primal_infeasibilities,
        Named("max_primal_infeasibility") = info.max_primal_infeasibility,
        Named("sum_primal_infeasibilities") = info.sum_primal_infeasibilities,
        Named("num_dual_infeasibilities") = info.num_dual_infeasibilities,
        Named("max_dual_infeasibility") = info.max_dual_infeasibility,
        Named("sum_dual_infeasibilities") = info.sum_dual_infeasibilities
    );
    return z;
}


// [[Rcpp::export]]
bool solver_get_bool_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    bool value;
    highs->getOptionValue(key, value);
    // HighsInt    
    return value;
}

// [[Rcpp::export]]
int32_t solver_get_int_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsInt value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
double_t solver_get_dbl_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    double_t value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
std::string solver_get_str_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    std::string value;
    highs->getOptionValue(key, value);   
    return value;
}


// [[Rcpp::export]]
int32_t solver_change_variable_bounds(SEXP hi, IntegerVector idx, NumericVector lower, NumericVector upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeColsBounds(idx.size(), &(idx[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_change_constraint_bounds(SEXP hi, IntegerVector idx, NumericVector lhs, NumericVector rhs) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeRowsBounds(idx.size(), &(idx[0]), &(lhs[0]), &(rhs[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_rows(SEXP hi, NumericVector lhs, NumericVector rhs,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addRows(
        static_cast<int32_t>(lhs.size()), //!< num_new_row = Number of new rows
        &(lhs[0]),  //!< lower = Array of size num_new_row with lower bounds
        &(rhs[0]),  //!< upper = Array of size num_new_row with upper bounds
        static_cast<int32_t>(value.size()), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the rows
        &(index[0]),   //!< indices = Array of size num_new_nz with column indices for all rows
        &(value[0])    //!< values = Array of size num_new_nz with column values for all rows
    );
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_cols(SEXP hi, NumericVector costs,
    NumericVector lower, NumericVector upper,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addCols(
        static_cast<int32_t>(lower.size()), //!< num_new_col = Number of new columns
        &(costs[0]),   //!< costs = Array of size num_new_col with costs
        &(lower[0]),   //!< lower = Array of size num_new_col with lower bounds
        &(upper[0]),   //!< upper = Array of size num_new_col with upper bounds
        static_cast<int32_t>(value.size()), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the columns
        &(index[0]),   //!< indices = Array of size num_new_nz with row indices for all columns
        &(value[0])    //!< values = Array of size num_new_nz with row values for all columns
    );
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_lp_costs(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumCol();
    NumericVector costs(nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        costs[i] = model.lp_.col_cost_[i];
    }
    return costs;
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_variable_bounds(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumCol();
    NumericVector vbounds(2 * nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        vbounds[i] = model.lp_.col_lower_[i];
        vbounds[nvar + i] = model.lp_.col_upper_[i];
    }
    return vbounds;
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_constraint_bounds(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumRow();
    NumericVector lhs_rhs(2 * nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        lhs_rhs[i] = model.lp_.row_lower_[i];
        lhs_rhs[nvar + i] = model.lp_.row_upper_[i];
    }
    return lhs_rhs;
}


// [[Rcpp::export]]
std::vector<double> solver_get_col_cost(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    return model.lp_.col_cost_;
}


// [[Rcpp::export]]
std::vector<double> solver_get_col_lower(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    return model.lp_.col_lower_;
}


// [[Rcpp::export]]
std::vector<double> solver_get_col_upper(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    return model.lp_.col_upper_;
}


// [[Rcpp::export]]
std::vector<double> solver_get_row_lower(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    return model.lp_.row_lower_;
}


// [[Rcpp::export]]
std::vector<double> solver_get_row_upper(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    return model.lp_.row_upper_;
}


// [[Rcpp::export]]
IntegerVector solver_get_integrality(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    int len = static_cast<int>(model.lp_.integrality_.size());
    IntegerVector vec(len);
    for (int32_t i = 0; i < len; i++) {
        vec[i] = static_cast<int>(model.lp_.integrality_[i]);
    }
    return vec;
}


// [[Rcpp::export]]
Rcpp::List solver_get_constraint_matrix(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    HighsSparseMatrix A = model.lp_.a_matrix_;
    List z = List::create(
        Named("format") = static_cast<int32_t>(A.format_),
        Named("nrow") = A.num_row_,
        Named("ncol") = A.num_col_,
        Named("start") = A.start_,
        Named("p_end") = A.p_end_,
        Named("index") = A.index_,
        Named("value") = A.value_
    );
    return z;
}


// [[Rcpp::export]]
Rcpp::IntegerVector solver_get_vartype(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsLp& lp = highs->getLp();
    int32_t m = static_cast<int32_t>(lp.integrality_.size());
    IntegerVector type(m);
    for (int32_t i = 0; i < type.size(); ++i) {
        type[i] = static_cast<int32_t>(lp.integrality_[i]);
    }
    return type;
}


// [[Rcpp::export]]
Rcpp::List solver_get_solution(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsSolution& solution = highs->getSolution();
    NumericVector col_value( solution.col_value.begin(), solution.col_value.end() );
    NumericVector col_dual( solution.col_dual.begin(), solution.col_dual.end() );
    NumericVector row_value( solution.row_value.begin(), solution.row_value.end() );
    NumericVector row_dual( solution.row_dual.begin(), solution.row_dual.end() );
    List z = List::create(
        Named("value_valid") = static_cast<bool>(solution.value_valid),
        Named("dual_valid") = static_cast<bool>(solution.dual_valid),
        Named("col_value") = col_value,
        Named("col_dual") = col_dual,
        Named("row_value") = row_value,
        Named("row_dual") = row_dual
    );
    return z;
}


Rcpp::List solver_get_solution__deprecated(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsSolution solution = highs->getSolution();
    List z = List::create(
        Named("value_valid") = static_cast<bool>(solution.value_valid),
        Named("dual_valid") = static_cast<bool>(solution.dual_valid),
        Named("col_value") = solution.col_value,
        Named("col_dual") = solution.col_dual,
        Named("row_value") = solution.row_value,
        Named("row_dual") = solution.row_dual
    );
    return z;
}


// [[Rcpp::export]]
int32_t solver_set_solution_obj(SEXP hi,
                                bool value_valid,
                                bool dual_valid,
                                std::vector<double> col_value,
                                std::vector<double> col_dual,
                                std::vector<double> row_value,
                                std::vector<double> row_dual) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsSolution solution = {value_valid, dual_valid, col_value, col_dual, row_value, row_dual};
    HighsStatus return_status = highs->setSolution(solution);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_set_solution_vec(SEXP hi, IntegerVector idx, NumericVector val) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->setSolution(static_cast<int32_t>(idx.size()), &(idx[0]), &(val[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
Rcpp::List solver_get_dual_ray(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    bool has_dual_ray = false;
    HighsInt num_row = highs->getNumRow();
    std::vector<double> dual_ray_value(num_row);
    HighsStatus return_status =
        highs->getDualRay(has_dual_ray, num_row > 0 ? dual_ray_value.data() : nullptr);

    Rcpp::List ret = Rcpp::List::create(
        Named("status") = static_cast<int32_t>(return_status),
        Named("has_dual_ray") = has_dual_ray
    );
    if (has_dual_ray && num_row > 0) {
        ret["dual_ray"] = Rcpp::NumericVector(dual_ray_value.begin(), dual_ray_value.end());
    } else {
        ret["dual_ray"] = R_NilValue;
    }
    return ret;
}


// [[Rcpp::export]]
Rcpp::List solver_get_primal_ray(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    bool has_primal_ray = false;
    HighsInt num_col = highs->getNumCol();
    std::vector<double> primal_ray_value(num_col);
    HighsStatus return_status =
        highs->getPrimalRay(has_primal_ray, num_col > 0 ? primal_ray_value.data() : nullptr);

    Rcpp::List ret = Rcpp::List::create(
        Named("status") = static_cast<int32_t>(return_status),
        Named("has_primal_ray") = has_primal_ray
    );
    if (has_primal_ray && num_col > 0) {
        ret["primal_ray"] = Rcpp::NumericVector(primal_ray_value.begin(), primal_ray_value.end());
    } else {
        ret["primal_ray"] = R_NilValue;
    }
    return ret;
}

// solver_set_hotstart
// HighsStatus return_status = highs->setHotStart(hotstart);


// [[Rcpp::export]]
Rcpp::List solver_get_basis(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsBasis& basis = highs->getBasis();
    int ncol = highs->getNumCol();
    int nrow = highs->getNumRow();
    IntegerVector col_status(ncol);
    IntegerVector row_status(nrow);
    for (int i = 0; i < ncol; i++) {
        col_status[i] = static_cast<int>(basis.col_status[i]);
    }
    for (int i = 0; i < nrow; i++) {
        row_status[i] = static_cast<int>(basis.row_status[i]);
    }
    List z = List::create(
        Named("valid") = basis.valid,
        Named("col_status") = col_status,
        Named("row_status") = row_status
    );
    return z;
}


// [[Rcpp::export]]
int32_t solver_set_basis(SEXP hi, IntegerVector col_status, IntegerVector row_status) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsBasis basis;
    basis.valid = true;
    basis.alien = true;
    basis.col_status.resize(col_status.size());
    for (int i = 0; i < col_status.size(); i++) {
        if (col_status[i] < 0 || col_status[i] > 4) {
            Rcpp::stop("Invalid column basis status %d at index %d (must be 0-4)",
                        col_status[i], i);
        }
        basis.col_status[i] = static_cast<HighsBasisStatus>(col_status[i]);
    }
    basis.row_status.resize(row_status.size());
    for (int i = 0; i < row_status.size(); i++) {
        if (row_status[i] < 0 || row_status[i] > 4) {
            Rcpp::stop("Invalid row basis status %d at index %d (must be 0-4)",
                        row_status[i], i);
        }
        basis.row_status[i] = static_cast<HighsBasisStatus>(row_status[i]);
    }
    HighsStatus return_status = highs->setBasis(basis, "R_highs_package");
    if (return_status == HighsStatus::kError) {
        Rcpp::stop("HiGHS setBasis failed");
    }
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_basis(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->setBasis();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
Rcpp::List solver_get_ranging(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsRanging ranging;
    HighsStatus status = highs->getRanging(ranging);
    if (status == HighsStatus::kError) {
        Rcpp::stop("HiGHS getRanging failed");
    }
    auto pack_record = [](const HighsRangingRecord& rec) {
        return List::create(
            Named("value") = rec.value_,
            Named("objective") = rec.objective_,
            Named("in_var") = rec.in_var_,
            Named("ou_var") = rec.ou_var_
        );
    };
    List z = List::create(
        Named("valid") = ranging.valid,
        Named("col_cost_up") = pack_record(ranging.col_cost_up),
        Named("col_cost_dn") = pack_record(ranging.col_cost_dn),
        Named("col_bound_up") = pack_record(ranging.col_bound_up),
        Named("col_bound_dn") = pack_record(ranging.col_bound_dn),
        Named("row_bound_up") = pack_record(ranging.row_bound_up),
        Named("row_bound_dn") = pack_record(ranging.row_bound_dn)
    );
    return z;
}


// [[Rcpp::export]]
std::string solver_version(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->version();
}


// [[Rcpp::export]]
double solver_get_run_time(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getRunTime();
}


// [[Rcpp::export]]
int32_t solver_presolve(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->presolve();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_postsolve_solution(SEXP hi,
                                   NumericVector col_value,
                                   NumericVector col_dual,
                                   NumericVector row_value,
                                   NumericVector row_dual) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = Rcpp::as<std::vector<double>>(col_value);
    solution.col_dual = Rcpp::as<std::vector<double>>(col_dual);
    solution.row_value = Rcpp::as<std::vector<double>>(row_value);
    solution.row_dual = Rcpp::as<std::vector<double>>(row_dual);
    HighsStatus return_status = highs->postsolve(solution);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_read_model(SEXP hi, std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->readModel(filename);
    if (return_status == HighsStatus::kError) {
        Rcpp::stop("HiGHS readModel failed for file: %s", filename.c_str());
    }
    return static_cast<int32_t>(return_status);
}
