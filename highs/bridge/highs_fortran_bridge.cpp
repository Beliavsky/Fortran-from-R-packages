// SPDX-License-Identifier: GPL-2.0-or-later
#include "Highs.h"
#include "highs_fortran_abi.h"
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#ifdef _WIN32
#define HF_EXPORT __declspec(dllexport)
#else
#define HF_EXPORT __attribute__((visibility("default")))
#endif

static int copy_string(const std::string& s, char* buffer, int n) {
  if (buffer && n > 0) {
    const int m = std::min<int>(n - 1, static_cast<int>(s.size()));
    if (m > 0) std::memcpy(buffer, s.data(), static_cast<size_t>(m));
    buffer[m] = '\0';
  }
  return static_cast<int>(s.size());
}

static Highs* as_highs(void* p) { return static_cast<Highs*>(p); }

extern "C" {

HF_EXPORT void* hf_backend_create(void) {
  try { return new Highs(); } catch (...) { return nullptr; }
}

HF_EXPORT void hf_backend_destroy(void* p) { delete as_highs(p); }

HF_EXPORT int hf_backend_pass_model(
    void* p, int num_col, int num_row, int num_nz, int a_format,
    int sense, double offset, const double* col_cost,
    const double* col_lower, const double* col_upper,
    const double* row_lower, const double* row_upper,
    const int* a_start, const int* a_index, const double* a_value,
    const int* integrality) {
  if (!p) return -1;
  return static_cast<int>(as_highs(p)->passModel(
      num_col, num_row, num_nz, a_format, sense, offset,
      col_cost, col_lower, col_upper, row_lower, row_upper,
      a_start, a_index, a_value, integrality));
}

HF_EXPORT int hf_backend_pass_hessian(void* p, int dim, int num_nz,
                                      int format, const int* start,
                                      const int* index, const double* value) {
  if (!p) return -1;
  return static_cast<int>(as_highs(p)->passHessian(dim, num_nz, format,
                                                   start, index, value));
}

HF_EXPORT int hf_backend_run(void* p) {
  return p ? static_cast<int>(as_highs(p)->run()) : -1;
}

HF_EXPORT int hf_backend_model_status(void* p) {
  return p ? static_cast<int>(as_highs(p)->getModelStatus()) : 0;
}

HF_EXPORT double hf_backend_objective_value(void* p) {
  return p ? as_highs(p)->getObjectiveValue() : 0.0;
}

HF_EXPORT int hf_backend_get_solution(void* p, double* col_value,
                                      double* col_dual, double* row_value,
                                      double* row_dual, int* value_valid,
                                      int* dual_valid) {
  if (!p) return -1;
  const HighsSolution& s = as_highs(p)->getSolution();
  if (value_valid) *value_valid = s.value_valid ? 1 : 0;
  if (dual_valid) *dual_valid = s.dual_valid ? 1 : 0;
  if (col_value && !s.col_value.empty())
    std::copy(s.col_value.begin(), s.col_value.end(), col_value);
  if (col_dual && !s.col_dual.empty())
    std::copy(s.col_dual.begin(), s.col_dual.end(), col_dual);
  if (row_value && !s.row_value.empty())
    std::copy(s.row_value.begin(), s.row_value.end(), row_value);
  if (row_dual && !s.row_dual.empty())
    std::copy(s.row_dual.begin(), s.row_dual.end(), row_dual);
  return 0;
}

HF_EXPORT int hf_backend_get_info(void* p, hf_info_c* out) {
  if (!p || !out) return -1;
  const HighsInfo& i = as_highs(p)->getInfo();
  out->valid = i.valid ? 1 : 0;
  out->mip_node_count = i.mip_node_count;
  out->simplex_iteration_count = i.simplex_iteration_count;
  out->ipm_iteration_count = i.ipm_iteration_count;
  out->crossover_iteration_count = i.crossover_iteration_count;
  out->qp_iteration_count = i.qp_iteration_count;
  out->primal_solution_status = i.primal_solution_status;
  out->dual_solution_status = i.dual_solution_status;
  out->basis_validity = i.basis_validity;
  out->objective_function_value = i.objective_function_value;
  out->mip_dual_bound = i.mip_dual_bound;
  out->mip_gap = i.mip_gap;
  out->max_integrality_violation = i.max_integrality_violation;
  out->num_primal_infeasibilities = i.num_primal_infeasibilities;
  out->max_primal_infeasibility = i.max_primal_infeasibility;
  out->sum_primal_infeasibilities = i.sum_primal_infeasibilities;
  out->num_dual_infeasibilities = i.num_dual_infeasibilities;
  out->max_dual_infeasibility = i.max_dual_infeasibility;
  out->sum_dual_infeasibilities = i.sum_dual_infeasibilities;
  out->run_time = as_highs(p)->getRunTime();
  return 0;
}

HF_EXPORT int hf_backend_status_message(void* p, char* buffer, int n) {
  if (!p) return copy_string("invalid solver", buffer, n);
  Highs* h = as_highs(p);
  return copy_string(h->modelStatusToString(h->getModelStatus()), buffer, n);
}

HF_EXPORT int hf_backend_version(char* buffer, int n) {
  Highs h;
  return copy_string(h.version(), buffer, n);
}

HF_EXPORT double hf_backend_infinity(void* p) {
  Highs local;
  return p ? as_highs(p)->getInfinity() : local.getInfinity();
}

HF_EXPORT int hf_backend_num_col(void* p) { return p ? as_highs(p)->getNumCol() : -1; }
HF_EXPORT int hf_backend_num_row(void* p) { return p ? as_highs(p)->getNumRow() : -1; }

HF_EXPORT int hf_backend_get_sense(void* p) {
  if (!p) return 0;
  ObjSense sense;
  if (as_highs(p)->getObjectiveSense(sense) == HighsStatus::kError) return 0;
  return static_cast<int>(sense);
}

HF_EXPORT int hf_backend_set_bool_option(void* p, const char* name, int value) {
  return p && name ? static_cast<int>(as_highs(p)->setOptionValue(name, value != 0)) : -1;
}
HF_EXPORT int hf_backend_set_int_option(void* p, const char* name, int value) {
  return p && name ? static_cast<int>(as_highs(p)->setOptionValue(name, value)) : -1;
}
HF_EXPORT int hf_backend_set_double_option(void* p, const char* name, double value) {
  return p && name ? static_cast<int>(as_highs(p)->setOptionValue(name, value)) : -1;
}
HF_EXPORT int hf_backend_set_string_option(void* p, const char* name, const char* value) {
  return p && name && value ? static_cast<int>(as_highs(p)->setOptionValue(name, value)) : -1;
}
HF_EXPORT int hf_backend_reset_options(void* p) {
  return p ? static_cast<int>(as_highs(p)->resetOptions()) : -1;
}
HF_EXPORT int hf_backend_clear_model(void* p) {
  return p ? static_cast<int>(as_highs(p)->clearModel()) : -1;
}
HF_EXPORT int hf_backend_clear_solver(void* p) {
  return p ? static_cast<int>(as_highs(p)->clearSolver()) : -1;
}
HF_EXPORT int hf_backend_presolve(void* p) {
  return p ? static_cast<int>(as_highs(p)->presolve()) : -1;
}

HF_EXPORT int hf_backend_change_cost(void* p, int n, const int* index,
                                     const double* value) {
  return p ? static_cast<int>(as_highs(p)->changeColsCost(n, index, value)) : -1;
}
HF_EXPORT int hf_backend_change_col_bounds(void* p, int n, const int* index,
                                           const double* lower,
                                           const double* upper) {
  return p ? static_cast<int>(as_highs(p)->changeColsBounds(n, index, lower, upper)) : -1;
}
HF_EXPORT int hf_backend_change_row_bounds(void* p, int n, const int* index,
                                           const double* lower,
                                           const double* upper) {
  return p ? static_cast<int>(as_highs(p)->changeRowsBounds(n, index, lower, upper)) : -1;
}
HF_EXPORT int hf_backend_change_coeff(void* p, int row, int col, double value) {
  return p ? static_cast<int>(as_highs(p)->changeCoeff(row, col, value)) : -1;
}
HF_EXPORT int hf_backend_change_integrality(void* p, int n, const int* index,
                                            const int* values) {
  if (!p) return -1;
  std::vector<HighsVarType> types(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) types[static_cast<size_t>(i)] = static_cast<HighsVarType>(values[i]);
  return static_cast<int>(as_highs(p)->changeColsIntegrality(n, index, types.data()));
}
HF_EXPORT int hf_backend_change_sense(void* p, int sense) {
  return p ? static_cast<int>(as_highs(p)->changeObjectiveSense(static_cast<ObjSense>(sense))) : -1;
}
HF_EXPORT int hf_backend_change_offset(void* p, double value) {
  return p ? static_cast<int>(as_highs(p)->changeObjectiveOffset(value)) : -1;
}
HF_EXPORT int hf_backend_write_model(void* p, const char* filename) {
  return p && filename ? static_cast<int>(as_highs(p)->writeModel(filename)) : -1;
}
HF_EXPORT int hf_backend_read_model(void* p, const char* filename) {
  return p && filename ? static_cast<int>(as_highs(p)->readModel(filename)) : -1;
}

HF_EXPORT int hf_backend_set_solution(void* p, const double* col_value,
                                      const double* col_dual,
                                      const double* row_value,
                                      const double* row_dual,
                                      int value_valid, int dual_valid) {
  if (!p) return -1;
  Highs* h = as_highs(p);
  HighsSolution s;
  s.value_valid = value_valid != 0;
  s.dual_valid = dual_valid != 0;
  if (col_value) s.col_value.assign(col_value, col_value + h->getNumCol());
  if (col_dual) s.col_dual.assign(col_dual, col_dual + h->getNumCol());
  if (row_value) s.row_value.assign(row_value, row_value + h->getNumRow());
  if (row_dual) s.row_dual.assign(row_dual, row_dual + h->getNumRow());
  return static_cast<int>(h->setSolution(s));
}

HF_EXPORT int hf_backend_get_basis(void* p, int* col_status,
                                   int* row_status, int* valid) {
  if (!p) return -1;
  Highs* h = as_highs(p);
  const HighsBasis& b = h->getBasis();
  if (valid) *valid = b.valid ? 1 : 0;
  for (int i = 0; col_status && i < h->getNumCol() && i < static_cast<int>(b.col_status.size()); ++i)
    col_status[i] = static_cast<int>(b.col_status[static_cast<size_t>(i)]);
  for (int i = 0; row_status && i < h->getNumRow() && i < static_cast<int>(b.row_status.size()); ++i)
    row_status[i] = static_cast<int>(b.row_status[static_cast<size_t>(i)]);
  return 0;
}

HF_EXPORT int hf_backend_set_basis(void* p, const int* col_status,
                                   const int* row_status) {
  if (!p || !col_status || !row_status) return -1;
  Highs* h = as_highs(p);
  HighsBasis b;
  b.valid = true;
  b.alien = true;
  b.col_status.resize(static_cast<size_t>(h->getNumCol()));
  b.row_status.resize(static_cast<size_t>(h->getNumRow()));
  for (int i = 0; i < h->getNumCol(); ++i)
    b.col_status[static_cast<size_t>(i)] = static_cast<HighsBasisStatus>(col_status[i]);
  for (int i = 0; i < h->getNumRow(); ++i)
    b.row_status[static_cast<size_t>(i)] = static_cast<HighsBasisStatus>(row_status[i]);
  return static_cast<int>(h->setBasis(b, "highs-fortran"));
}

HF_EXPORT int hf_backend_clear_basis(void* p) {
  return p ? static_cast<int>(as_highs(p)->setBasis()) : -1;
}

HF_EXPORT int hf_backend_get_dual_ray(void* p, double* value, int* has_ray) {
  if (!p) return -1;
  bool has = false;
  const int status = static_cast<int>(as_highs(p)->getDualRay(has, value));
  if (has_ray) *has_ray = has ? 1 : 0;
  return status;
}
HF_EXPORT int hf_backend_get_primal_ray(void* p, double* value, int* has_ray) {
  if (!p) return -1;
  bool has = false;
  const int status = static_cast<int>(as_highs(p)->getPrimalRay(has, value));
  if (has_ray) *has_ray = has ? 1 : 0;
  return status;
}

}  // extern "C"
