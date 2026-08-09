// This file is part of PIQP-R.
//
// Copyright (c) 2023 piqp-r authors
//
// This source code is licensed under the BSD 2-Clause License found in the
// LICENSE file in the root directory of this source tree.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "piqp_types.h"

#include <piqp/piqp.hpp>

// [[Rcpp::depends(RcppEigen)]]

// Update settings
void piqp_update_settings(piqp::Settings<double>& s, const Rcpp::List rs) {
  if (rs.containsElementNamed("rho_init"))
    s.rho_init = Rcpp::as<double>(rs["rho_init"]);
  if (rs.containsElementNamed("delta_init"))
    s.delta_init = Rcpp::as<double>(rs["delta_init"]);
  if (rs.containsElementNamed("eps_abs"))
    s.eps_abs = Rcpp::as<double>(rs["eps_abs"]);
  if (rs.containsElementNamed("eps_rel"))
    s.eps_rel = Rcpp::as<double>(rs["eps_rel"]);
  if (rs.containsElementNamed("check_duality_gap"))
    s.check_duality_gap = Rcpp::as<bool>(rs["check_duality_gap"]);
  if (rs.containsElementNamed("eps_duality_gap_abs"))
    s.eps_duality_gap_abs = Rcpp::as<double>(rs["eps_duality_gap_abs"]);
  if (rs.containsElementNamed("eps_duality_gap_rel"))
    s.eps_duality_gap_rel = Rcpp::as<double>(rs["eps_duality_gap_rel"]);
  if (rs.containsElementNamed("infeasibility_threshold"))
    s.infeasibility_threshold = Rcpp::as<double>(rs["infeasibility_threshold"]);
  if (rs.containsElementNamed("reg_lower_limit"))
    s.reg_lower_limit = Rcpp::as<double>(rs["reg_lower_limit"]);
  if (rs.containsElementNamed("reg_finetune_lower_limit"))
    s.reg_finetune_lower_limit = Rcpp::as<double>(rs["reg_finetune_lower_limit"]);
  if (rs.containsElementNamed("reg_finetune_primal_update_threshold"))
    s.reg_finetune_primal_update_threshold = Rcpp::as<piqp::isize>(rs["reg_finetune_primal_update_threshold"]);
  if (rs.containsElementNamed("reg_finetune_dual_update_threshold"))
    s.reg_finetune_dual_update_threshold = Rcpp::as<piqp::isize>(rs["reg_finetune_dual_update_threshold"]);
  if (rs.containsElementNamed("max_iter"))
    s.max_iter = Rcpp::as<piqp::isize>(rs["max_iter"]);
  if (rs.containsElementNamed("max_factor_retires"))
    s.max_factor_retires = Rcpp::as<piqp::isize>(rs["max_factor_retires"]);
  if (rs.containsElementNamed("preconditioner_scale_cost"))
    s.preconditioner_scale_cost = Rcpp::as<bool>(rs["preconditioner_scale_cost"]);
  if (rs.containsElementNamed("preconditioner_reuse_on_update"))
    s.preconditioner_reuse_on_update = Rcpp::as<bool>(rs["preconditioner_reuse_on_update"]);
  if (rs.containsElementNamed("preconditioner_iter"))
    s.preconditioner_iter = Rcpp::as<piqp::isize>(rs["preconditioner_iter"]);
  if (rs.containsElementNamed("tau"))
    s.tau = Rcpp::as<double>(rs["tau"]);
  if (rs.containsElementNamed("iterative_refinement_always_enabled"))
    s.iterative_refinement_always_enabled = Rcpp::as<bool>(rs["iterative_refinement_always_enabled"]);
  if (rs.containsElementNamed("iterative_refinement_eps_abs"))
    s.iterative_refinement_eps_abs = Rcpp::as<double>(rs["iterative_refinement_eps_abs"]);
  if (rs.containsElementNamed("iterative_refinement_eps_rel"))
    s.iterative_refinement_eps_rel = Rcpp::as<double>(rs["iterative_refinement_eps_rel"]);
  if (rs.containsElementNamed("iterative_refinement_max_iter"))
    s.iterative_refinement_max_iter = Rcpp::as<piqp::isize>(rs["iterative_refinement_max_iter"]);
  if (rs.containsElementNamed("iterative_refinement_min_improvement_rate"))
    s.iterative_refinement_min_improvement_rate = Rcpp::as<double>(rs["iterative_refinement_min_improvement_rate"]);
  if (rs.containsElementNamed("iterative_refinement_static_regularization_eps"))
    s.iterative_refinement_static_regularization_eps = Rcpp::as<double>(rs["iterative_refinement_static_regularization_eps"]);
  if (rs.containsElementNamed("iterative_refinement_static_regularization_rel"))
    s.iterative_refinement_static_regularization_rel = Rcpp::as<double>(rs["iterative_refinement_static_regularization_rel"]);
  if (rs.containsElementNamed("verbose"))
    s.verbose = Rcpp::as<bool>(rs["verbose"]);
  if (rs.containsElementNamed("compute_timings"))
    s.compute_timings = Rcpp::as<bool>(rs["compute_timings"]);
}



// [[Rcpp::export]]
Rcpp::List get_settings(SEXP solver_p, bool dense_backend) {
  piqp::Settings<double> s = dense_backend ?
    (Rcpp::as<Rcpp::XPtr<piqp::DenseSolver<double>>>(solver_p))->settings() :
    (Rcpp::as<Rcpp::XPtr<piqp::SparseSolver<double>>>(solver_p))->settings();

  Rcpp::List result = Rcpp::List::create();
  result["rho_init"] = s.rho_init;
  result["delta_init"] = s.delta_init;
  result["eps_abs"] = s.eps_abs;
  result["eps_rel"] = s.eps_rel;
  result["check_duality_gap"] = s.check_duality_gap;
  result["eps_duality_gap_abs"] = s.eps_duality_gap_abs;
  result["eps_duality_gap_rel"] = s.eps_duality_gap_rel;
  result["infeasibility_threshold"] = s.infeasibility_threshold;
  result["reg_lower_limit"] = s.reg_lower_limit;
  result["reg_finetune_lower_limit"] = s.reg_finetune_lower_limit;
  result["reg_finetune_primal_update_threshold"] = s.reg_finetune_primal_update_threshold;
  result["reg_finetune_dual_update_threshold"] = s.reg_finetune_dual_update_threshold;
  result["max_iter"] = s.max_iter;
  result["max_factor_retires"] = s.max_factor_retires;
  result["preconditioner_scale_cost"] = s.preconditioner_scale_cost;
  result["preconditioner_reuse_on_update"] = s.preconditioner_reuse_on_update;
  result["preconditioner_iter"] = s.preconditioner_iter;
  result["tau"] = s.tau;
  result["iterative_refinement_always_enabled"] = s.iterative_refinement_always_enabled;
  result["iterative_refinement_eps_abs"] = s.iterative_refinement_eps_abs;
  result["iterative_refinement_eps_rel"] = s.iterative_refinement_eps_rel;
  result["iterative_refinement_max_iter"] = s.iterative_refinement_max_iter;
  result["iterative_refinement_min_improvement_rate"] = s.iterative_refinement_min_improvement_rate;
  result["iterative_refinement_static_regularization_eps"] = s.iterative_refinement_static_regularization_eps;
  result["iterative_refinement_static_regularization_rel"] = s.iterative_refinement_static_regularization_rel;
  result["verbose"] = s.verbose;
  result["compute_timings"] = s.compute_timings;
  return result;
}

// [[Rcpp::export]]
void update_settings(SEXP solver_p, bool dense_backend, const Rcpp::List& settings) {
    if (settings.size() > 0) {
      if (dense_backend) {
        piqp_update_settings((Rcpp::as<Rcpp::XPtr<piqp::DenseSolver<double>>>(solver_p))->settings(), settings);
      } else {
        piqp_update_settings((Rcpp::as<Rcpp::XPtr<piqp::SparseSolver<double>>>(solver_p))->settings(), settings);
      }
    }
}

// [[Rcpp::export]]
SEXP piqp_dense_setup (Eigen::Map<Eigen::MatrixXd> P,
		       Eigen::Map<Eigen::VectorXd> c,
		       Eigen::Map<Eigen::MatrixXd> A,
		       Eigen::Map<Eigen::VectorXd> b,
		       Eigen::Map<Eigen::MatrixXd> G,
		       Eigen::Map<Eigen::VectorXd> h_l,
		       Eigen::Map<Eigen::VectorXd> h_u,
		       Eigen::Map<Eigen::VectorXd> x_l,
		       Eigen::Map<Eigen::VectorXd> x_u,
		       Rcpp::List settings) {
  piqp::DenseSolver<double>* solver = new piqp::DenseSolver<double>();

  if (settings.size()) piqp_update_settings(solver->settings(), settings);
  solver->setup(P, c, A, b, G, h_l, h_u, x_l, x_u);

  Rcpp::XPtr<piqp::DenseSolver<double>> ptr(solver);
  return ptr;
}

// [[Rcpp::export]]
SEXP piqp_sparse_setup (Eigen::Map<Eigen::SparseMatrix<double>> P,
			Eigen::Map<Eigen::VectorXd> c,
			Eigen::Map<Eigen::SparseMatrix<double>> A,
			Eigen::Map<Eigen::VectorXd> b,
			Eigen::Map<Eigen::SparseMatrix<double>> G,
			Eigen::Map<Eigen::VectorXd> h_l,
			Eigen::Map<Eigen::VectorXd> h_u,
			Eigen::Map<Eigen::VectorXd> x_l,
			Eigen::Map<Eigen::VectorXd> x_u,
			Rcpp::List settings) {
  piqp::SparseSolver<double>* solver = new piqp::SparseSolver<double>();

  if (settings.size() > 0) piqp_update_settings(solver->settings(), settings);
  solver->setup(P, c, A, b, G, h_l, h_u, x_l, x_u);

  Rcpp::XPtr<piqp::SparseSolver<double>> ptr(solver);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List solve_model(SEXP solver_p, bool dense_backend) {
  piqp::Result<double> result;
  if (dense_backend) {
    (Rcpp::as<Rcpp::XPtr<piqp::DenseSolver<double>>>(solver_p))->solve();
    result = (Rcpp::as<Rcpp::XPtr<piqp::DenseSolver<double>>>(solver_p))->result();
  } else {
    (Rcpp::as<Rcpp::XPtr<piqp::SparseSolver<double>>>(solver_p))->solve();
    result = (Rcpp::as<Rcpp::XPtr<piqp::SparseSolver<double>>>(solver_p))->result();
  }
  // Return solver result as R list
  Rcpp::List info = Rcpp::List::create();
  info["status_desc"] = piqp::status_to_string(result.info.status);
  info["iter"] = result.info.iter;
  info["rho"] = result.info.rho;
  info["delta"] = result.info.delta;
  info["mu"] = result.info.mu;
  info["sigma"] = result.info.sigma;
  info["primal_step"] = result.info.primal_step;
  info["dual_step"] = result.info.dual_step;
  info["primal_res"] = result.info.primal_res;
  info["primal_res_rel"] = result.info.primal_res_rel;
  info["dual_res"] = result.info.dual_res;
  info["dual_res_rel"] = result.info.dual_res_rel;
  info["primal_res_reg"] = result.info.primal_res_reg;
  info["primal_res_reg_rel"] = result.info.primal_res_reg_rel;
  info["dual_res_reg"] = result.info.dual_res_reg;
  info["dual_res_reg_rel"] = result.info.dual_res_reg_rel;
  info["primal_prox_inf"] = result.info.primal_prox_inf;
  info["dual_prox_inf"] = result.info.dual_prox_inf;
  info["primal_obj"] = result.info.primal_obj;
  info["dual_obj"] = result.info.dual_obj;
  info["duality_gap"] = result.info.duality_gap;
  info["duality_gap_rel"] = result.info.duality_gap_rel;
  info["factor_retires"] = result.info.factor_retires;
  info["reg_limit"] = result.info.reg_limit;
  info["no_primal_update"] = result.info.no_primal_update;
  info["no_dual_update"] = result.info.no_dual_update;
  info["setup_time"] = result.info.setup_time;
  info["update_time"] = result.info.update_time;
  info["solve_time"] = result.info.solve_time;
  info["kkt_factor_time"] = result.info.kkt_factor_time;
  info["kkt_solve_time"] = result.info.kkt_solve_time;
  info["run_time"] = result.info.run_time;

  return Rcpp::List::create(
			    Rcpp::_["status"] = (int) result.info.status,
			    Rcpp::_["x"] = result.x,
			    Rcpp::_["y"] = result.y,
			    Rcpp::_["z_l"] = result.z_l,
			    Rcpp::_["z_u"] = result.z_u,
			    Rcpp::_["z_bl"] = result.z_bl,
			    Rcpp::_["z_bu"] = result.z_bu,
			    Rcpp::_["s_l"] = result.s_l,
			    Rcpp::_["s_u"] = result.s_u,
			    Rcpp::_["s_bl"] = result.s_bl,
			    Rcpp::_["s_bu"] = result.s_bu,
			    Rcpp::_["info"] = info
			    );
}

template<typename T>
piqp::optional<T> nullable2optional(Rcpp::Nullable<T> data) {
    if (data.isNotNull()) {
        return piqp::optional<T>(Rcpp::as<T>(data.get()));
    }
    return piqp::optional<T>();
}

// [[Rcpp::export]]
void piqp_update_dense(SEXP solver_p,
                       Rcpp::Nullable<Eigen::Map<Mat>> P,
                       Rcpp::Nullable<Eigen::Map<Vec>> c,
                       Rcpp::Nullable<Eigen::Map<Mat>> A,
                       Rcpp::Nullable<Eigen::Map<Vec>> b,
                       Rcpp::Nullable<Eigen::Map<Mat>> G,
                       Rcpp::Nullable<Eigen::Map<Vec>> h_l,
                       Rcpp::Nullable<Eigen::Map<Vec>> h_u,
                       Rcpp::Nullable<Eigen::Map<Vec>> x_l,
                       Rcpp::Nullable<Eigen::Map<Vec>> x_u) {
    auto solver = Rcpp::as<Rcpp::XPtr<piqp::DenseSolver<double>>>(solver_p);
    solver->update(nullable2optional(P),
                   nullable2optional(c),
                   nullable2optional(A),
                   nullable2optional(b),
                   nullable2optional(G),
                   nullable2optional(h_l),
                   nullable2optional(h_u),
                   nullable2optional(x_l),
                   nullable2optional(x_u));
}

// [[Rcpp::export]]
void piqp_update_sparse(SEXP solver_p,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> P,
                        Rcpp::Nullable<Eigen::Map<Vec>> c,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> A,
                        Rcpp::Nullable<Eigen::Map<Vec>> b,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> G,
                        Rcpp::Nullable<Eigen::Map<Vec>> h_l,
                        Rcpp::Nullable<Eigen::Map<Vec>> h_u,
                        Rcpp::Nullable<Eigen::Map<Vec>> x_l,
                        Rcpp::Nullable<Eigen::Map<Vec>> x_u) {
    auto solver = Rcpp::as<Rcpp::XPtr<piqp::SparseSolver<double>>>(solver_p);
    solver->update(nullable2optional(P),
                   nullable2optional(c),
                   nullable2optional(A),
                   nullable2optional(b),
                   nullable2optional(G),
                   nullable2optional(h_l),
                   nullable2optional(h_u),
                   nullable2optional(x_l),
                   nullable2optional(x_u));
}
