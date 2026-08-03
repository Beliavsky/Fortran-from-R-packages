// SPDX-License-Identifier: Apache-2.0
#![allow(non_snake_case)]

use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettings, DefaultSolver, IPSolver, SupportedConeT,
    ExponentialConeT, GenPowerConeT, NonnegativeConeT, PowerConeT,
    PSDTriangleConeT, SecondOrderConeT, ZeroConeT,
};
use std::ffi::{c_char, c_void};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::ptr;
use std::slice;

#[repr(C)]
pub struct ClarabelCscC {
    nrows: usize,
    ncols: usize,
    nnz: usize,
    colptr: *const usize,
    rowind: *const usize,
    values: *const f64,
}

#[repr(C)]
pub struct ClarabelConeC {
    tag: u8,
    dim: usize,
    parameter: f64,
    alpha: *const f64,
    alpha_len: usize,
}

#[repr(C)]
pub struct ClarabelSettingsC {
    max_iter: i32,
    time_limit: f64,
    verbose: i32,
    max_step_fraction: f64,
    tol_gap_abs: f64,
    tol_gap_rel: f64,
    tol_feas: f64,
    tol_infeas_abs: f64,
    tol_infeas_rel: f64,
    tol_ktratio: f64,
    reduced_tol_gap_abs: f64,
    reduced_tol_gap_rel: f64,
    reduced_tol_feas: f64,
    reduced_tol_infeas_abs: f64,
    reduced_tol_infeas_rel: f64,
    reduced_tol_ktratio: f64,
    equilibrate_enable: i32,
    equilibrate_max_iter: i32,
    equilibrate_min_scaling: f64,
    equilibrate_max_scaling: f64,
    linesearch_backtrack_step: f64,
    min_switch_step_length: f64,
    min_terminate_step_length: f64,
    max_threads: i32,
    direct_kkt_solver: i32,
    direct_solve_method: i32,
    static_regularization_enable: i32,
    static_regularization_constant: f64,
    static_regularization_proportional: f64,
    dynamic_regularization_enable: i32,
    dynamic_regularization_eps: f64,
    dynamic_regularization_delta: f64,
    iterative_refinement_enable: i32,
    iterative_refinement_reltol: f64,
    iterative_refinement_abstol: f64,
    iterative_refinement_max_iter: i32,
    iterative_refinement_stop_ratio: f64,
    presolve_enable: i32,
    input_sparse_dropzeros: i32,
    chordal_decomposition_enable: i32,
    chordal_decomposition_merge_method: i32,
    chordal_decomposition_compact: i32,
    chordal_decomposition_complete_dual: i32,
}

#[repr(C)]
pub struct ClarabelResultC {
    status: i32,
    iterations: i32,
    obj_val: f64,
    obj_val_dual: f64,
    solve_time: f64,
    r_prim: f64,
    r_dual: f64,
    mu: f64,
    sigma: f64,
    step_length: f64,
    cost_primal: f64,
    cost_dual: f64,
    res_primal: f64,
    res_dual: f64,
    res_primal_inf: f64,
    res_dual_inf: f64,
    gap_abs: f64,
    gap_rel: f64,
    ktratio: f64,
    linear_solver_threads: usize,
    linear_solver_nnz_a: usize,
    linear_solver_nnz_l: usize,
}

fn bool_i32(x: bool) -> i32 { if x { 1 } else { 0 } }
fn bool_from_i32(x: i32) -> bool { x != 0 }

fn settings_to_c(s: &DefaultSettings<f64>) -> ClarabelSettingsC {
    ClarabelSettingsC {
        max_iter: s.max_iter as i32,
        time_limit: s.time_limit,
        verbose: bool_i32(s.verbose),
        max_step_fraction: s.max_step_fraction,
        tol_gap_abs: s.tol_gap_abs,
        tol_gap_rel: s.tol_gap_rel,
        tol_feas: s.tol_feas,
        tol_infeas_abs: s.tol_infeas_abs,
        tol_infeas_rel: s.tol_infeas_rel,
        tol_ktratio: s.tol_ktratio,
        reduced_tol_gap_abs: s.reduced_tol_gap_abs,
        reduced_tol_gap_rel: s.reduced_tol_gap_rel,
        reduced_tol_feas: s.reduced_tol_feas,
        reduced_tol_infeas_abs: s.reduced_tol_infeas_abs,
        reduced_tol_infeas_rel: s.reduced_tol_infeas_rel,
        reduced_tol_ktratio: s.reduced_tol_ktratio,
        equilibrate_enable: bool_i32(s.equilibrate_enable),
        equilibrate_max_iter: s.equilibrate_max_iter as i32,
        equilibrate_min_scaling: s.equilibrate_min_scaling,
        equilibrate_max_scaling: s.equilibrate_max_scaling,
        linesearch_backtrack_step: s.linesearch_backtrack_step,
        min_switch_step_length: s.min_switch_step_length,
        min_terminate_step_length: s.min_terminate_step_length,
        max_threads: s.max_threads as i32,
        direct_kkt_solver: bool_i32(s.direct_kkt_solver),
        direct_solve_method: match s.direct_solve_method.as_str() {
            "qdldl" => 1, "faer" => 2, "mkl" => 3, "panua" => 4, _ => 0
        },
        static_regularization_enable: bool_i32(s.static_regularization_enable),
        static_regularization_constant: s.static_regularization_constant,
        static_regularization_proportional: s.static_regularization_proportional,
        dynamic_regularization_enable: bool_i32(s.dynamic_regularization_enable),
        dynamic_regularization_eps: s.dynamic_regularization_eps,
        dynamic_regularization_delta: s.dynamic_regularization_delta,
        iterative_refinement_enable: bool_i32(s.iterative_refinement_enable),
        iterative_refinement_reltol: s.iterative_refinement_reltol,
        iterative_refinement_abstol: s.iterative_refinement_abstol,
        iterative_refinement_max_iter: s.iterative_refinement_max_iter as i32,
        iterative_refinement_stop_ratio: s.iterative_refinement_stop_ratio,
        presolve_enable: bool_i32(s.presolve_enable),
        input_sparse_dropzeros: bool_i32(s.input_sparse_dropzeros),
        chordal_decomposition_enable: bool_i32(s.chordal_decomposition_enable),
        chordal_decomposition_merge_method: match s.chordal_decomposition_merge_method.as_str() {
            "parent_child" => 1, "clique_graph" => 2, _ => 0
        },
        chordal_decomposition_compact: bool_i32(s.chordal_decomposition_compact),
        chordal_decomposition_complete_dual: bool_i32(s.chordal_decomposition_complete_dual),
    }
}

fn settings_from_c(s: &ClarabelSettingsC) -> DefaultSettings<f64> {
    let mut out = DefaultSettings::<f64>::default();
    out.max_iter = s.max_iter.max(0) as u32;
    out.time_limit = s.time_limit;
    out.verbose = bool_from_i32(s.verbose);
    out.max_step_fraction = s.max_step_fraction;
    out.tol_gap_abs = s.tol_gap_abs;
    out.tol_gap_rel = s.tol_gap_rel;
    out.tol_feas = s.tol_feas;
    out.tol_infeas_abs = s.tol_infeas_abs;
    out.tol_infeas_rel = s.tol_infeas_rel;
    out.tol_ktratio = s.tol_ktratio;
    out.reduced_tol_gap_abs = s.reduced_tol_gap_abs;
    out.reduced_tol_gap_rel = s.reduced_tol_gap_rel;
    out.reduced_tol_feas = s.reduced_tol_feas;
    out.reduced_tol_infeas_abs = s.reduced_tol_infeas_abs;
    out.reduced_tol_infeas_rel = s.reduced_tol_infeas_rel;
    out.reduced_tol_ktratio = s.reduced_tol_ktratio;
    out.equilibrate_enable = bool_from_i32(s.equilibrate_enable);
    out.equilibrate_max_iter = s.equilibrate_max_iter.max(0) as u32;
    out.equilibrate_min_scaling = s.equilibrate_min_scaling;
    out.equilibrate_max_scaling = s.equilibrate_max_scaling;
    out.linesearch_backtrack_step = s.linesearch_backtrack_step;
    out.min_switch_step_length = s.min_switch_step_length;
    out.min_terminate_step_length = s.min_terminate_step_length;
    out.max_threads = s.max_threads.max(0) as u32;
    out.direct_kkt_solver = bool_from_i32(s.direct_kkt_solver);
    out.direct_solve_method = match s.direct_solve_method {
        1 => "qdldl", 2 => "faer", 3 => "mkl", 4 => "panua", _ => "auto"
    }.to_string();
    out.static_regularization_enable = bool_from_i32(s.static_regularization_enable);
    out.static_regularization_constant = s.static_regularization_constant;
    out.static_regularization_proportional = s.static_regularization_proportional;
    out.dynamic_regularization_enable = bool_from_i32(s.dynamic_regularization_enable);
    out.dynamic_regularization_eps = s.dynamic_regularization_eps;
    out.dynamic_regularization_delta = s.dynamic_regularization_delta;
    out.iterative_refinement_enable = bool_from_i32(s.iterative_refinement_enable);
    out.iterative_refinement_reltol = s.iterative_refinement_reltol;
    out.iterative_refinement_abstol = s.iterative_refinement_abstol;
    out.iterative_refinement_max_iter = s.iterative_refinement_max_iter.max(0) as u32;
    out.iterative_refinement_stop_ratio = s.iterative_refinement_stop_ratio;
    out.presolve_enable = bool_from_i32(s.presolve_enable);
    out.input_sparse_dropzeros = bool_from_i32(s.input_sparse_dropzeros);
    out.chordal_decomposition_enable = bool_from_i32(s.chordal_decomposition_enable);
    out.chordal_decomposition_merge_method = match s.chordal_decomposition_merge_method {
        1 => "parent_child", 2 => "clique_graph", _ => "none"
    }.to_string();
    out.chordal_decomposition_compact = bool_from_i32(s.chordal_decomposition_compact);
    out.chordal_decomposition_complete_dual = bool_from_i32(s.chordal_decomposition_complete_dual);
    out
}

unsafe fn slice_or_empty<'a, T>(p: *const T, len: usize) -> Result<&'a [T], String> {
    if len == 0 { Ok(&[]) }
    else if p.is_null() { Err("null pointer with nonzero length".to_string()) }
    else { Ok(slice::from_raw_parts(p, len)) }
}

unsafe fn csc_from_c(c: &ClarabelCscC) -> Result<CscMatrix<f64>, String> {
    let colptr = slice_or_empty(c.colptr, c.ncols + 1)?.to_vec();
    let rowind = slice_or_empty(c.rowind, c.nnz)?.to_vec();
    let values = slice_or_empty(c.values, c.nnz)?.to_vec();
    let out = CscMatrix::new(c.nrows, c.ncols, colptr, rowind, values);
    out.check_format().map_err(|e| format!("CSC format error: {e:?}"))?;
    Ok(out)
}

unsafe fn cones_from_c(p: *const ClarabelConeC, n: usize) -> Result<Vec<SupportedConeT<f64>>, String> {
    let input = slice_or_empty(p, n)?;
    let mut out = Vec::with_capacity(n);
    for c in input {
        let cone = match c.tag {
            0 => ZeroConeT(c.dim),
            1 => NonnegativeConeT(c.dim),
            2 => SecondOrderConeT(c.dim),
            3 => ExponentialConeT(),
            4 => PowerConeT(c.parameter),
            5 => GenPowerConeT(slice_or_empty(c.alpha, c.alpha_len)?.to_vec(), c.dim),
            6 => PSDTriangleConeT(c.dim),
            _ => return Err(format!("unknown cone tag {}", c.tag)),
        };
        out.push(cone);
    }
    Ok(out)
}

unsafe fn write_error(buffer: *mut c_char, capacity: usize, message: &str) {
    if buffer.is_null() || capacity == 0 { return; }
    let bytes = message.as_bytes();
    let n = bytes.len().min(capacity - 1);
    ptr::copy_nonoverlapping(bytes.as_ptr(), buffer as *mut u8, n);
    *buffer.add(n) = 0;
}

fn panic_message(p: Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = p.downcast_ref::<&str>() { (*s).to_string() }
    else if let Some(s) = p.downcast_ref::<String>() { s.clone() }
    else { "Rust panic in Clarabel bridge".to_string() }
}

#[no_mangle]
pub extern "C" fn clarabel_settings_default(settings: *mut ClarabelSettingsC) {
    if settings.is_null() { return; }
    unsafe { *settings = settings_to_c(&DefaultSettings::<f64>::default()); }
}

#[no_mangle]
pub extern "C" fn clarabel_solver_create(
    P: *const ClarabelCscC, q: *const f64, q_len: usize,
    A: *const ClarabelCscC, b: *const f64, b_len: usize,
    cones: *const ClarabelConeC, ncones: usize,
    settings: *const ClarabelSettingsC, out: *mut *mut c_void,
    error: *mut c_char, error_capacity: usize,
) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<(), String> {
        if P.is_null() || A.is_null() || settings.is_null() || out.is_null() {
            return Err("null required argument".to_string());
        }
        let pm = unsafe { csc_from_c(&*P)? };
        let am = unsafe { csc_from_c(&*A)? };
        let qv = unsafe { slice_or_empty(q, q_len)?.to_vec() };
        let bv = unsafe { slice_or_empty(b, b_len)?.to_vec() };
        let cv = unsafe { cones_from_c(cones, ncones)? };
        let st = unsafe { settings_from_c(&*settings) };
        let solver = DefaultSolver::new(&pm, &qv, &am, &bv, &cv, st).map_err(|e| e.to_string())?;
        unsafe { *out = Box::into_raw(Box::new(solver)) as *mut c_void; }
        Ok(())
    }));
    match result {
        Ok(Ok(())) => 0,
        Ok(Err(msg)) => { unsafe { write_error(error, error_capacity, &msg) }; -1 }
        Err(p) => { unsafe { write_error(error, error_capacity, &panic_message(p)) }; -2 }
    }
}

#[no_mangle]
pub extern "C" fn clarabel_solver_solve(
    solver: *mut c_void, x: *mut f64, x_len: usize,
    z: *mut f64, z_len: usize, s: *mut f64, s_len: usize,
    result: *mut ClarabelResultC, error: *mut c_char, error_capacity: usize,
) -> i32 {
    let call = catch_unwind(AssertUnwindSafe(|| -> Result<(), String> {
        if solver.is_null() || result.is_null() { return Err("null solver or result".to_string()); }
        let solver = unsafe { &mut *(solver as *mut DefaultSolver<f64>) };
        solver.solve();
        if x_len != solver.solution.x.len() || z_len != solver.solution.z.len() || s_len != solver.solution.s.len() {
            return Err("output array length mismatch".to_string());
        }
        if x_len > 0 { unsafe { ptr::copy_nonoverlapping(solver.solution.x.as_ptr(), x, x_len) } }
        if z_len > 0 { unsafe { ptr::copy_nonoverlapping(solver.solution.z.as_ptr(), z, z_len) } }
        if s_len > 0 { unsafe { ptr::copy_nonoverlapping(solver.solution.s.as_ptr(), s, s_len) } }
        let sol = &solver.solution;
        let info = &solver.info;
        unsafe {
            *result = ClarabelResultC {
                status: sol.status as i32,
                iterations: sol.iterations as i32,
                obj_val: sol.obj_val,
                obj_val_dual: sol.obj_val_dual,
                solve_time: sol.solve_time,
                r_prim: sol.r_prim,
                r_dual: sol.r_dual,
                mu: info.mu,
                sigma: info.sigma,
                step_length: info.step_length,
                cost_primal: info.cost_primal,
                cost_dual: info.cost_dual,
                res_primal: info.res_primal,
                res_dual: info.res_dual,
                res_primal_inf: info.res_primal_inf,
                res_dual_inf: info.res_dual_inf,
                gap_abs: info.gap_abs,
                gap_rel: info.gap_rel,
                ktratio: info.ktratio,
                linear_solver_threads: info.linsolver.threads,
                linear_solver_nnz_a: info.linsolver.nnzA,
                linear_solver_nnz_l: info.linsolver.nnzL,
            };
        }
        Ok(())
    }));
    match call {
        Ok(Ok(())) => 0,
        Ok(Err(msg)) => { unsafe { write_error(error, error_capacity, &msg) }; -1 }
        Err(p) => { unsafe { write_error(error, error_capacity, &panic_message(p)) }; -2 }
    }
}

#[no_mangle]
pub extern "C" fn clarabel_solver_update(
    solver: *mut c_void,
    p_values: *const f64, p_len: usize,
    a_values: *const f64, a_len: usize,
    q: *const f64, q_len: usize,
    b: *const f64, b_len: usize,
    error: *mut c_char, error_capacity: usize,
) -> i32 {
    let call = catch_unwind(AssertUnwindSafe(|| -> Result<(), String> {
        if solver.is_null() { return Err("null solver".to_string()); }
        let solver = unsafe { &mut *(solver as *mut DefaultSolver<f64>) };
        if p_len > 0 { solver.update_P(unsafe { slice_or_empty(p_values, p_len)? }).map_err(|e| e.to_string())?; }
        if a_len > 0 { solver.update_A(unsafe { slice_or_empty(a_values, a_len)? }).map_err(|e| e.to_string())?; }
        if q_len > 0 { solver.update_q(unsafe { slice_or_empty(q, q_len)? }).map_err(|e| e.to_string())?; }
        if b_len > 0 { solver.update_b(unsafe { slice_or_empty(b, b_len)? }).map_err(|e| e.to_string())?; }
        Ok(())
    }));
    match call {
        Ok(Ok(())) => 0,
        Ok(Err(msg)) => { unsafe { write_error(error, error_capacity, &msg) }; -1 }
        Err(p) => { unsafe { write_error(error, error_capacity, &panic_message(p)) }; -2 }
    }
}

#[no_mangle]
pub extern "C" fn clarabel_solver_is_update_allowed(solver: *const c_void) -> i32 {
    if solver.is_null() { return 0; }
    let solver = unsafe { &*(solver as *const DefaultSolver<f64>) };
    bool_i32(solver.is_data_update_allowed())
}

#[no_mangle]
pub extern "C" fn clarabel_solver_free(solver: *mut c_void) {
    if solver.is_null() { return; }
    unsafe { drop(Box::from_raw(solver as *mut DefaultSolver<f64>)); }
}
