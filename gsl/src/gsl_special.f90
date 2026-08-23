module gsl_special
  use iso_c_binding, only: c_ptr, c_loc, c_double, c_int, c_size_t
  implicit none
  private

  public :: airy_ai
  public :: airy_bi
  public :: airy_ai_scaled
  public :: airy_bi_scaled
  public :: airy_ai_deriv
  public :: airy_bi_deriv
  public :: airy_ai_deriv_scaled
  public :: airy_bi_deriv_scaled
  public :: airy_zero_ai
  public :: airy_zero_bi
  public :: airy_zero_ai_deriv
  public :: airy_zero_bi_deriv
  public :: bessel_cyl_j0
  public :: bessel_cyl_j1
  public :: bessel_cyl_jn
  public :: bessel_cyl_jn_array
  public :: bessel_cyl_y0
  public :: bessel_cyl_y1
  public :: bessel_cyl_yn
  public :: bessel_cyl_yn_array
  public :: bessel_mod_i0
  public :: bessel_mod_i1
  public :: bessel_mod_in
  public :: bessel_mod_in_array
  public :: bessel_mod_i0_scaled
  public :: bessel_mod_i1_scaled
  public :: bessel_mod_in_scaled
  public :: bessel_mod_in_scaled_array
  public :: bessel_mod_k0
  public :: bessel_mod_k1
  public :: bessel_mod_kn
  public :: bessel_mod_kn_array
  public :: bessel_mod_k0_scaled
  public :: bessel_mod_k1_scaled
  public :: bessel_mod_kn_scaled
  public :: bessel_mod_kn_scaled_array
  public :: bessel_sph_j0
  public :: bessel_sph_j1
  public :: bessel_sph_j2
  public :: bessel_sph_jl
  public :: bessel_sph_jl_array
  public :: bessel_sph_jl_steed_array
  public :: bessel_sph_y0
  public :: bessel_sph_y1
  public :: bessel_sph_y2
  public :: bessel_sph_yl
  public :: bessel_sph_yl_array
  public :: bessel_sph_i0_scaled
  public :: bessel_sph_i1_scaled
  public :: bessel_sph_i2_scaled
  public :: bessel_sph_il_scaled
  public :: bessel_sph_il_scaled_array
  public :: bessel_sph_k0_scaled
  public :: bessel_sph_k1_scaled
  public :: bessel_sph_k2_scaled
  public :: bessel_sph_kl_scaled
  public :: bessel_sph_kl_scaled_array
  public :: bessel_cyl_jnu
  public :: bessel_sequence_jnu
  public :: bessel_cyl_ynu
  public :: bessel_mod_inu
  public :: bessel_mod_inu_scaled
  public :: bessel_mod_knu
  public :: bessel_lnknu
  public :: bessel_mod_knu_scaled
  public :: bessel_zero_j0
  public :: bessel_zero_j1
  public :: bessel_zero_jnu
  public :: clausen
  public :: hydrogenicr_1
  public :: hydrogenicr
  public :: coulomb_wave_fg
  public :: coulomb_wave_f_array
  public :: coulomb_wave_fg_array
  public :: coulomb_wave_fgp_array
  public :: coulomb_wave_sphf_array
  public :: coulomb_cl
  public :: coulomb_cl_array
  public :: coupling_3j
  public :: coupling_6j
  public :: coupling_9j
  public :: dawson
  public :: debye_1
  public :: debye_2
  public :: debye_3
  public :: debye_4
  public :: dilog
  public :: complex_dilog
  public :: ellint_kcomp
  public :: ellint_ecomp
  public :: ellint_f
  public :: ellint_e
  public :: ellint_p
  public :: ellint_d
  public :: ellint_rc
  public :: ellint_rd
  public :: ellint_rf
  public :: ellint_rj
  public :: elljac
  public :: erf
  public :: erfc
  public :: log_erfc
  public :: erf_z
  public :: erf_q
  public :: hazard
  public :: expint_e1
  public :: expint_e2
  public :: expint_en
  public :: expint_ei
  public :: shi
  public :: chi
  public :: expint_3
  public :: si
  public :: ci
  public :: atanint
  public :: fermi_dirac_m1
  public :: fermi_dirac_0
  public :: fermi_dirac_1
  public :: fermi_dirac_2
  public :: fermi_dirac_int
  public :: fermi_dirac_mhalf
  public :: fermi_dirac_half
  public :: fermi_dirac_3half
  public :: fermi_dirac_inc_0
  public :: gsl_sf_gamma
  public :: lngamma
  public :: lngamma_sgn
  public :: gammastar
  public :: gammainv
  public :: lngamma_complex
  public :: taylorcoeff
  public :: fact
  public :: doublefact
  public :: lnfact
  public :: lndoublefact
  public :: gsl_sf_choose
  public :: lnchoose
  public :: poch
  public :: lnpoch
  public :: lnpoch_sgn
  public :: pochrel
  public :: gamma_inc_p
  public :: gamma_inc_q
  public :: gamma_inc
  public :: gsl_sf_beta
  public :: lnbeta
  public :: beta_inc
  public :: gegenpoly_1
  public :: gegenpoly_2
  public :: gegenpoly_3
  public :: gegenpoly_n
  public :: gegenpoly_array
  public :: hyperg_0f1
  public :: hyperg_1f1_int
  public :: hyperg_1f1
  public :: hyperg_u_int
  public :: hyperg_u
  public :: hyperg_2f1
  public :: hyperg_2f1_conj
  public :: hyperg_2f1_renorm
  public :: hyperg_2f1_conj_renorm
  public :: hyperg_2f0
  public :: laguerre_1
  public :: laguerre_2
  public :: laguerre_3
  public :: laguerre_n
  public :: lambert_w0
  public :: lambert_wm1
  public :: legendre_p1_raw
  public :: legendre_p1
  public :: legendre_p2_raw
  public :: legendre_p2
  public :: legendre_p3_raw
  public :: legendre_p3
  public :: legendre_pl_raw
  public :: legendre_pl
  public :: legendre_pl_array
  public :: legendre_q0_raw
  public :: legendre_q0
  public :: legendre_q1_raw
  public :: legendre_q1
  public :: legendre_ql_raw
  public :: legendre_ql
  public :: legendre_array_n
  public :: legendre_array_index
  public :: legendre_array
  public :: legendre_deriv_array
  public :: legendre_deriv_alt_array
  public :: legendre_deriv2_array
  public :: legendre_deriv2_alt_array
  public :: legendre_plm_raw
  public :: legendre_plm
  public :: legendre_sphplm_raw
  public :: legendre_sphplm
  public :: conicalp_half
  public :: conicalp_mhalf
  public :: conicalp_0
  public :: conicalp_1
  public :: conicalp_sph_reg
  public :: conicalp_cyl_reg
  public :: legendre_h3d_0
  public :: legendre_h3d_1
  public :: legendre_h3d
  public :: legendre_h3d_array
  public :: gsl_sf_log
  public :: log_abs
  public :: complex_log
  public :: log_1plusx
  public :: log_1plusx_mx
  public :: gsl_poly_c
  public :: pow_int
  public :: psi_int
  public :: psi
  public :: psi_1piy
  public :: psi_1_int
  public :: psi_1
  public :: psi_n
  public :: synchrotron_1
  public :: synchrotron_2
  public :: transport_2
  public :: transport_3
  public :: transport_4
  public :: transport_5
  public :: gsl_sf_sin
  public :: gsl_sf_cos
  public :: hypot
  public :: sinc
  public :: complex_sin
  public :: complex_cos
  public :: complex_logsin
  public :: lnsinh
  public :: lncosh
  public :: zeta_int
  public :: zeta
  public :: zetam1_int
  public :: zetam1
  public :: hzeta
  public :: eta_int
  public :: eta

  interface
    subroutine cshim_001( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Ai_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_001
    subroutine cshim_002( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Bi_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_002
    subroutine cshim_003( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Ai_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_003
    subroutine cshim_004( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Bi_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_004
    subroutine cshim_005( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Ai_deriv_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_005
    subroutine cshim_006( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Bi_deriv_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_006
    subroutine cshim_007( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Ai_deriv_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_007
    subroutine cshim_008( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="airy_Bi_deriv_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_008
    subroutine cshim_009( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="airy_zero_Ai_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_009
    subroutine cshim_010( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="airy_zero_Bi_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_010
    subroutine cshim_011( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="airy_zero_Ai_deriv_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_011
    subroutine cshim_012( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="airy_zero_Bi_deriv_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_012
    subroutine cshim_013( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_J0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_013
    subroutine cshim_014( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_J1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_014
    subroutine cshim_015( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Jn_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_015
    subroutine cshim_016( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_Jn_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_016
    subroutine cshim_017( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Y0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_017
    subroutine cshim_018( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Y1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_018
    subroutine cshim_019( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Yn_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_019
    subroutine cshim_020( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_Yn_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_020
    subroutine cshim_021( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_I0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_021
    subroutine cshim_022( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_I1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_022
    subroutine cshim_023( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_In_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_023
    subroutine cshim_024( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_In_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_024
    subroutine cshim_025( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_I0_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_025
    subroutine cshim_026( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_I1_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_026
    subroutine cshim_027( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_In_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_027
    subroutine cshim_028( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_In_scaled_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_028
    subroutine cshim_029( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_K0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_029
    subroutine cshim_030( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_K1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_030
    subroutine cshim_031( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Kn_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_031
    subroutine cshim_032( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_Kn_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_032
    subroutine cshim_033( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_K0_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_033
    subroutine cshim_034( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_K1_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_034
    subroutine cshim_035( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Kn_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_035
    subroutine cshim_036( &
      & nmin, &
      & nmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_Kn_scaled_array_e")
      import :: c_ptr
      type(c_ptr), value :: nmin
      type(c_ptr), value :: nmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_036
    subroutine cshim_037( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_j0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_037
    subroutine cshim_038( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_j1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_038
    subroutine cshim_039( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_j2_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_039
    subroutine cshim_040( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_jl_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_040
    subroutine cshim_041( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_jl_array_e")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_041
    subroutine cshim_042( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_jl_steed_array_e")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_042
    subroutine cshim_043( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_y0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_043
    subroutine cshim_044( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_y1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_044
    subroutine cshim_045( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_y2_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_045
    subroutine cshim_046( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_yl_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_046
    subroutine cshim_047( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_yl_array_e")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_047
    subroutine cshim_048( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_i0_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_048
    subroutine cshim_049( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_i1_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_049
    subroutine cshim_050( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_i2_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_050
    subroutine cshim_051( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_il_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_051
    subroutine cshim_052( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_il_scaled_array_e")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_052
    subroutine cshim_053( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_k0_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_053
    subroutine cshim_054( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_k1_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_054
    subroutine cshim_055( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_k2_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_055
    subroutine cshim_056( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_kl_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_056
    subroutine cshim_057( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="bessel_kl_scaled_array_e")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_057
    subroutine cshim_058( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Jnu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_058
    subroutine cshim_059( &
      & nu, &
      & v, &
      & nv, &
      & mode, &
      & status) bind(C, name="bessel_sequence_Jnu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: v
      type(c_ptr), value :: nv
      type(c_ptr), value :: mode
      type(c_ptr), value :: status
    end subroutine cshim_059
    subroutine cshim_060( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Ynu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_060
    subroutine cshim_061( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Inu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_061
    subroutine cshim_062( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Inu_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_062
    subroutine cshim_063( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Knu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_063
    subroutine cshim_064( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_lnKnu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_064
    subroutine cshim_065( &
      & nu, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_Knu_scaled_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_065
    subroutine cshim_066( &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_zero_J0_e")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_066
    subroutine cshim_067( &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_zero_J1_e")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_067
    subroutine cshim_068( &
      & nu, &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="bessel_zero_Jnu_e")
      import :: c_ptr
      type(c_ptr), value :: nu
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_068
    subroutine cshim_069( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="clausen_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_069
    subroutine cshim_070( &
      & Z, &
      & r, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hydrogenicR_1")
      import :: c_ptr
      type(c_ptr), value :: Z
      type(c_ptr), value :: r
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_070
    subroutine cshim_071( &
      & n, &
      & l, &
      & Z, &
      & r, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hydrogenicR")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: l
      type(c_ptr), value :: Z
      type(c_ptr), value :: r
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_071
    subroutine cshim_072( &
      & eta, &
      & x, &
      & L_F, &
      & k, &
      & len, &
      & val_F, &
      & err_F, &
      & val_Fp, &
      & err_Fp, &
      & val_G, &
      & err_G, &
      & val_Gp, &
      & err_Gp, &
      & exp_F, &
      & exp_G, &
      & status) bind(C, name="coulomb_wave_FG")
      import :: c_ptr
      type(c_ptr), value :: eta
      type(c_ptr), value :: x
      type(c_ptr), value :: L_F
      type(c_ptr), value :: k
      type(c_ptr), value :: len
      type(c_ptr), value :: val_F
      type(c_ptr), value :: err_F
      type(c_ptr), value :: val_Fp
      type(c_ptr), value :: err_Fp
      type(c_ptr), value :: val_G
      type(c_ptr), value :: err_G
      type(c_ptr), value :: val_Gp
      type(c_ptr), value :: err_Gp
      type(c_ptr), value :: exp_F
      type(c_ptr), value :: exp_G
      type(c_ptr), value :: status
    end subroutine cshim_072
    subroutine cshim_073( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & len, &
      & fc_array, &
      & F_exponent, &
      & status) bind(C, name="coulomb_wave_F_array")
      import :: c_ptr
      type(c_ptr), value :: L_min
      type(c_ptr), value :: kmax
      type(c_ptr), value :: eta
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: fc_array
      type(c_ptr), value :: F_exponent
      type(c_ptr), value :: status
    end subroutine cshim_073
    subroutine cshim_074( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & len, &
      & fc_array, &
      & gc_array, &
      & F_exponent, &
      & G_exponent, &
      & status) bind(C, name="coulomb_wave_FG_array")
      import :: c_ptr
      type(c_ptr), value :: L_min
      type(c_ptr), value :: kmax
      type(c_ptr), value :: eta
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: fc_array
      type(c_ptr), value :: gc_array
      type(c_ptr), value :: F_exponent
      type(c_ptr), value :: G_exponent
      type(c_ptr), value :: status
    end subroutine cshim_074
    subroutine cshim_075( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & len, &
      & fc_array, &
      & fcp_array, &
      & gc_array, &
      & gcp_array, &
      & F_exponent, &
      & G_exponent, &
      & status) bind(C, name="coulomb_wave_FGp_array")
      import :: c_ptr
      type(c_ptr), value :: L_min
      type(c_ptr), value :: kmax
      type(c_ptr), value :: eta
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: fc_array
      type(c_ptr), value :: fcp_array
      type(c_ptr), value :: gc_array
      type(c_ptr), value :: gcp_array
      type(c_ptr), value :: F_exponent
      type(c_ptr), value :: G_exponent
      type(c_ptr), value :: status
    end subroutine cshim_075
    subroutine cshim_076( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & len, &
      & fc_array, &
      & F_exponent, &
      & status) bind(C, name="coulomb_wave_sphF_array")
      import :: c_ptr
      type(c_ptr), value :: L_min
      type(c_ptr), value :: kmax
      type(c_ptr), value :: eta
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: fc_array
      type(c_ptr), value :: F_exponent
      type(c_ptr), value :: status
    end subroutine cshim_076
    subroutine cshim_077( &
      & L, &
      & eta, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="coulomb_CL")
      import :: c_ptr
      type(c_ptr), value :: L
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_077
    subroutine cshim_078( &
      & L_min, &
      & kmax, &
      & eta, &
      & len, &
      & cl, &
      & status) bind(C, name="coulomb_CL_array")
      import :: c_ptr
      type(c_ptr), value :: L_min
      type(c_ptr), value :: kmax
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: cl
      type(c_ptr), value :: status
    end subroutine cshim_078
    subroutine cshim_079( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_ma, &
      & two_mb, &
      & two_mc, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="coupling_3j")
      import :: c_ptr
      type(c_ptr), value :: two_ja
      type(c_ptr), value :: two_jb
      type(c_ptr), value :: two_jc
      type(c_ptr), value :: two_ma
      type(c_ptr), value :: two_mb
      type(c_ptr), value :: two_mc
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_079
    subroutine cshim_080( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_jd, &
      & two_je, &
      & two_jf, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="coupling_6j")
      import :: c_ptr
      type(c_ptr), value :: two_ja
      type(c_ptr), value :: two_jb
      type(c_ptr), value :: two_jc
      type(c_ptr), value :: two_jd
      type(c_ptr), value :: two_je
      type(c_ptr), value :: two_jf
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_080
    subroutine cshim_081( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_jd, &
      & two_je, &
      & two_jf, &
      & two_jg, &
      & two_jh, &
      & two_ji, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="coupling_9j")
      import :: c_ptr
      type(c_ptr), value :: two_ja
      type(c_ptr), value :: two_jb
      type(c_ptr), value :: two_jc
      type(c_ptr), value :: two_jd
      type(c_ptr), value :: two_je
      type(c_ptr), value :: two_jf
      type(c_ptr), value :: two_jg
      type(c_ptr), value :: two_jh
      type(c_ptr), value :: two_ji
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_081
    subroutine cshim_082( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="dawson")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_082
    subroutine cshim_083( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="debye_1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_083
    subroutine cshim_084( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="debye_2")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_084
    subroutine cshim_085( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="debye_3")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_085
    subroutine cshim_086( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="debye_4")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_086
    subroutine cshim_087( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="dilog_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_087
    subroutine cshim_088( &
      & r, &
      & theta, &
      & len, &
      & val_re, &
      & val_im, &
      & err_re, &
      & err_im, &
      & status) bind(C, name="complex_dilog_e")
      import :: c_ptr
      type(c_ptr), value :: r
      type(c_ptr), value :: theta
      type(c_ptr), value :: len
      type(c_ptr), value :: val_re
      type(c_ptr), value :: val_im
      type(c_ptr), value :: err_re
      type(c_ptr), value :: err_im
      type(c_ptr), value :: status
    end subroutine cshim_088
    subroutine cshim_089( &
      & k, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_Kcomp_e")
      import :: c_ptr
      type(c_ptr), value :: k
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_089
    subroutine cshim_090( &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_Ecomp_e")
      import :: c_ptr
      type(c_ptr), value :: k
      type(c_ptr), value :: nk
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_090
    subroutine cshim_091( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_F_e")
      import :: c_ptr
      type(c_ptr), value :: phi
      type(c_ptr), value :: k
      type(c_ptr), value :: nk
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_091
    subroutine cshim_092( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_E_e")
      import :: c_ptr
      type(c_ptr), value :: phi
      type(c_ptr), value :: k
      type(c_ptr), value :: nk
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_092
    subroutine cshim_093( &
      & phi, &
      & k, &
      & n, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_P_e")
      import :: c_ptr
      type(c_ptr), value :: phi
      type(c_ptr), value :: k
      type(c_ptr), value :: n
      type(c_ptr), value :: nk
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_093
    subroutine cshim_094( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_D_e")
      import :: c_ptr
      type(c_ptr), value :: phi
      type(c_ptr), value :: k
      type(c_ptr), value :: nk
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_094
    subroutine cshim_095( &
      & x, &
      & y, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_RC_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: y
      type(c_ptr), value :: nx
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_095
    subroutine cshim_096( &
      & x, &
      & y, &
      & z, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_RD_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: y
      type(c_ptr), value :: z
      type(c_ptr), value :: nx
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_096
    subroutine cshim_097( &
      & x, &
      & y, &
      & z, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_RF_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: y
      type(c_ptr), value :: z
      type(c_ptr), value :: nx
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_097
    subroutine cshim_098( &
      & x, &
      & y, &
      & z, &
      & p, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="ellint_RJ_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: y
      type(c_ptr), value :: z
      type(c_ptr), value :: p
      type(c_ptr), value :: nx
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_098
    subroutine cshim_099( &
      & u, &
      & m, &
      & len, &
      & sn, &
      & cn, &
      & dn, &
      & status) bind(C, name="elljac_e")
      import :: c_ptr
      type(c_ptr), value :: u
      type(c_ptr), value :: m
      type(c_ptr), value :: len
      type(c_ptr), value :: sn
      type(c_ptr), value :: cn
      type(c_ptr), value :: dn
      type(c_ptr), value :: status
    end subroutine cshim_099
    subroutine cshim_100( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="erf_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_100
    subroutine cshim_101( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="erfc_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_101
    subroutine cshim_102( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="log_erfc_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_102
    subroutine cshim_103( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="erf_Z_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_103
    subroutine cshim_104( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="erf_Q_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_104
    subroutine cshim_105( &
      & x, &
      & len, &
      & mode, &
      & val, &
      & err, &
      & status) bind(C, name="hazard_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: mode
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_105
    subroutine cshim_106( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="expint_E1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_106
    subroutine cshim_107( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="expint_E2_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_107
    subroutine cshim_108( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="expint_En_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_108
    subroutine cshim_109( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="expint_Ei_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_109
    subroutine cshim_110( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="Shi_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_110
    subroutine cshim_111( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="Chi_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_111
    subroutine cshim_112( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="expint_3_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_112
    subroutine cshim_113( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="Si_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_113
    subroutine cshim_114( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="Ci_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_114
    subroutine cshim_115( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="atanint_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_115
    subroutine cshim_116( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_m1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_116
    subroutine cshim_117( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_0")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_117
    subroutine cshim_118( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_118
    subroutine cshim_119( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_2")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_119
    subroutine cshim_120( &
      & j, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_int")
      import :: c_ptr
      type(c_ptr), value :: j
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_120
    subroutine cshim_121( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_mhalf")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_121
    subroutine cshim_122( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_half")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_122
    subroutine cshim_123( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_3half")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_123
    subroutine cshim_124( &
      & x, &
      & b, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fermi_dirac_inc_0")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: b
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_124
    subroutine cshim_125( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gamma_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_125
    subroutine cshim_126( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lngamma_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_126
    subroutine cshim_127( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status, &
      & sgn) bind(C, name="lngamma_sgn_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
      type(c_ptr), value :: sgn
    end subroutine cshim_127
    subroutine cshim_128( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gammastar_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_128
    subroutine cshim_129( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gammainv_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_129
    subroutine cshim_130( &
      & zr, &
      & zi, &
      & len, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status) bind(C, name="lngamma_complex_e")
      import :: c_ptr
      type(c_ptr), value :: zr
      type(c_ptr), value :: zi
      type(c_ptr), value :: len
      type(c_ptr), value :: val_lnr
      type(c_ptr), value :: val_arg
      type(c_ptr), value :: err_lnr
      type(c_ptr), value :: err_arg
      type(c_ptr), value :: status
    end subroutine cshim_130
    subroutine cshim_131( &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="taylorcoeff_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_131
    subroutine cshim_132( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="fact_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_132
    subroutine cshim_133( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="doublefact_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_133
    subroutine cshim_134( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lnfact_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_134
    subroutine cshim_135( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lndoublefact_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_135
    subroutine cshim_136( &
      & n, &
      & m, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="choose_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: m
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_136
    subroutine cshim_137( &
      & n, &
      & m, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lnchoose_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: m
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_137
    subroutine cshim_138( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="poch_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_138
    subroutine cshim_139( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lnpoch_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_139
    subroutine cshim_140( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status, &
      & sgn) bind(C, name="lnpoch_sgn_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
      type(c_ptr), value :: sgn
    end subroutine cshim_140
    subroutine cshim_141( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="pochrel_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_141
    subroutine cshim_142( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gamma_inc_P_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_142
    subroutine cshim_143( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gamma_inc_Q_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_143
    subroutine cshim_144( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gamma_inc_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_144
    subroutine cshim_145( &
      & a, &
      & b, &
      & nb, &
      & val, &
      & err, &
      & status) bind(C, name="beta_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: nb
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_145
    subroutine cshim_146( &
      & a, &
      & b, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lnbeta_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_146
    subroutine cshim_147( &
      & a, &
      & b, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="beta_inc_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_147
    subroutine cshim_148( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gegenpoly_1_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_148
    subroutine cshim_149( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gegenpoly_2_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_149
    subroutine cshim_150( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gegenpoly_3_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_150
    subroutine cshim_151( &
      & n, &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="gegenpoly_n_e")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_151
    subroutine cshim_152( &
      & nmax, &
      & lambda, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="gegenpoly_array")
      import :: c_ptr
      type(c_ptr), value :: nmax
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_152
    subroutine cshim_153( &
      & c, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_0F1_e")
      import :: c_ptr
      type(c_ptr), value :: c
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_153
    subroutine cshim_154( &
      & m, &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_1F1_int_e")
      import :: c_ptr
      type(c_ptr), value :: m
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_154
    subroutine cshim_155( &
      & a, &
      & b, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_1F1_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_155
    subroutine cshim_156( &
      & m, &
      & n, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_U_int_e")
      import :: c_ptr
      type(c_ptr), value :: m
      type(c_ptr), value :: n
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_156
    subroutine cshim_157( &
      & a, &
      & b, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_U_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_157
    subroutine cshim_158( &
      & a, &
      & b, &
      & c, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_2F1_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: c
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_158
    subroutine cshim_159( &
      & aR, &
      & aI, &
      & c, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_2F1_conj_e")
      import :: c_ptr
      type(c_ptr), value :: aR
      type(c_ptr), value :: aI
      type(c_ptr), value :: c
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_159
    subroutine cshim_160( &
      & a, &
      & b, &
      & c, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_2F1_renorm_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: c
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_160
    subroutine cshim_161( &
      & aR, &
      & aI, &
      & c, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_2F1_conj_renorm_e")
      import :: c_ptr
      type(c_ptr), value :: aR
      type(c_ptr), value :: aI
      type(c_ptr), value :: c
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_161
    subroutine cshim_162( &
      & a, &
      & b, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hyperg_2F0_e")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: b
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_162
    subroutine cshim_163( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="laguerre_1")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_163
    subroutine cshim_164( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="laguerre_2")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_164
    subroutine cshim_165( &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="laguerre_3")
      import :: c_ptr
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_165
    subroutine cshim_166( &
      & n, &
      & a, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="laguerre_n")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: a
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_166
    subroutine cshim_167( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lambert_W0")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_167
    subroutine cshim_168( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lambert_Wm1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_168
    subroutine cshim_169( &
      & x, &
      & len, &
      & out) bind(C, name="legendre_P1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_169
    subroutine cshim_170( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_P1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_170
    subroutine cshim_171( &
      & x, &
      & len, &
      & out) bind(C, name="legendre_P2")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_171
    subroutine cshim_172( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_P2_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_172
    subroutine cshim_173( &
      & x, &
      & len, &
      & out) bind(C, name="legendre_P3")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_173
    subroutine cshim_174( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_P3_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_174
    subroutine cshim_175( &
      & l, &
      & x, &
      & len, &
      & out) bind(C, name="legendre_Pl")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_175
    subroutine cshim_176( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_Pl_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_176
    subroutine cshim_177( &
      & lmax, &
      & x, &
      & len, &
      & out, &
      & status) bind(C, name="legendre_Pl_array")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_177
    subroutine cshim_178( &
      & x, &
      & len, &
      & out) bind(C, name="legendre_Q0")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_178
    subroutine cshim_179( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_Q0_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_179
    subroutine cshim_180( &
      & x, &
      & len, &
      & out) bind(C, name="legendre_Q1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_180
    subroutine cshim_181( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_Q1_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_181
    subroutine cshim_182( &
      & l, &
      & x, &
      & len, &
      & out) bind(C, name="legendre_Ql")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_182
    subroutine cshim_183( &
      & l, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_Ql_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_183
    subroutine cshim_184( &
      & lmax, &
      & ans) bind(C, name="legendre_array_n")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: ans
    end subroutine cshim_184
    subroutine cshim_185( &
      & l, &
      & m, &
      & ans) bind(C, name="legendre_array_index")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: m
      type(c_ptr), value :: ans
    end subroutine cshim_185
    subroutine cshim_186( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & out) bind(C, name="legendre_array")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: lmax
      type(c_ptr), value :: norm
      type(c_ptr), value :: csphase
      type(c_ptr), value :: result_array
      type(c_ptr), value :: out
    end subroutine cshim_186
    subroutine cshim_187( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & out) bind(C, name="legendre_deriv_array")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: lmax
      type(c_ptr), value :: norm
      type(c_ptr), value :: csphase
      type(c_ptr), value :: result_array
      type(c_ptr), value :: result_deriv_array
      type(c_ptr), value :: out
    end subroutine cshim_187
    subroutine cshim_188( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & out) bind(C, name="legendre_deriv_alt_array")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: lmax
      type(c_ptr), value :: norm
      type(c_ptr), value :: csphase
      type(c_ptr), value :: result_array
      type(c_ptr), value :: result_deriv_array
      type(c_ptr), value :: out
    end subroutine cshim_188
    subroutine cshim_189( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & result_deriv2_array, &
      & out) bind(C, name="legendre_deriv2_array")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: lmax
      type(c_ptr), value :: norm
      type(c_ptr), value :: csphase
      type(c_ptr), value :: result_array
      type(c_ptr), value :: result_deriv_array
      type(c_ptr), value :: result_deriv2_array
      type(c_ptr), value :: out
    end subroutine cshim_189
    subroutine cshim_190( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & result_deriv2_array, &
      & out) bind(C, name="legendre_deriv2_alt_array")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: lmax
      type(c_ptr), value :: norm
      type(c_ptr), value :: csphase
      type(c_ptr), value :: result_array
      type(c_ptr), value :: result_deriv_array
      type(c_ptr), value :: result_deriv2_array
      type(c_ptr), value :: out
    end subroutine cshim_190
    subroutine cshim_191( &
      & l, &
      & m, &
      & x, &
      & len, &
      & out) bind(C, name="legendre_Plm")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: m
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_191
    subroutine cshim_192( &
      & l, &
      & m, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_Plm_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: m
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_192
    subroutine cshim_193( &
      & l, &
      & m, &
      & x, &
      & len, &
      & out) bind(C, name="legendre_sphPlm")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: m
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: out
    end subroutine cshim_193
    subroutine cshim_194( &
      & l, &
      & m, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_sphPlm_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: m
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_194
    subroutine cshim_195( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_half_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_195
    subroutine cshim_196( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_mhalf_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_196
    subroutine cshim_197( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_0_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_197
    subroutine cshim_198( &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_1_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_198
    subroutine cshim_199( &
      & l, &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_sph_reg_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_199
    subroutine cshim_200( &
      & m, &
      & lambda, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="conicalP_cyl_reg_e")
      import :: c_ptr
      type(c_ptr), value :: m
      type(c_ptr), value :: lambda
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_200
    subroutine cshim_201( &
      & lambda, &
      & eta, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_H3d_0_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_201
    subroutine cshim_202( &
      & lambda, &
      & eta, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_H3d_1_e")
      import :: c_ptr
      type(c_ptr), value :: lambda
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_202
    subroutine cshim_203( &
      & l, &
      & lambda, &
      & eta, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="legendre_H3d_e")
      import :: c_ptr
      type(c_ptr), value :: l
      type(c_ptr), value :: lambda
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_203
    subroutine cshim_204( &
      & lmax, &
      & lambda, &
      & eta, &
      & len, &
      & out, &
      & status) bind(C, name="legendre_H3d_array")
      import :: c_ptr
      type(c_ptr), value :: lmax
      type(c_ptr), value :: lambda
      type(c_ptr), value :: eta
      type(c_ptr), value :: len
      type(c_ptr), value :: out
      type(c_ptr), value :: status
    end subroutine cshim_204
    subroutine cshim_205( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="log_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_205
    subroutine cshim_206( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="log_abs_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_206
    subroutine cshim_207( &
      & zr, &
      & zi, &
      & len, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status) bind(C, name="complex_log_e")
      import :: c_ptr
      type(c_ptr), value :: zr
      type(c_ptr), value :: zi
      type(c_ptr), value :: len
      type(c_ptr), value :: val_lnr
      type(c_ptr), value :: val_arg
      type(c_ptr), value :: err_lnr
      type(c_ptr), value :: err_arg
      type(c_ptr), value :: status
    end subroutine cshim_207
    subroutine cshim_208( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="log_1plusx_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_208
    subroutine cshim_209( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="log_1plusx_mx_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_209
    subroutine cshim_210( &
      & c, &
      & len, &
      & x, &
      & lenx, &
      & ans) bind(C, name="gsl_poly")
      import :: c_ptr
      type(c_ptr), value :: c
      type(c_ptr), value :: len
      type(c_ptr), value :: x
      type(c_ptr), value :: lenx
      type(c_ptr), value :: ans
    end subroutine cshim_210
    subroutine cshim_211( &
      & x, &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="pow_int")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_211
    subroutine cshim_212( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi_int")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_212
    subroutine cshim_213( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_213
    subroutine cshim_214( &
      & y, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi_1piy")
      import :: c_ptr
      type(c_ptr), value :: y
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_214
    subroutine cshim_215( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi_1_int")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_215
    subroutine cshim_216( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi_1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_216
    subroutine cshim_217( &
      & m, &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="psi_n")
      import :: c_ptr
      type(c_ptr), value :: m
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_217
    subroutine cshim_218( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="synchrotron_1")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_218
    subroutine cshim_219( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="synchrotron_2")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_219
    subroutine cshim_220( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="transport_2")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_220
    subroutine cshim_221( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="transport_3")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_221
    subroutine cshim_222( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="transport_4")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_222
    subroutine cshim_223( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="transport_5")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_223
    subroutine cshim_224( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="sin_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_224
    subroutine cshim_225( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="cos_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_225
    subroutine cshim_226( &
      & x, &
      & y, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hypot_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: y
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_226
    subroutine cshim_227( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="sinc_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_227
    subroutine cshim_228( &
      & zr, &
      & zi, &
      & len, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status) bind(C, name="complex_sin_e")
      import :: c_ptr
      type(c_ptr), value :: zr
      type(c_ptr), value :: zi
      type(c_ptr), value :: len
      type(c_ptr), value :: val_lnr
      type(c_ptr), value :: val_arg
      type(c_ptr), value :: err_lnr
      type(c_ptr), value :: err_arg
      type(c_ptr), value :: status
    end subroutine cshim_228
    subroutine cshim_229( &
      & zr, &
      & zi, &
      & len, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status) bind(C, name="complex_cos_e")
      import :: c_ptr
      type(c_ptr), value :: zr
      type(c_ptr), value :: zi
      type(c_ptr), value :: len
      type(c_ptr), value :: val_lnr
      type(c_ptr), value :: val_arg
      type(c_ptr), value :: err_lnr
      type(c_ptr), value :: err_arg
      type(c_ptr), value :: status
    end subroutine cshim_229
    subroutine cshim_230( &
      & zr, &
      & zi, &
      & len, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status) bind(C, name="complex_logsin_e")
      import :: c_ptr
      type(c_ptr), value :: zr
      type(c_ptr), value :: zi
      type(c_ptr), value :: len
      type(c_ptr), value :: val_lnr
      type(c_ptr), value :: val_arg
      type(c_ptr), value :: err_lnr
      type(c_ptr), value :: err_arg
      type(c_ptr), value :: status
    end subroutine cshim_230
    subroutine cshim_231( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lnsinh_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_231
    subroutine cshim_232( &
      & x, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="lncosh_e")
      import :: c_ptr
      type(c_ptr), value :: x
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_232
    subroutine cshim_233( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="zeta_int")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_233
    subroutine cshim_234( &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="zeta")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_234
    subroutine cshim_235( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="zetam1_int")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_235
    subroutine cshim_236( &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="zetam1")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_236
    subroutine cshim_237( &
      & s, &
      & q, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="hzeta")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: q
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_237
    subroutine cshim_238( &
      & n, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="eta_int")
      import :: c_ptr
      type(c_ptr), value :: n
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_238
    subroutine cshim_239( &
      & s, &
      & len, &
      & val, &
      & err, &
      & status) bind(C, name="eta")
      import :: c_ptr
      type(c_ptr), value :: s
      type(c_ptr), value :: len
      type(c_ptr), value :: val
      type(c_ptr), value :: err
      type(c_ptr), value :: status
    end subroutine cshim_239
  end interface

contains

  subroutine airy_ai( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_001( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_ai

  subroutine airy_bi( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_002( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_bi

  subroutine airy_ai_scaled( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_003( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_ai_scaled

  subroutine airy_bi_scaled( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_004( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_bi_scaled

  subroutine airy_ai_deriv( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_005( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_ai_deriv

  subroutine airy_bi_deriv( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_006( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_bi_deriv

  subroutine airy_ai_deriv_scaled( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_007( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_ai_deriv_scaled

  subroutine airy_bi_deriv_scaled( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_008( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_bi_deriv_scaled

  subroutine airy_zero_ai( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_009( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_zero_ai

  subroutine airy_zero_bi( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_010( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_zero_bi

  subroutine airy_zero_ai_deriv( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_011( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_zero_ai_deriv

  subroutine airy_zero_bi_deriv( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from airy.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_012( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine airy_zero_bi_deriv

  subroutine bessel_cyl_j0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_013( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_j0

  subroutine bessel_cyl_j1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_014( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_j1

  subroutine bessel_cyl_jn( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_015( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_jn

  subroutine bessel_cyl_jn_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_016( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_jn_array

  subroutine bessel_cyl_y0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_017( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_y0

  subroutine bessel_cyl_y1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_018( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_y1

  subroutine bessel_cyl_yn( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_019( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_yn

  subroutine bessel_cyl_yn_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_020( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_yn_array

  subroutine bessel_mod_i0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_021( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_i0

  subroutine bessel_mod_i1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_022( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_i1

  subroutine bessel_mod_in( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_023( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_in

  subroutine bessel_mod_in_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_024( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_in_array

  subroutine bessel_mod_i0_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_025( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_i0_scaled

  subroutine bessel_mod_i1_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_026( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_i1_scaled

  subroutine bessel_mod_in_scaled( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_027( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_in_scaled

  subroutine bessel_mod_in_scaled_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_028( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_in_scaled_array

  subroutine bessel_mod_k0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_029( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_k0

  subroutine bessel_mod_k1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_030( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_k1

  subroutine bessel_mod_kn( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_031( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_kn

  subroutine bessel_mod_kn_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_032( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_kn_array

  subroutine bessel_mod_k0_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_033( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_k0_scaled

  subroutine bessel_mod_k1_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_034( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_k1_scaled

  subroutine bessel_mod_kn_scaled( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_035( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_kn_scaled

  subroutine bessel_mod_kn_scaled_array( &
      & nmin, &
      & nmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: nmin
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_036( &
      & c_loc(nmin), &
      & c_loc(nmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_kn_scaled_array

  subroutine bessel_sph_j0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_037( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_j0

  subroutine bessel_sph_j1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_038( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_j1

  subroutine bessel_sph_j2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_039( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_j2

  subroutine bessel_sph_jl( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_040( &
      & c_loc(l(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_jl

  subroutine bessel_sph_jl_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_041( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_jl_array

  subroutine bessel_sph_jl_steed_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_042( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_jl_steed_array

  subroutine bessel_sph_y0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_043( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_y0

  subroutine bessel_sph_y1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_044( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_y1

  subroutine bessel_sph_y2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_045( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_y2

  subroutine bessel_sph_yl( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_046( &
      & c_loc(l(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_yl

  subroutine bessel_sph_yl_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_047( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_yl_array

  subroutine bessel_sph_i0_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_048( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_i0_scaled

  subroutine bessel_sph_i1_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_049( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_i1_scaled

  subroutine bessel_sph_i2_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_050( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_i2_scaled

  subroutine bessel_sph_il_scaled( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_051( &
      & c_loc(l(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_il_scaled

  subroutine bessel_sph_il_scaled_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_052( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_il_scaled_array

  subroutine bessel_sph_k0_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_053( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_k0_scaled

  subroutine bessel_sph_k1_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_054( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_k1_scaled

  subroutine bessel_sph_k2_scaled( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_055( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_k2_scaled

  subroutine bessel_sph_kl_scaled( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_056( &
      & c_loc(l(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_kl_scaled

  subroutine bessel_sph_kl_scaled_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_057( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine bessel_sph_kl_scaled_array

  subroutine bessel_cyl_jnu( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_058( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_jnu

  subroutine bessel_sequence_jnu( &
      & nu, &
      & v, &
      & nv, &
      & mode, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu
    real(c_double), intent(in), target :: v(:)
    integer(c_int), intent(in), target :: nv
    integer(c_int), intent(in), target :: mode
    integer(c_int), intent(out), target :: status
    call cshim_059( &
      & c_loc(nu), &
      & c_loc(v(1)), &
      & c_loc(nv), &
      & c_loc(mode), &
      & c_loc(status))
  end subroutine bessel_sequence_jnu

  subroutine bessel_cyl_ynu( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_060( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_cyl_ynu

  subroutine bessel_mod_inu( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_061( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_inu

  subroutine bessel_mod_inu_scaled( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_062( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_inu_scaled

  subroutine bessel_mod_knu( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_063( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_knu

  subroutine bessel_lnknu( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_064( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_lnknu

  subroutine bessel_mod_knu_scaled( &
      & nu, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_065( &
      & c_loc(nu(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_mod_knu_scaled

  subroutine bessel_zero_j0( &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_066( &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_zero_j0

  subroutine bessel_zero_j1( &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    integer(c_int), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_067( &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_zero_j1

  subroutine bessel_zero_jnu( &
      & nu, &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from bessel.c.
    real(c_double), intent(in), target :: nu(:)
    integer(c_int), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(nu), c_int)
    call cshim_068( &
      & c_loc(nu(1)), &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine bessel_zero_jnu

  subroutine clausen( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from clausen.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_069( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine clausen

  subroutine hydrogenicr_1( &
      & Z, &
      & r, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: Z(:)
    real(c_double), intent(in), target :: r(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(Z), c_int)
    call cshim_070( &
      & c_loc(Z(1)), &
      & c_loc(r(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hydrogenicr_1

  subroutine hydrogenicr( &
      & n, &
      & l, &
      & Z, &
      & r, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    integer(c_int), intent(in), target :: n(:)
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: Z(:)
    real(c_double), intent(in), target :: r(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_071( &
      & c_loc(n(1)), &
      & c_loc(l(1)), &
      & c_loc(Z(1)), &
      & c_loc(r(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hydrogenicr

  subroutine coulomb_wave_fg( &
      & eta, &
      & x, &
      & L_F, &
      & k, &
      & val_F, &
      & err_F, &
      & val_Fp, &
      & err_Fp, &
      & val_G, &
      & err_G, &
      & val_Gp, &
      & err_Gp, &
      & exp_F, &
      & exp_G, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: L_F(:)
    integer(c_int), intent(in), target :: k(:)
    real(c_double), intent(in), target :: val_F(:)
    real(c_double), intent(in), target :: err_F(:)
    real(c_double), intent(in), target :: val_Fp(:)
    real(c_double), intent(in), target :: err_Fp(:)
    real(c_double), intent(in), target :: val_G(:)
    real(c_double), intent(in), target :: err_G(:)
    real(c_double), intent(in), target :: val_Gp(:)
    real(c_double), intent(in), target :: err_Gp(:)
    real(c_double), intent(in), target :: exp_F
    real(c_double), intent(in), target :: exp_G
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_072( &
      & c_loc(eta(1)), &
      & c_loc(x(1)), &
      & c_loc(L_F(1)), &
      & c_loc(k(1)), &
      & c_loc(n_local), &
      & c_loc(val_F(1)), &
      & c_loc(err_F(1)), &
      & c_loc(val_Fp(1)), &
      & c_loc(err_Fp(1)), &
      & c_loc(val_G(1)), &
      & c_loc(err_G(1)), &
      & c_loc(val_Gp(1)), &
      & c_loc(err_Gp(1)), &
      & c_loc(exp_F), &
      & c_loc(exp_G), &
      & c_loc(status(1)))
  end subroutine coulomb_wave_fg

  subroutine coulomb_wave_f_array( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & fc_array, &
      & F_exponent, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L_min
    integer(c_int), intent(in), target :: kmax
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: fc_array(:)
    real(c_double), intent(out), target :: F_exponent(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_073( &
      & c_loc(L_min), &
      & c_loc(kmax), &
      & c_loc(eta(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(fc_array(1)), &
      & c_loc(F_exponent(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_wave_f_array

  subroutine coulomb_wave_fg_array( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & fc_array, &
      & gc_array, &
      & F_exponent, &
      & G_exponent, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L_min
    integer(c_int), intent(in), target :: kmax
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: fc_array(:)
    real(c_double), intent(out), target :: gc_array(:)
    real(c_double), intent(out), target :: F_exponent(:)
    real(c_double), intent(out), target :: G_exponent(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_074( &
      & c_loc(L_min), &
      & c_loc(kmax), &
      & c_loc(eta(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(fc_array(1)), &
      & c_loc(gc_array(1)), &
      & c_loc(F_exponent(1)), &
      & c_loc(G_exponent(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_wave_fg_array

  subroutine coulomb_wave_fgp_array( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & fc_array, &
      & fcp_array, &
      & gc_array, &
      & gcp_array, &
      & F_exponent, &
      & G_exponent, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L_min
    integer(c_int), intent(in), target :: kmax
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: fc_array(:)
    real(c_double), intent(out), target :: fcp_array(:)
    real(c_double), intent(out), target :: gc_array(:)
    real(c_double), intent(out), target :: gcp_array(:)
    real(c_double), intent(out), target :: F_exponent(:)
    real(c_double), intent(out), target :: G_exponent(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_075( &
      & c_loc(L_min), &
      & c_loc(kmax), &
      & c_loc(eta(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(fc_array(1)), &
      & c_loc(fcp_array(1)), &
      & c_loc(gc_array(1)), &
      & c_loc(gcp_array(1)), &
      & c_loc(F_exponent(1)), &
      & c_loc(G_exponent(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_wave_fgp_array

  subroutine coulomb_wave_sphf_array( &
      & L_min, &
      & kmax, &
      & eta, &
      & x, &
      & fc_array, &
      & F_exponent, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L_min
    integer(c_int), intent(in), target :: kmax
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: fc_array(:)
    real(c_double), intent(out), target :: F_exponent(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_076( &
      & c_loc(L_min), &
      & c_loc(kmax), &
      & c_loc(eta(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(fc_array(1)), &
      & c_loc(F_exponent(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_wave_sphf_array

  subroutine coulomb_cl( &
      & L, &
      & eta, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L(:)
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(L), c_int)
    call cshim_077( &
      & c_loc(L(1)), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_cl

  subroutine coulomb_cl_array( &
      & L_min, &
      & kmax, &
      & eta, &
      & cl, &
      & status)
    ! Vector-oriented wrapper translated from coulomb.c.
    real(c_double), intent(in), target :: L_min
    integer(c_int), intent(in), target :: kmax
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: cl(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(eta), c_int)
    call cshim_078( &
      & c_loc(L_min), &
      & c_loc(kmax), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(cl(1)), &
      & c_loc(status(1)))
  end subroutine coulomb_cl_array

  subroutine coupling_3j( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_ma, &
      & two_mb, &
      & two_mc, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coupling.c.
    integer(c_int), intent(in), target :: two_ja(:)
    integer(c_int), intent(in), target :: two_jb(:)
    integer(c_int), intent(in), target :: two_jc(:)
    integer(c_int), intent(in), target :: two_ma(:)
    integer(c_int), intent(in), target :: two_mb(:)
    integer(c_int), intent(in), target :: two_mc(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(two_ja), c_int)
    call cshim_079( &
      & c_loc(two_ja(1)), &
      & c_loc(two_jb(1)), &
      & c_loc(two_jc(1)), &
      & c_loc(two_ma(1)), &
      & c_loc(two_mb(1)), &
      & c_loc(two_mc(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine coupling_3j

  subroutine coupling_6j( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_jd, &
      & two_je, &
      & two_jf, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coupling.c.
    integer(c_int), intent(in), target :: two_ja(:)
    integer(c_int), intent(in), target :: two_jb(:)
    integer(c_int), intent(in), target :: two_jc(:)
    integer(c_int), intent(in), target :: two_jd(:)
    integer(c_int), intent(in), target :: two_je(:)
    integer(c_int), intent(in), target :: two_jf(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(two_ja), c_int)
    call cshim_080( &
      & c_loc(two_ja(1)), &
      & c_loc(two_jb(1)), &
      & c_loc(two_jc(1)), &
      & c_loc(two_jd(1)), &
      & c_loc(two_je(1)), &
      & c_loc(two_jf(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine coupling_6j

  subroutine coupling_9j( &
      & two_ja, &
      & two_jb, &
      & two_jc, &
      & two_jd, &
      & two_je, &
      & two_jf, &
      & two_jg, &
      & two_jh, &
      & two_ji, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from coupling.c.
    integer(c_int), intent(in), target :: two_ja(:)
    integer(c_int), intent(in), target :: two_jb(:)
    integer(c_int), intent(in), target :: two_jc(:)
    integer(c_int), intent(in), target :: two_jd(:)
    integer(c_int), intent(in), target :: two_je(:)
    integer(c_int), intent(in), target :: two_jf(:)
    integer(c_int), intent(in), target :: two_jg(:)
    integer(c_int), intent(in), target :: two_jh(:)
    integer(c_int), intent(in), target :: two_ji(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(two_ja), c_int)
    call cshim_081( &
      & c_loc(two_ja(1)), &
      & c_loc(two_jb(1)), &
      & c_loc(two_jc(1)), &
      & c_loc(two_jd(1)), &
      & c_loc(two_je(1)), &
      & c_loc(two_jf(1)), &
      & c_loc(two_jg(1)), &
      & c_loc(two_jh(1)), &
      & c_loc(two_ji(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine coupling_9j

  subroutine dawson( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from dawson.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_082( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine dawson

  subroutine debye_1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from debye.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_083( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine debye_1

  subroutine debye_2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from debye.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_084( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine debye_2

  subroutine debye_3( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from debye.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_085( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine debye_3

  subroutine debye_4( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from debye.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_086( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine debye_4

  subroutine dilog( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from dilog.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_087( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine dilog

  subroutine complex_dilog( &
      & r, &
      & theta, &
      & val_re, &
      & val_im, &
      & err_re, &
      & err_im, &
      & status)
    ! Vector-oriented wrapper translated from dilog.c.
    real(c_double), intent(in), target :: r(:)
    real(c_double), intent(in), target :: theta(:)
    real(c_double), intent(in), target :: val_re(:)
    real(c_double), intent(in), target :: val_im(:)
    real(c_double), intent(in), target :: err_re(:)
    real(c_double), intent(in), target :: err_im(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(r), c_int)
    call cshim_088( &
      & c_loc(r(1)), &
      & c_loc(theta(1)), &
      & c_loc(n_local), &
      & c_loc(val_re(1)), &
      & c_loc(val_im(1)), &
      & c_loc(err_re(1)), &
      & c_loc(err_im(1)), &
      & c_loc(status(1)))
  end subroutine complex_dilog

  subroutine ellint_kcomp( &
      & k, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: k(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(k), c_int)
    call cshim_089( &
      & c_loc(k(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_kcomp

  subroutine ellint_ecomp( &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: k(:)
    integer(c_int), intent(in), target :: nk
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_090( &
      & c_loc(k(1)), &
      & c_loc(nk), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_ecomp

  subroutine ellint_f( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: phi(:)
    real(c_double), intent(in), target :: k(:)
    integer(c_int), intent(in), target :: nk
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_091( &
      & c_loc(phi(1)), &
      & c_loc(k(1)), &
      & c_loc(nk), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_f

  subroutine ellint_e( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: phi(:)
    real(c_double), intent(in), target :: k(:)
    integer(c_int), intent(in), target :: nk
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_092( &
      & c_loc(phi(1)), &
      & c_loc(k(1)), &
      & c_loc(nk), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_e

  subroutine ellint_p( &
      & phi, &
      & k, &
      & n, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: phi(:)
    real(c_double), intent(in), target :: k(:)
    real(c_double), intent(in), target :: n(:)
    integer(c_int), intent(in), target :: nk
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_093( &
      & c_loc(phi(1)), &
      & c_loc(k(1)), &
      & c_loc(n(1)), &
      & c_loc(nk), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_p

  subroutine ellint_d( &
      & phi, &
      & k, &
      & nk, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: phi(:)
    real(c_double), intent(in), target :: k(:)
    integer(c_int), intent(in), target :: nk
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_094( &
      & c_loc(phi(1)), &
      & c_loc(k(1)), &
      & c_loc(nk), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_d

  subroutine ellint_rc( &
      & x, &
      & y, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: y(:)
    integer(c_int), intent(in), target :: nx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_095( &
      & c_loc(x(1)), &
      & c_loc(y(1)), &
      & c_loc(nx), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_rc

  subroutine ellint_rd( &
      & x, &
      & y, &
      & z, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: y(:)
    real(c_double), intent(in), target :: z(:)
    integer(c_int), intent(in), target :: nx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_096( &
      & c_loc(x(1)), &
      & c_loc(y(1)), &
      & c_loc(z(1)), &
      & c_loc(nx), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_rd

  subroutine ellint_rf( &
      & x, &
      & y, &
      & z, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: y(:)
    real(c_double), intent(in), target :: z(:)
    integer(c_int), intent(in), target :: nx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_097( &
      & c_loc(x(1)), &
      & c_loc(y(1)), &
      & c_loc(z(1)), &
      & c_loc(nx), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_rf

  subroutine ellint_rj( &
      & x, &
      & y, &
      & z, &
      & p, &
      & nx, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from ellint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: y(:)
    real(c_double), intent(in), target :: z(:)
    real(c_double), intent(in), target :: p
    integer(c_int), intent(in), target :: nx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_098( &
      & c_loc(x(1)), &
      & c_loc(y(1)), &
      & c_loc(z(1)), &
      & c_loc(p), &
      & c_loc(nx), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ellint_rj

  subroutine elljac( &
      & u, &
      & m, &
      & sn, &
      & cn, &
      & dn, &
      & status)
    ! Vector-oriented wrapper translated from elljac.c.
    real(c_double), intent(in), target :: u(:)
    real(c_double), intent(in), target :: m(:)
    real(c_double), intent(out), target :: sn(:)
    real(c_double), intent(out), target :: cn(:)
    real(c_double), intent(out), target :: dn(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(u), c_int)
    call cshim_099( &
      & c_loc(u(1)), &
      & c_loc(m(1)), &
      & c_loc(n_local), &
      & c_loc(sn(1)), &
      & c_loc(cn(1)), &
      & c_loc(dn(1)), &
      & c_loc(status(1)))
  end subroutine elljac

  subroutine erf( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_100( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine erf

  subroutine erfc( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_101( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine erfc

  subroutine log_erfc( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_102( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine log_erfc

  subroutine erf_z( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_103( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine erf_z

  subroutine erf_q( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_104( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine erf_q

  subroutine hazard( &
      & x, &
      & mode, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from error.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_105( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(mode), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hazard

  subroutine expint_e1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_106( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine expint_e1

  subroutine expint_e2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_107( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine expint_e2

  subroutine expint_en( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_108( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine expint_en

  subroutine expint_ei( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_109( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine expint_ei

  subroutine shi( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_110( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine shi

  subroutine chi( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_111( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine chi

  subroutine expint_3( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_112( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine expint_3

  subroutine si( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_113( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine si

  subroutine ci( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_114( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine ci

  subroutine atanint( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from expint.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_115( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine atanint

  subroutine fermi_dirac_m1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_116( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_m1

  subroutine fermi_dirac_0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_117( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_0

  subroutine fermi_dirac_1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_118( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_1

  subroutine fermi_dirac_2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_119( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_2

  subroutine fermi_dirac_int( &
      & j, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    integer(c_int), intent(in), target :: j(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(j), c_int)
    call cshim_120( &
      & c_loc(j(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_int

  subroutine fermi_dirac_mhalf( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_121( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_mhalf

  subroutine fermi_dirac_half( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_122( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_half

  subroutine fermi_dirac_3half( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_123( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_3half

  subroutine fermi_dirac_inc_0( &
      & x, &
      & b, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from fermi_dirac.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_124( &
      & c_loc(x(1)), &
      & c_loc(b(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fermi_dirac_inc_0

  subroutine gsl_sf_gamma( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_125( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_gamma

  subroutine lngamma( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_126( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lngamma

  subroutine lngamma_sgn( &
      & x, &
      & val, &
      & err, &
      & status, &
      & sgn)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    real(c_double), intent(out), target :: sgn(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_127( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)), &
      & c_loc(sgn(1)))
  end subroutine lngamma_sgn

  subroutine gammastar( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_128( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gammastar

  subroutine gammainv( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_129( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gammainv

  subroutine lngamma_complex( &
      & zr, &
      & zi, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: zr(:)
    real(c_double), intent(in), target :: zi(:)
    real(c_double), intent(in), target :: val_lnr(:)
    real(c_double), intent(in), target :: val_arg(:)
    real(c_double), intent(in), target :: err_lnr(:)
    real(c_double), intent(in), target :: err_arg(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(zr), c_int)
    call cshim_130( &
      & c_loc(zr(1)), &
      & c_loc(zi(1)), &
      & c_loc(n_local), &
      & c_loc(val_lnr(1)), &
      & c_loc(val_arg(1)), &
      & c_loc(err_lnr(1)), &
      & c_loc(err_arg(1)), &
      & c_loc(status(1)))
  end subroutine lngamma_complex

  subroutine taylorcoeff( &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_131( &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine taylorcoeff

  subroutine fact( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_132( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine fact

  subroutine doublefact( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_133( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine doublefact

  subroutine lnfact( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_134( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lnfact

  subroutine lndoublefact( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_135( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lndoublefact

  subroutine gsl_sf_choose( &
      & n, &
      & m, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    integer(c_int), intent(in), target :: m(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_136( &
      & c_loc(n(1)), &
      & c_loc(m(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_choose

  subroutine lnchoose( &
      & n, &
      & m, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    integer(c_int), intent(in), target :: n(:)
    integer(c_int), intent(in), target :: m(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_137( &
      & c_loc(n(1)), &
      & c_loc(m(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lnchoose

  subroutine poch( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_138( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine poch

  subroutine lnpoch( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_139( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lnpoch

  subroutine lnpoch_sgn( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status, &
      & sgn)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    real(c_double), intent(out), target :: sgn(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_140( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)), &
      & c_loc(sgn(1)))
  end subroutine lnpoch_sgn

  subroutine pochrel( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_141( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine pochrel

  subroutine gamma_inc_p( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_142( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gamma_inc_p

  subroutine gamma_inc_q( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_143( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gamma_inc_q

  subroutine gamma_inc( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_144( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gamma_inc

  subroutine gsl_sf_beta( &
      & a, &
      & b, &
      & nb, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    integer(c_int), intent(in), target :: nb
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    call cshim_145( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(nb), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_beta

  subroutine lnbeta( &
      & a, &
      & b, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_146( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lnbeta

  subroutine beta_inc( &
      & a, &
      & b, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gamma.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_147( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine beta_inc

  subroutine gegenpoly_1( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gegenbauer.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_148( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gegenpoly_1

  subroutine gegenpoly_2( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gegenbauer.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_149( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gegenpoly_2

  subroutine gegenpoly_3( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gegenbauer.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_150( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gegenpoly_3

  subroutine gegenpoly_n( &
      & n, &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from gegenbauer.c.
    integer(c_int), intent(in), target :: n
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_151( &
      & c_loc(n), &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gegenpoly_n

  subroutine gegenpoly_array( &
      & nmax, &
      & lambda, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from gegenbauer.c.
    integer(c_int), intent(in), target :: nmax
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_152( &
      & c_loc(nmax), &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine gegenpoly_array

  subroutine hyperg_0f1( &
      & c, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: c(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(c), c_int)
    call cshim_153( &
      & c_loc(c(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_0f1

  subroutine hyperg_1f1_int( &
      & m, &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    integer(c_int), intent(in), target :: m(:)
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(m), c_int)
    call cshim_154( &
      & c_loc(m(1)), &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_1f1_int

  subroutine hyperg_1f1( &
      & a, &
      & b, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_155( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_1f1

  subroutine hyperg_u_int( &
      & m, &
      & n, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    integer(c_int), intent(in), target :: m(:)
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(m), c_int)
    call cshim_156( &
      & c_loc(m(1)), &
      & c_loc(n(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_u_int

  subroutine hyperg_u( &
      & a, &
      & b, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_157( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_u

  subroutine hyperg_2f1( &
      & a, &
      & b, &
      & c, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: c(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_158( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(c(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_2f1

  subroutine hyperg_2f1_conj( &
      & aR, &
      & aI, &
      & c, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: aR(:)
    real(c_double), intent(in), target :: aI(:)
    real(c_double), intent(in), target :: c(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(aR), c_int)
    call cshim_159( &
      & c_loc(aR(1)), &
      & c_loc(aI(1)), &
      & c_loc(c(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_2f1_conj

  subroutine hyperg_2f1_renorm( &
      & a, &
      & b, &
      & c, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: c(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_160( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(c(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_2f1_renorm

  subroutine hyperg_2f1_conj_renorm( &
      & aR, &
      & aI, &
      & c, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: aR(:)
    real(c_double), intent(in), target :: aI(:)
    real(c_double), intent(in), target :: c(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(aR), c_int)
    call cshim_161( &
      & c_loc(aR(1)), &
      & c_loc(aI(1)), &
      & c_loc(c(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_2f1_conj_renorm

  subroutine hyperg_2f0( &
      & a, &
      & b, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from hyperg.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: b(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_162( &
      & c_loc(a(1)), &
      & c_loc(b(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hyperg_2f0

  subroutine laguerre_1( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from laguerre.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_163( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine laguerre_1

  subroutine laguerre_2( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from laguerre.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_164( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine laguerre_2

  subroutine laguerre_3( &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from laguerre.c.
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_165( &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine laguerre_3

  subroutine laguerre_n( &
      & n, &
      & a, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from laguerre.c.
    integer(c_int), intent(in), target :: n
    real(c_double), intent(in), target :: a(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(a), c_int)
    call cshim_166( &
      & c_loc(n), &
      & c_loc(a(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine laguerre_n

  subroutine lambert_w0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from lambert.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_167( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lambert_w0

  subroutine lambert_wm1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from lambert.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_168( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lambert_wm1

  subroutine legendre_p1_raw( &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_169( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_p1_raw

  subroutine legendre_p1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_170( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_p1

  subroutine legendre_p2_raw( &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_171( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_p2_raw

  subroutine legendre_p2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_172( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_p2

  subroutine legendre_p3_raw( &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_173( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_p3_raw

  subroutine legendre_p3( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_174( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_p3

  subroutine legendre_pl_raw( &
      & l, &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_175( &
      & c_loc(l), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_pl_raw

  subroutine legendre_pl( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_176( &
      & c_loc(l), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_pl

  subroutine legendre_pl_array( &
      & lmax, &
      & x, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_177( &
      & c_loc(lmax), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine legendre_pl_array

  subroutine legendre_q0_raw( &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_178( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_q0_raw

  subroutine legendre_q0( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_179( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_q0

  subroutine legendre_q1_raw( &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_180( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_q1_raw

  subroutine legendre_q1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_181( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_q1

  subroutine legendre_ql_raw( &
      & l, &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_182( &
      & c_loc(l), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_ql_raw

  subroutine legendre_ql( &
      & l, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_183( &
      & c_loc(l), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_ql

  subroutine legendre_array_n( &
      & lmax, &
      & ans)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(out), target :: ans
    call cshim_184( &
      & c_loc(lmax), &
      & c_loc(ans))
  end subroutine legendre_array_n

  subroutine legendre_array_index( &
      & l, &
      & m, &
      & ans)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    integer(c_int), intent(in), target :: m
    integer(c_int), intent(out), target :: ans
    call cshim_185( &
      & c_loc(l), &
      & c_loc(m), &
      & c_loc(ans))
  end subroutine legendre_array_index

  subroutine legendre_array( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(in), target :: norm
    real(c_double), intent(in), target :: csphase
    real(c_double), intent(out), target :: result_array(:)
    real(c_double), intent(out), target :: out(:)
    call cshim_186( &
      & c_loc(x), &
      & c_loc(lmax), &
      & c_loc(norm), &
      & c_loc(csphase), &
      & c_loc(result_array(1)), &
      & c_loc(out(1)))
  end subroutine legendre_array

  subroutine legendre_deriv_array( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(in), target :: norm
    real(c_double), intent(in), target :: csphase
    real(c_double), intent(out), target :: result_array(:)
    real(c_double), intent(out), target :: result_deriv_array(:)
    real(c_double), intent(out), target :: out(:)
    call cshim_187( &
      & c_loc(x), &
      & c_loc(lmax), &
      & c_loc(norm), &
      & c_loc(csphase), &
      & c_loc(result_array(1)), &
      & c_loc(result_deriv_array(1)), &
      & c_loc(out(1)))
  end subroutine legendre_deriv_array

  subroutine legendre_deriv_alt_array( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(in), target :: norm
    real(c_double), intent(in), target :: csphase
    real(c_double), intent(out), target :: result_array(:)
    real(c_double), intent(out), target :: result_deriv_array(:)
    real(c_double), intent(out), target :: out(:)
    call cshim_188( &
      & c_loc(x), &
      & c_loc(lmax), &
      & c_loc(norm), &
      & c_loc(csphase), &
      & c_loc(result_array(1)), &
      & c_loc(result_deriv_array(1)), &
      & c_loc(out(1)))
  end subroutine legendre_deriv_alt_array

  subroutine legendre_deriv2_array( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & result_deriv2_array, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(in), target :: norm
    real(c_double), intent(in), target :: csphase
    real(c_double), intent(out), target :: result_array(:)
    real(c_double), intent(out), target :: result_deriv_array(:)
    real(c_double), intent(out), target :: result_deriv2_array(:)
    real(c_double), intent(out), target :: out(:)
    call cshim_189( &
      & c_loc(x), &
      & c_loc(lmax), &
      & c_loc(norm), &
      & c_loc(csphase), &
      & c_loc(result_array(1)), &
      & c_loc(result_deriv_array(1)), &
      & c_loc(result_deriv2_array(1)), &
      & c_loc(out(1)))
  end subroutine legendre_deriv2_array

  subroutine legendre_deriv2_alt_array( &
      & x, &
      & lmax, &
      & norm, &
      & csphase, &
      & result_array, &
      & result_deriv_array, &
      & result_deriv2_array, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: x
    integer(c_int), intent(in), target :: lmax
    integer(c_int), intent(in), target :: norm
    real(c_double), intent(in), target :: csphase
    real(c_double), intent(out), target :: result_array(:)
    real(c_double), intent(out), target :: result_deriv_array(:)
    real(c_double), intent(out), target :: result_deriv2_array(:)
    real(c_double), intent(out), target :: out(:)
    call cshim_190( &
      & c_loc(x), &
      & c_loc(lmax), &
      & c_loc(norm), &
      & c_loc(csphase), &
      & c_loc(result_array(1)), &
      & c_loc(result_deriv_array(1)), &
      & c_loc(result_deriv2_array(1)), &
      & c_loc(out(1)))
  end subroutine legendre_deriv2_alt_array

  subroutine legendre_plm_raw( &
      & l, &
      & m, &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    integer(c_int), intent(in), target :: m
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_191( &
      & c_loc(l), &
      & c_loc(m), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_plm_raw

  subroutine legendre_plm( &
      & l, &
      & m, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    integer(c_int), intent(in), target :: m
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_192( &
      & c_loc(l), &
      & c_loc(m), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_plm

  subroutine legendre_sphplm_raw( &
      & l, &
      & m, &
      & x, &
      & out)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    integer(c_int), intent(in), target :: m
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_193( &
      & c_loc(l), &
      & c_loc(m), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)))
  end subroutine legendre_sphplm_raw

  subroutine legendre_sphplm( &
      & l, &
      & m, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l
    integer(c_int), intent(in), target :: m
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_194( &
      & c_loc(l), &
      & c_loc(m), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_sphplm

  subroutine conicalp_half( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_195( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_half

  subroutine conicalp_mhalf( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_196( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_mhalf

  subroutine conicalp_0( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_197( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_0

  subroutine conicalp_1( &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_198( &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_1

  subroutine conicalp_sph_reg( &
      & l, &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_199( &
      & c_loc(l(1)), &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_sph_reg

  subroutine conicalp_cyl_reg( &
      & m, &
      & lambda, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: m(:)
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(m), c_int)
    call cshim_200( &
      & c_loc(m(1)), &
      & c_loc(lambda(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine conicalp_cyl_reg

  subroutine legendre_h3d_0( &
      & lambda, &
      & eta, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_201( &
      & c_loc(lambda(1)), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_h3d_0

  subroutine legendre_h3d_1( &
      & lambda, &
      & eta, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_202( &
      & c_loc(lambda(1)), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_h3d_1

  subroutine legendre_h3d( &
      & l, &
      & lambda, &
      & eta, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: l(:)
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(l), c_int)
    call cshim_203( &
      & c_loc(l(1)), &
      & c_loc(lambda(1)), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine legendre_h3d

  subroutine legendre_h3d_array( &
      & lmax, &
      & lambda, &
      & eta, &
      & out, &
      & status)
    ! Vector-oriented wrapper translated from legendre.c.
    integer(c_int), intent(in), target :: lmax
    real(c_double), intent(in), target :: lambda(:)
    real(c_double), intent(in), target :: eta(:)
    real(c_double), intent(out), target :: out(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(lambda), c_int)
    call cshim_204( &
      & c_loc(lmax), &
      & c_loc(lambda(1)), &
      & c_loc(eta(1)), &
      & c_loc(n_local), &
      & c_loc(out(1)), &
      & c_loc(status(1)))
  end subroutine legendre_h3d_array

  subroutine gsl_sf_log( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from log.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_205( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_log

  subroutine log_abs( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from log.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_206( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine log_abs

  subroutine complex_log( &
      & zr, &
      & zi, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status)
    ! Vector-oriented wrapper translated from log.c.
    real(c_double), intent(in), target :: zr(:)
    real(c_double), intent(in), target :: zi(:)
    real(c_double), intent(in), target :: val_lnr(:)
    real(c_double), intent(in), target :: val_arg(:)
    real(c_double), intent(in), target :: err_lnr(:)
    real(c_double), intent(in), target :: err_arg(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(zr), c_int)
    call cshim_207( &
      & c_loc(zr(1)), &
      & c_loc(zi(1)), &
      & c_loc(n_local), &
      & c_loc(val_lnr(1)), &
      & c_loc(val_arg(1)), &
      & c_loc(err_lnr(1)), &
      & c_loc(err_arg(1)), &
      & c_loc(status(1)))
  end subroutine complex_log

  subroutine log_1plusx( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from log.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_208( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine log_1plusx

  subroutine log_1plusx_mx( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from log.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_209( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine log_1plusx_mx

  subroutine gsl_poly_c( &
      & c, &
      & x, &
      & lenx, &
      & ans)
    ! Vector-oriented wrapper translated from poly.c.
    real(c_double), intent(in), target :: c
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: lenx
    real(c_double), intent(out), target :: ans(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_210( &
      & c_loc(c), &
      & c_loc(n_local), &
      & c_loc(x(1)), &
      & c_loc(lenx), &
      & c_loc(ans(1)))
  end subroutine gsl_poly_c

  subroutine pow_int( &
      & x, &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from pow_int.c.
    real(c_double), intent(in), target :: x(:)
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_211( &
      & c_loc(x(1)), &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine pow_int

  subroutine psi_int( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_212( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi_int

  subroutine psi( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_213( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi

  subroutine psi_1piy( &
      & y, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    real(c_double), intent(in), target :: y(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(y), c_int)
    call cshim_214( &
      & c_loc(y(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi_1piy

  subroutine psi_1_int( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_215( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi_1_int

  subroutine psi_1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_216( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi_1

  subroutine psi_n( &
      & m, &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from psi.c.
    integer(c_int), intent(in), target :: m(:)
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(m), c_int)
    call cshim_217( &
      & c_loc(m(1)), &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine psi_n

  subroutine synchrotron_1( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from synchrotron.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_218( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine synchrotron_1

  subroutine synchrotron_2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from synchrotron.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_219( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine synchrotron_2

  subroutine transport_2( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from transport.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_220( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine transport_2

  subroutine transport_3( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from transport.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_221( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine transport_3

  subroutine transport_4( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from transport.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_222( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine transport_4

  subroutine transport_5( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from transport.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_223( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine transport_5

  subroutine gsl_sf_sin( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_224( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_sin

  subroutine gsl_sf_cos( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_225( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine gsl_sf_cos

  subroutine hypot( &
      & x, &
      & y, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(in), target :: y(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_226( &
      & c_loc(x(1)), &
      & c_loc(y(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hypot

  subroutine sinc( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_227( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine sinc

  subroutine complex_sin( &
      & zr, &
      & zi, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: zr(:)
    real(c_double), intent(in), target :: zi(:)
    real(c_double), intent(in), target :: val_lnr(:)
    real(c_double), intent(in), target :: val_arg(:)
    real(c_double), intent(in), target :: err_lnr(:)
    real(c_double), intent(in), target :: err_arg(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(zr), c_int)
    call cshim_228( &
      & c_loc(zr(1)), &
      & c_loc(zi(1)), &
      & c_loc(n_local), &
      & c_loc(val_lnr(1)), &
      & c_loc(val_arg(1)), &
      & c_loc(err_lnr(1)), &
      & c_loc(err_arg(1)), &
      & c_loc(status(1)))
  end subroutine complex_sin

  subroutine complex_cos( &
      & zr, &
      & zi, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: zr(:)
    real(c_double), intent(in), target :: zi(:)
    real(c_double), intent(in), target :: val_lnr(:)
    real(c_double), intent(in), target :: val_arg(:)
    real(c_double), intent(in), target :: err_lnr(:)
    real(c_double), intent(in), target :: err_arg(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(zr), c_int)
    call cshim_229( &
      & c_loc(zr(1)), &
      & c_loc(zi(1)), &
      & c_loc(n_local), &
      & c_loc(val_lnr(1)), &
      & c_loc(val_arg(1)), &
      & c_loc(err_lnr(1)), &
      & c_loc(err_arg(1)), &
      & c_loc(status(1)))
  end subroutine complex_cos

  subroutine complex_logsin( &
      & zr, &
      & zi, &
      & val_lnr, &
      & val_arg, &
      & err_lnr, &
      & err_arg, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: zr(:)
    real(c_double), intent(in), target :: zi(:)
    real(c_double), intent(in), target :: val_lnr(:)
    real(c_double), intent(in), target :: val_arg(:)
    real(c_double), intent(in), target :: err_lnr(:)
    real(c_double), intent(in), target :: err_arg(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(zr), c_int)
    call cshim_230( &
      & c_loc(zr(1)), &
      & c_loc(zi(1)), &
      & c_loc(n_local), &
      & c_loc(val_lnr(1)), &
      & c_loc(val_arg(1)), &
      & c_loc(err_lnr(1)), &
      & c_loc(err_arg(1)), &
      & c_loc(status(1)))
  end subroutine complex_logsin

  subroutine lnsinh( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_231( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lnsinh

  subroutine lncosh( &
      & x, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from trig.c.
    real(c_double), intent(in), target :: x(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(x), c_int)
    call cshim_232( &
      & c_loc(x(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine lncosh

  subroutine zeta_int( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_233( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine zeta_int

  subroutine zeta( &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    real(c_double), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_234( &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine zeta

  subroutine zetam1_int( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_235( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine zetam1_int

  subroutine zetam1( &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    real(c_double), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_236( &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine zetam1

  subroutine hzeta( &
      & s, &
      & q, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    real(c_double), intent(in), target :: s(:)
    real(c_double), intent(in), target :: q(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_237( &
      & c_loc(s(1)), &
      & c_loc(q(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine hzeta

  subroutine eta_int( &
      & n, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    integer(c_int), intent(in), target :: n(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(n), c_int)
    call cshim_238( &
      & c_loc(n(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine eta_int

  subroutine eta( &
      & s, &
      & val, &
      & err, &
      & status)
    ! Vector-oriented wrapper translated from zeta.c.
    real(c_double), intent(in), target :: s(:)
    real(c_double), intent(out), target :: val(:)
    real(c_double), intent(out), target :: err(:)
    integer(c_int), intent(out), target :: status(:)
    integer(c_int), target :: n_local
    n_local = int(size(s), c_int)
    call cshim_239( &
      & c_loc(s(1)), &
      & c_loc(n_local), &
      & c_loc(val(1)), &
      & c_loc(err(1)), &
      & c_loc(status(1)))
  end subroutine eta

end module gsl_special
