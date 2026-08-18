! SPDX-License-Identifier: GPL-2.0-or-later
module desolve
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, ode_result, complex_ode_rhs, complex_ode_result, &
      ode_root, lsodar_result, dae_residual
  use desolve_odepack, only : lsoda, lsode, vode
  use desolve_roots_sparse, only : lsodes, lsodar
  use desolve_stiff, only : radau, daspk, zvode
  use desolve_rk, only : rk_method, rk_method_by_name, rk_integrate, euler, rk4
  use desolve_utilities
  use desolve_dde, only : dde_rhs, dede_rk4
  implicit none
  private
  public :: dp, ode_rhs, ode_result, complex_ode_rhs, complex_ode_result, ode_root, lsodar_result, dae_residual
  public :: lsoda,lsode,lsodes,lsodar,vode,radau,daspk,zvode
  public :: rk_method,rk_method_by_name,rk_integrate,euler,rk4
  public :: forcing_table,event_table,history_buffer,sparsity_pattern,brent_root,iterate_map,sparsity_2d,sparsity_3d
  public :: hermite_value,hermite_deriv,nearest_event,clean_event_times,dde_rhs,dede_rk4
  public :: ode
contains
  function ode(rhs,y0,times,method,rtol,atol,h) result(sol)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::y0(:),times(:)
    character(len=*),intent(in),optional::method
    real(dp),intent(in),optional::rtol,atol,h
    type(ode_result)::sol
    character(len=:),allocatable::m
    type(rk_method)::rkm
    m='lsoda';if(present(method))m=lower(trim(method))
    select case(m)
    case('lsoda')
      sol=lsoda(rhs,y0,times,rtol=opt_real(rtol,1e-8_dp),atol=opt_real(atol,1e-10_dp))
    case('lsode')
      sol=lsode(rhs,y0,times,rtol=opt_real(rtol,1e-8_dp),atol=opt_real(atol,1e-10_dp))
    case('lsodes')
      sol=lsodes(rhs,y0,times,rtol=opt_real(rtol,1e-8_dp),atol=opt_real(atol,1e-10_dp))
    case('vode')
      sol=vode(rhs,y0,times,rtol=opt_real(rtol,1e-8_dp),atol=opt_real(atol,1e-10_dp))
    case('radau','radau5')
      sol=radau(rhs,y0,times,rtol=opt_real(rtol,1e-8_dp),atol=opt_real(atol,1e-10_dp))
    case default
      rkm=rk_method_by_name(m)
      if(present(h))then
        sol=rk_integrate(rhs,y0,times,rkm,h=h,rtol=opt_real(rtol,1e-6_dp),atol=opt_real(atol,1e-8_dp))
      else
        sol=rk_integrate(rhs,y0,times,rkm,rtol=opt_real(rtol,1e-6_dp),atol=opt_real(atol,1e-8_dp))
      end if
    end select
  end function ode
  pure real(dp) function opt_real(x,default)result(v)
    real(dp),intent(in),optional::x;real(dp),intent(in)::default;v=default;if(present(x))v=x
  end function opt_real
  pure function lower(s)result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,k;out=s
    do i=1,len(s);k=iachar(out(i:i));if(k>=65.and.k<=90)out(i:i)=achar(k+32);end do
  end function lower
end module desolve
