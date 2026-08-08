! SPDX-License-Identifier: CECILL-2.0
! Derived from the R package neldermead 1.0-13 and its Scilab lineage.
! See LICENSE and UPSTREAM_PROVENANCE.md.

module neldermead_frontends
  use neldermead_kinds, only : dp
  use neldermead_types
  use neldermead_core, only : neldermead_search
  implicit none
  private
  public :: fminsearch, fminbnd, fmin_gridsearch, grid_result
  public :: fminsearch_options, fminbnd_options

  type :: grid_result
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: f(:)
    logical, allocatable :: feasible(:)
  end type grid_result
contains
  function fminsearch_options(n) result(o)
    integer,intent(in)::n
    type(nm_options)::o
    o%method=nm_method_variable;o%simplex0_method=nm_simplex_pfeffer
    o%simplex0_delta_usual=0.05_dp;o%simplex0_delta_zero=0.0075_dp
    o%max_iter=200*n;o%max_fun_evals=200*n;o%tol_x_method=.false.;o%tol_fun_method=.false.
    o%tol_size_delta_f_method=.true.;o%tol_simplex_size_method=.false.
    o%tol_delta_f=1.0e-4_dp;o%tol_simplex_size_absolute=1.0e-4_dp
  end function fminsearch_options

  function fminbnd_options(n) result(o)
    integer,intent(in)::n
    type(nm_options)::o
    o%method=nm_method_box;o%simplex0_method=nm_simplex_randbounds;o%box_npoints=2*n
    o%max_iter=200*n;o%max_fun_evals=200*n;o%tol_x_method=.false.;o%tol_fun_method=.false.
    o%tol_size_delta_f_method=.false.;o%tol_simplex_size_method=.false.;o%box_termination=.true.
    o%box_tol_f = 1.0e-4_dp
    o%box_nb_match = 5
    o%box_bounds_alpha = 1.0e-6_dp
    o%box_ineq_scaling = 0.5_dp
    o%gui_alpha_min = 1.0e-5_dp
  end function fminbnd_options

  subroutine fminsearch(fn,x0,result,options,output)
    procedure(nm_objective)::fn
    real(dp),intent(in)::x0(:)
    type(nm_result),intent(out)::result
    type(nm_options),intent(in),optional::options
    procedure(nm_output_callback),optional::output
    type(nm_options)::o
    if(present(options))then;o=options;else;o=fminsearch_options(size(x0));end if
    if (present(output)) then
      call neldermead_search(fn,x0,o,result,output=output)
    else
      call neldermead_search(fn,x0,o,result)
    end if
  end subroutine fminsearch

  subroutine fminbnd(fn,x0,lower,upper,result,options,constraint,ncons,output)
    procedure(nm_objective)::fn
    real(dp),intent(in)::x0(:),lower(:),upper(:)
    type(nm_result),intent(out)::result
    type(nm_options),intent(in),optional::options
    procedure(nm_constraints),optional::constraint
    integer,intent(in),optional::ncons
    procedure(nm_output_callback),optional::output
    type(nm_options)::o;integer::nc
    if(present(options))then;o=options;else;o=fminbnd_options(size(x0));end if
    nc=0;if(present(ncons))nc=ncons
    if(present(constraint))then
      if(present(output))then;call neldermead_search(fn,x0,o,result,lower,upper,constraint,nc,output=output)
      else;call neldermead_search(fn,x0,o,result,lower,upper,constraint,nc);end if
    else
      if(present(output))then;call neldermead_search(fn,x0,o,result,lower,upper,ncons=0,output=output)
      else;call neldermead_search(fn,x0,o,result,lower,upper,ncons=0);end if
    end if
  end subroutine fminbnd

  subroutine fmin_gridsearch(fn,x0,npts,alpha,grid,lower,upper)
    procedure(nm_objective)::fn
    real(dp),intent(in)::x0(:)
    integer,intent(in)::npts
    real(dp),intent(in),optional::alpha(:),lower(:),upper(:)
    type(grid_result),intent(out)::grid
    integer::n,total,k,j,t,ii
    real(dp),allocatable::lo(:),hi(:),a(:),x(:)
    integer,allocatable::idx(:),ord(:)
    real(dp)::tmpf;real(dp),allocatable::tmpx(:)
    logical::tmpl
    if(npts<3)error stop 'npts must be at least 3'
    n=size(x0);total=npts**n;allocate(lo(n),hi(n),a(n),x(n),idx(n),ord(total),tmpx(n))
    if(present(lower).neqv.present(upper))error stop 'gridsearch requires both lower and upper bounds or neither'
    if(present(lower))then
      lo=lower;hi=upper
    else
      a=10.0_dp;if(present(alpha))then
        do j=1,n;a(j)=alpha(1+mod(j-1,size(alpha)));end do
      end if
      do j=1,n
        lo(j)=min(x0(j)/a(j),x0(j)*a(j));hi(j)=max(x0(j)/a(j),x0(j)*a(j))
        if(abs(x0(j))<=tiny(1.0_dp))then;lo(j)=-a(j);hi(j)=a(j);end if
      end do
    end if
    allocate(grid%x(n,total),grid%f(total),grid%feasible(total))
    do k=1,total
      t=k-1
      do j=1,n;idx(j)=mod(t,npts);t=t/npts;x(j)=lo(j)+real(idx(j),dp)*(hi(j)-lo(j))/real(npts-1,dp);end do
      grid%x(:,k)=x;grid%f(k)=fn(x);grid%feasible(k)=.true.;ord(k)=k
    end do
    do k=1,total-1
      ii=k;do j=k+1,total;if(grid%f(j)<grid%f(ii))ii=j;end do
      if(ii/=k)then
        tmpf=grid%f(k);grid%f(k)=grid%f(ii);grid%f(ii)=tmpf
        tmpx=grid%x(:,k);grid%x(:,k)=grid%x(:,ii);grid%x(:,ii)=tmpx
        tmpl=grid%feasible(k);grid%feasible(k)=grid%feasible(ii);grid%feasible(ii)=tmpl
      end if
    end do
  end subroutine fmin_gridsearch
end module neldermead_frontends
