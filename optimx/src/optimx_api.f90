! SPDX-License-Identifier: GPL-2.0-only
module optimx_api
  use optimx_kinds, only: dp
  use optimx_types
  use optimx_eval, only: initialize_problem, grfwd, grback, grcentral, grnd, grpracma
  use optimx_solvers
  use optimx_checks
  implicit none
  private
  public :: dp, optimx_problem, optimx_control, optimx_result, optimx_multi_result
  public :: derivative_check, hessian_check, kkt_result, bounds_result, scale_result, optsp_context
  public :: optsp
  public :: initialize_problem, ctrldefault, dispdefault, checksolver, checkallsolvers
  public :: optimr, optimx, opm, multistart, polyopt, proptimr
  public :: rvmmin, rvmminu, rvmminb, rcgmin, rcgminu, rcgminb
  public :: nvm, ncg, hjn, snewton, snewtm, tn, tnbc
  public :: axsearch, bmstep, bmchk, fnchk, grchk, hesschk, kktchk, optchk
  public :: gHgen, gHgenb, grfwd, grback, grcentral, grnd, grpracma
  public :: pd_check, scalechk, opm2optimr, optimr2opm
  public :: OPTIMX_SUCCESS, OPTIMX_MAXIT, OPTIMX_SMALL_GRADIENT, OPTIMX_BAD_UPDATE
  public :: OPTIMX_INVALID_INPUT, OPTIMX_BAD_EVALUATION, OPTIMX_LINESEARCH_FAILED

  type(optsp_context), save :: optsp

contains
  function ctrldefault(npar) result(control)
    integer,intent(in),optional::npar
    type(optimx_control)::control
    if(present(npar))then
      control%maxit=500+2*max(0,npar)
      control%maxfeval=3000+10*max(0,npar)
    end if
  end function ctrldefault

  subroutine dispdefault(control,unit)
    type(optimx_control),intent(in),optional::control
    integer,intent(in),optional::unit
    type(optimx_control)::c
    integer::u
    c=optimx_control();if(present(control))c=control;u=6;if(present(unit))u=unit
    write(u,'(a,i0)')'maxit: ',c%maxit
    write(u,'(a,i0)')'maxfeval: ',c%maxfeval
    write(u,'(a,es12.4)')'reltol: ',c%reltol
    write(u,'(a,es12.4)')'gradtol: ',c%gradtol
    write(u,'(a,l1)')'maximize: ',c%maximize
  end subroutine dispdefault

  pure function lower_string(s) result(t)
    character(len=*),intent(in)::s
    character(len=len(s))::t
    integer::i,k
    do i=1,len(s)
      k=iachar(s(i:i));if(k>=iachar('A').and.k<=iachar('Z'))then;t(i:i)=achar(k+32);else;t(i:i)=s(i:i);end if
    end do
  end function lower_string

  pure logical function checksolver(method) result(ok)
    character(len=*),intent(in)::method
    character(len=:),allocatable::m
    m=trim(adjustl(lower_string(method)))
    select case(m)
    case('rvmmin','rvmminu','rvmminb','bfgs','l-bfgs-b','nvm','nlminb')
      ok=.true.
    case('rcgmin','rcgminu','rcgminb','cg','ncg')
      ok=.true.
    case('hjn','hooke-jeeves','nelder-mead','anms','nmsimplex')
      ok=.true.
    case('snewton','snewtm','tn','tnbc','newton','spg','ucminf')
      ok=.true.
    case default
      ok=.false.
    end select
  end function checksolver

  subroutine checkallsolvers(names)
    character(len=24),allocatable,intent(out)::names(:)
    allocate(names(18))
    names=[character(len=24)::'Rvmmin','Rvmminu','Rvmminb','BFGS','L-BFGS-B','nvm', &
      'Rcgmin','Rcgminu','Rcgminb','CG','ncg','hjn','Hooke-Jeeves','Nelder-Mead', &
      'snewton','snewtm','tn','tnbc']
  end subroutine checkallsolvers

  subroutine optimr(problem,x0,method,control,result)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    character(len=*),intent(in),optional::method
    type(optimx_control),intent(in),optional::control
    type(optimx_result),intent(out)::result
    character(len=32)::m
    m='Rvmmin';if(present(method))m=method
    select case(trim(adjustl(lower_string(m))))
    case('rvmmin','rvmminu','rvmminb','bfgs','l-bfgs-b','nvm','nlminb')
      call bfgs_solve(problem,x0,control,result,m)
    case('rcgmin','rcgminu','rcgminb','cg','ncg')
      call cg_solve(problem,x0,control,result,m)
    case('hjn','hooke-jeeves')
      call hj_solve(problem,x0,control,result,m)
    case('nelder-mead','anms','nmsimplex')
      call nelder_mead_solve(problem,x0,control,result,m)
    case('snewton','snewtm','tn','tnbc','newton','spg','ucminf')
      call newton_solve(problem,x0,control,result,.true.,m)
    case default
      call bfgs_solve(problem,x0,control,result,m)
      result%message='unknown method mapped to Rvmmin/BFGS'
    end select
  end subroutine optimr

  subroutine rvmmin(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call bfgs_solve(problem,x0,control,result,'Rvmmin')
  end subroutine rvmmin
  subroutine rvmminu(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call bfgs_solve(problem,x0,control,result,'Rvmminu')
  end subroutine rvmminu
  subroutine rvmminb(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call bfgs_solve(problem,x0,control,result,'Rvmminb')
  end subroutine rvmminb
  subroutine rcgmin(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call cg_solve(problem,x0,control,result,'Rcgmin')
  end subroutine rcgmin
  subroutine rcgminu(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call cg_solve(problem,x0,control,result,'Rcgminu')
  end subroutine rcgminu
  subroutine rcgminb(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call cg_solve(problem,x0,control,result,'Rcgminb')
  end subroutine rcgminb
  subroutine nvm(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call bfgs_solve(problem,x0,control,result,'nvm')
  end subroutine nvm
  subroutine ncg(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call cg_solve(problem,x0,control,result,'ncg')
  end subroutine ncg
  subroutine hjn(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call hj_solve(problem,x0,control,result,'hjn')
  end subroutine hjn
  subroutine snewton(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call newton_solve(problem,x0,control,result,.false.,'snewton')
  end subroutine snewton
  subroutine snewtm(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call newton_solve(problem,x0,control,result,.true.,'snewtm')
  end subroutine snewtm
  subroutine tn(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call newton_solve(problem,x0,control,result,.false.,'tn')
  end subroutine tn
  subroutine tnbc(problem,x0,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control;type(optimx_result),intent(out)::result
    call newton_solve(problem,x0,control,result,.true.,'tnbc')
  end subroutine tnbc

  subroutine opm(problem,x0,methods,control,result)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    character(len=*),intent(in)::methods(:)
    type(optimx_control),intent(in),optional::control
    type(optimx_multi_result),intent(out)::result
    integer::i
    allocate(result%runs(size(methods)))
    do i=1,size(methods)
      call optimr(problem,x0,methods(i),control,result%runs(i))
      if(result%best==0)then
        result%best=i
      else if(result%runs(i)%value<result%runs(result%best)%value)then
        result%best=i
      end if
    end do
  end subroutine opm

  subroutine optimx(problem,x0,methods,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    character(len=*),intent(in)::methods(:);type(optimx_control),intent(in),optional::control
    type(optimx_multi_result),intent(out)::result
    call opm(problem,x0,methods,control,result)
  end subroutine optimx

  subroutine multistart(problem,starts,method,control,result)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::starts(:,:)
    character(len=*),intent(in),optional::method
    type(optimx_control),intent(in),optional::control
    type(optimx_multi_result),intent(out)::result
    integer::i
    allocate(result%runs(size(starts,2)))
    do i=1,size(starts,2)
      call optimr(problem,starts(:,i),method,control,result%runs(i))
      if(result%best==0)then;result%best=i
      else if(result%runs(i)%value<result%runs(result%best)%value)then;result%best=i;end if
    end do
  end subroutine multistart

  subroutine polyopt(problem,x0,methods,control,result)
    type(optimx_problem),intent(in)::problem;real(dp),intent(in)::x0(:)
    character(len=*),intent(in)::methods(:);type(optimx_control),intent(in),optional::control
    type(optimx_multi_result),intent(out)::result
    real(dp)::x(size(x0));integer::i
    allocate(result%runs(size(methods)));x=x0
    do i=1,size(methods)
      call optimr(problem,x,methods(i),control,result%runs(i));x=result%runs(i)%par
    end do
    result%best=size(methods)
  end subroutine polyopt

  integer function proptimr(result,nlim) result(best)
    type(optimx_multi_result),intent(in)::result
    integer,intent(in),optional::nlim
    integer::i,limit
    best=0;if(.not.allocated(result%runs))return
    limit=size(result%runs);if(present(nlim))limit=min(limit,nlim)
    do i = 1, limit
      if (.not. result%runs(i)%converged) cycle
      if (best == 0) then
        best = i
      else if (result%runs(i)%value < result%runs(best)%value) then
        best = i
      end if
    end do
    if(best==0)best=result%best
  end function proptimr

  subroutine opm2optimr(multi,index,single)
    type(optimx_multi_result),intent(in)::multi;integer,intent(in)::index
    type(optimx_result),intent(out)::single
    if(index>=1 .and. index<=size(multi%runs))single=multi%runs(index)
  end subroutine opm2optimr

  subroutine optimr2opm(single,multi)
    type(optimx_result),intent(in)::single;type(optimx_multi_result),intent(out)::multi
    allocate(multi%runs(1));multi%runs(1)=single;multi%best=1
  end subroutine optimr2opm
end module optimx_api
