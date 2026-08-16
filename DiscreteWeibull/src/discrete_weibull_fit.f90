module discrete_weibull_fit
   use rsolnp, only : solnp_problem, solnp_result, solnp_control, csolnp
   use discrete_weibull_kinds, only : dp, i64
   use discrete_weibull_numerics, only : logistic, logit, nelder_mead_2d, invert_2x2
   use discrete_weibull_dist, only : log_ddweibull, log_ddweibull3, &
      Edweibull, E2dweibull, Edweibull3, E2dweibull3
   implicit none
   private

   type, public :: dweibull_fit_result
      real(dp) :: pars(2) = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: status = 0
      integer :: iterations = 0
      character(len=120) :: message = ""
   end type dweibull_fit_result

   type, public :: fisher_result
      real(dp) :: information(2,2) = 0.0_dp
      real(dp) :: inverse(2,2) = 0.0_dp
      real(dp) :: mle(2) = 0.0_dp
      integer :: status = 0
   end type fisher_result

   type :: type1_moment_context
      real(dp) :: m1 = 0.0_dp
      real(dp) :: m2 = 0.0_dp
      real(dp) :: eps = 1.0e-4_dp
      integer :: nmax = 1000
      logical :: zero = .false.
   end type type1_moment_context

   type :: fit_context
      integer(i64), allocatable :: x(:)
      logical :: zero = .false.
      real(dp) :: eps = 1.0e-4_dp
   end type fit_context

   public :: loglikedw, lossdw, estdweibull, varFisher
   public :: loglikedw3, lossdw3, estdweibull3

contains

   real(dp) function loglikedw(par,x,zero) result(value)
      real(dp), intent(in) :: par(2)
      integer(i64), intent(in) :: x(:)
      logical, intent(in), optional :: zero
      logical :: z
      integer :: i
      real(dp) :: lp
      z = .false.
      if (present(zero)) z = zero
      if (par(1) < 0.0_dp .or. par(1) >= 1.0_dp .or. par(2) <= 0.0_dp) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      value = 0.0_dp
      do i = 1, size(x)
         lp = log_ddweibull(x(i),par(1),par(2),z)
         if (lp <= -0.5_dp*huge(1.0_dp)) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         value = value-lp
      end do
   end function loglikedw

   real(dp) function lossdw(par,x,zero,eps,nmax) result(value)
      real(dp), intent(in) :: par(2)
      integer(i64), intent(in) :: x(:)
      logical, intent(in), optional :: zero
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      real(dp) :: m1,m2,t1,t2,tol
      integer :: mx
      logical :: z
      z = .false.
      tol = 1.0e-4_dp
      mx = 1000
      if (present(zero)) z = zero
      if (present(eps)) tol = eps
      if (present(nmax)) mx = nmax
      if (par(1) < 0.0_dp .or. par(1) >= 1.0_dp .or. par(2) <= 0.0_dp) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      m1 = sum(real(x,dp))/real(size(x),dp)
      m2 = sum(real(x,dp)**2)/real(size(x),dp)
      t1 = Edweibull(par(1),par(2),tol,mx,z)
      t2 = E2dweibull(par(1),par(2),tol,mx,z)
      if (t1 >= 0.1_dp*huge(1.0_dp) .or. t2 >= 0.1_dp*huge(1.0_dp)) then
         value = huge(1.0_dp)/100.0_dp
      else
         value = (m1-t1)**2+(m2-t2)**2
      end if
   end function lossdw

   subroutine type1_moment_objective(par,value,data)
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      real(dp) :: t1,t2
      if (.not. present(data)) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      select type(ctx => data)
      type is(type1_moment_context)
         if (par(1) < 0.0_dp .or. par(1) >= 1.0_dp .or. par(2) <= 0.0_dp) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         t1 = Edweibull(par(1),par(2),ctx%eps,ctx%nmax,ctx%zero)
         t2 = E2dweibull(par(1),par(2),ctx%eps,ctx%nmax,ctx%zero)
         if (t1 >= 0.1_dp*huge(1.0_dp) .or. t2 >= 0.1_dp*huge(1.0_dp)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = (ctx%m1-t1)**2+(ctx%m2-t2)**2
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end subroutine type1_moment_objective

   real(dp) function type1_ml_objective(u,data) result(v)
      real(dp), intent(in) :: u(2)
      class(*), intent(in) :: data
      real(dp) :: p(2)
      select type(ctx => data)
      type is(fit_context)
         p = [logistic(u(1)),exp(u(2))]
         v = loglikedw(p,ctx%x,ctx%zero)
      class default
         v = huge(1.0_dp)/100.0_dp
      end select
   end function type1_ml_objective

   real(dp) function type3_ml_objective(u,data) result(v)
      real(dp), intent(in) :: u(2)
      class(*), intent(in) :: data
      real(dp) :: p(2)
      select type(ctx => data)
      type is(fit_context)
         p = [exp(u(1)),-1.0_dp+exp(u(2))]
         v = loglikedw3(p,ctx%x)
      class default
         v = huge(1.0_dp)/100.0_dp
      end select
   end function type3_ml_objective

   real(dp) function type3_mom_objective(u,data) result(v)
      real(dp), intent(in) :: u(2)
      class(*), intent(in) :: data
      real(dp) :: p(2)
      select type(ctx => data)
      type is(fit_context)
         p = [exp(u(1)),-1.0_dp+exp(u(2))]
         v = lossdw3(p,ctx%x,ctx%eps)
      class default
         v = huge(1.0_dp)/100.0_dp
      end select
   end function type3_mom_objective

   function estdweibull(x,method,zero,eps,nmax) result(res)
      integer(i64), intent(in) :: x(:)
      character(len=*), intent(in), optional :: method
      logical, intent(in), optional :: zero
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      type(dweibull_fit_result) :: res
      character(len=8) :: meth
      logical :: z
      integer(i64) :: lower
      integer :: n,y,zcount,mx,nmstat
      real(dp) :: m1,m2,q0,beta0,tol,start(2),uopt(2),fopt,qhat,bhat
      type(solnp_problem) :: problem
      type(solnp_result) :: sol
      type(solnp_control) :: control
      type(fit_context) :: ctxfit

      meth = "ML"
      if (present(method)) meth = adjustl(method)
      z = .false.
      if (present(zero)) z = zero
      tol = 1.0e-4_dp
      mx = 1000
      if (present(eps)) tol = eps
      if (present(nmax)) mx = nmax
      n = size(x)
      if (n == 0) then
         res%status = 10
         res%message = "empty sample"
         return
      end if
      lower = merge(0_i64,1_i64,z)
      if (any(x < lower)) then
         res%status = 11
         res%message = "sample outside support"
         return
      end if
      m1 = sum(real(x,dp))/real(n,dp)
      m2 = sum(real(x,dp)**2)/real(n,dp)
      beta0 = 1.0_dp
      if (z) then
         q0 = m1/(m1+1.0_dp)
      else
         if (m1 <= 1.0_dp) then
            q0 = 0.5_dp
         else
            q0 = (m1-1.0_dp)/m1
         end if
      end if
      q0 = min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,q0))

      if (trim(meth) == "P") then
         y = count(x == lower)
         if (y == 0) then
            res%status = 2
            res%message = "proportion method cannot estimate q"
            return
         end if
         qhat = 1.0_dp-real(y,dp)/real(n,dp)
         zcount = count(x == lower+1_i64)
         res%pars(1) = qhat
         if (zcount+y == n .or. zcount == 0 .or. qhat <= 0.0_dp) then
            res%status = 3
            res%pars(2) = huge(1.0_dp)
            res%message = "proportion method cannot estimate beta"
            return
         end if
         bhat = log(log(qhat-real(zcount,dp)/real(n,dp))/log(qhat))/log(2.0_dp)
         if (bhat <= 0.0_dp) then
            res%status = 4
            res%message = "invalid proportion beta estimate"
            return
         end if
         res%pars(2) = bhat
         res%objective = loglikedw(res%pars,x,z)
         return
      end if

      if (count(x <= lower+1_i64) == n) then
         res%status = 5
         res%message = "method not applicable: sample has only first two support values"
         return
      end if

      if (trim(meth) == "M") then
         problem%name = "DiscreteWeibull type-I moment fit"
         problem%n = 2
         problem%fn => type1_moment_objective
         allocate(problem%start(2),problem%lower(2),problem%upper(2))
         problem%start = [q0,beta0]
         problem%lower = [0.0_dp,1.0e-8_dp]
         problem%upper = [1.0_dp-1.0e-10_dp,100.0_dp]
         allocate(type1_moment_context :: problem%data)
         select type(ctx => problem%data)
         type is(type1_moment_context)
            ctx%m1 = m1
            ctx%m2 = m2
            ctx%eps = tol
            ctx%nmax = mx
            ctx%zero = z
         end select
         control%trace = 0
         control%max_iter = 250
         control%tol = 1.0e-8_dp
         call csolnp(problem,sol,control)
         if (.not. allocated(sol%pars)) then
            res%status = 6
            res%message = "Rsolnp returned no parameters"
            return
         end if
         res%pars = sol%pars
         res%objective = sol%objective
         res%iterations = sol%out_iterations
         res%status = sol%convergence
         res%message = trim(sol%message)
         return
      end if

      if (trim(meth) == "ML") then
         start = [logit(q0),log(beta0)]
         ctxfit%x = x
         ctxfit%zero = z
         call nelder_mead_2d(type1_ml_objective,start,uopt,fopt,res%iterations, &
              nmstat,ctxfit)
         res%pars = [logistic(uopt(1)),exp(uopt(2))]
         res%objective = fopt
         res%status = nmstat
         if (nmstat == 0) then
            res%message = "converged"
         else
            res%message = "Nelder-Mead iteration limit"
         end if
         return
      end if

      res%status = 7
      res%message = "unknown method"

   end function estdweibull

   real(dp) function loglikedw3(par,x) result(value)
      real(dp), intent(in) :: par(2)
      integer(i64), intent(in) :: x(:)
      integer :: i
      real(dp) :: lp
      if (par(1) <= 0.0_dp .or. par(2) < -1.0_dp) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      value = 0.0_dp
      do i = 1, size(x)
         lp = log_ddweibull3(x(i),par(1),par(2))
         if (lp <= -0.5_dp*huge(1.0_dp)) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         value = value-lp
      end do
   end function loglikedw3

   real(dp) function lossdw3(par,x,eps) result(value)
      real(dp), intent(in) :: par(2)
      integer(i64), intent(in) :: x(:)
      real(dp), intent(in), optional :: eps
      real(dp) :: tol,m1,m2,t1,t2
      tol = 1.0e-4_dp
      if (present(eps)) tol = eps
      if (par(1) <= 0.0_dp .or. par(2) < -1.0_dp) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      m1 = sum(real(x,dp))/real(size(x),dp)
      m2 = sum(real(x,dp)**2)/real(size(x),dp)
      t1 = Edweibull3(par(1),par(2),tol)
      t2 = E2dweibull3(par(1),par(2),tol)
      if (t1 >= 0.1_dp*huge(1.0_dp) .or. t2 >= 0.1_dp*huge(1.0_dp)) then
         value = huge(1.0_dp)/100.0_dp
      else
         value = (m1-t1)**2+(m2-t2)**2
      end if
   end function lossdw3

   function estdweibull3(x,method,eps) result(res)
      integer(i64), intent(in) :: x(:)
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: eps
      type(dweibull_fit_result) :: res
      character(len=8) :: meth
      integer :: n,y,zcount,nmstat
      real(dp) :: tol,chat,bhat,start(2),uopt(2),fopt
      type(fit_context) :: ctxfit
      meth = "P"
      if (present(method)) meth = adjustl(method)
      tol = 1.0e-4_dp
      if (present(eps)) tol = eps
      n = size(x)
      if (n == 0 .or. any(x < 0_i64)) then
         res%status = 10
         res%message = "empty sample or values outside support"
         return
      end if

      y = count(x == 0_i64)
      zcount = count(x == 1_i64)

      if (trim(meth) == "P") then
         if (y == 0 .or. y == n) then
            res%status = 2
            res%message = "proportion method cannot estimate c"
            return
         end if
         chat = -log(1.0_dp-real(y,dp)/real(n,dp))
         res%pars(1) = chat
         if (y+zcount == n .or. zcount == 0) then
            res%status = 3
            res%pars(2) = huge(1.0_dp)
            res%message = "proportion method cannot estimate beta"
            return
         end if
         bhat = log(log(1.0_dp-real(y+zcount,dp)/real(n,dp))/ &
                log(1.0_dp-real(y,dp)/real(n,dp))-1.0_dp)/log(2.0_dp)
         if (bhat < -1.0_dp) then
            res%status = 4
            res%message = "proportion beta estimate outside proper-distribution domain"
            return
         end if
         res%pars(2) = bhat
         res%objective = loglikedw3(res%pars,x)
         return
      end if

      if (count(x > 1_i64) == 0) then
         res%status = 5
         res%message = "method not applicable: no observations above 1"
         return
      end if

      start = [0.0_dp,0.0_dp]  ! c=1, beta=0 after transformation.
      ctxfit%x = x
      ctxfit%eps = tol
      if (trim(meth) == "ML") then
         call nelder_mead_2d(type3_ml_objective,start,uopt,fopt,res%iterations, &
              nmstat,ctxfit)
      else if (trim(meth) == "M") then
         call nelder_mead_2d(type3_mom_objective,start,uopt,fopt,res%iterations, &
              nmstat,ctxfit)
      else
         res%status = 6
         res%message = "unknown method"
         return
      end if

      res%pars = [exp(uopt(1)),-1.0_dp+exp(uopt(2))]
      res%objective = fopt
      res%status = nmstat
      if (nmstat == 0) then
         res%message = "converged"
      else
         res%message = "Nelder-Mead iteration limit"
      end if

   end function estdweibull3

   function varFisher(x,zero) result(fr)
      integer(i64), intent(in) :: x(:)
      logical, intent(in), optional :: zero
      type(fisher_result) :: fr
      type(dweibull_fit_result) :: fit
      logical :: z
      real(dp) :: q,b,hq,hb,f0,fp,fm,fpp,fpm,fmp,fmm
      real(dp) :: pp(2)
      integer :: invstat,n

      z = .false.
      if (present(zero)) z = zero
      fit = estdweibull(x,"ML",z)
      if (fit%status /= 0) then
         fr%status = fit%status
         return
      end if
      fr%mle = fit%pars
      q = fit%pars(1)
      b = fit%pars(2)
      n = size(x)

      hq = max(1.0e-5_dp,1.0e-4_dp*min(q,1.0_dp-q))
      hb = max(1.0e-5_dp,1.0e-4_dp*b)
      hq = min(hq,0.25_dp*min(q,1.0_dp-q))
      if (hq <= 0.0_dp) hq = 1.0e-7_dp

      pp = [q,b]
      f0 = loglikedw(pp,x,z)/real(n,dp)
      pp = [q+hq,b]
      fp = loglikedw(pp,x,z)/real(n,dp)
      pp = [q-hq,b]
      fm = loglikedw(pp,x,z)/real(n,dp)
      fr%information(1,1) = (fp-2.0_dp*f0+fm)/(hq*hq)

      pp = [q,b+hb]
      fp = loglikedw(pp,x,z)/real(n,dp)
      pp = [q,b-hb]
      fm = loglikedw(pp,x,z)/real(n,dp)
      fr%information(2,2) = (fp-2.0_dp*f0+fm)/(hb*hb)

      pp = [q+hq,b+hb]
      fpp = loglikedw(pp,x,z)/real(n,dp)
      pp = [q+hq,b-hb]
      fpm = loglikedw(pp,x,z)/real(n,dp)
      pp = [q-hq,b+hb]
      fmp = loglikedw(pp,x,z)/real(n,dp)
      pp = [q-hq,b-hb]
      fmm = loglikedw(pp,x,z)/real(n,dp)
      fr%information(1,2) = (fpp-fpm-fmp+fmm)/(4.0_dp*hq*hb)
      fr%information(2,1) = fr%information(1,2)

      call invert_2x2(fr%information,fr%inverse,invstat)
      fr%status = invstat
   end function varFisher

end module discrete_weibull_fit
