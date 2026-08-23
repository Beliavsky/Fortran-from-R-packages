module rfast2_regression
   use rfast_special, only : dp, normal_cdf
   use rfast_linalg, only : solve_linear, inverse_matrix
   use rfast_regression, only : regression_result, lmfit, glm_logistic, glm_poisson, weighted_least_squares
   use rfast_regression_v02, only : gamma_regression, multinomial_regression, multinomial_result
   use rfast_regression_v03, only : weibull_regression, weibull_regression_result, gamma_regs, multinomial_regs, &
                                    poisson_regs, normlog_regs, tobit_mle, tobit_result
   use rfast_mle, only : ztp_mle, mle_result
   use rfast2_types, only : constrained_ls_result, scan_result
   implicit none
   private

   public :: cls_fit, robust_lm, heteroscedastic_lm
   public :: gamma_reg, gamma_reg_scan, weibull_reg, weibull_reg_scan
   public :: multinom_reg, multinom_reg_scan, binom_reg, ztp_reg
   public :: poisson_reg, logistic_reg, normlog_reg_scan, sp_logiregs
   public :: batch_logistic, tobit_reg

contains

   function add_intercept(x) result(xx)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: xx(:,:)

      allocate(xx(size(x,1),size(x,2)+1))
      xx(:,1) = 1.0_dp
      xx(:,2:) = x
   end function add_intercept

   function cls_fit(y,x,rmat,ca) result(res)
      real(dp), intent(in) :: y(:),x(:,:),rmat(:,:),ca(:)
      type(constrained_ls_result) :: res
      real(dp), allocatable :: xtx(:,:),inv(:,:),rhs(:),tmpm(:,:),mid(:,:),midinv(:,:),delta(:),corr(:)
      integer :: info,p,q

      p = size(x,2)
      q = size(rmat,1)
      allocate(xtx(p,p),inv(p,p),rhs(p),tmpm(q,p),mid(q,q),midinv(q,q),delta(q),corr(p))
      xtx = matmul(transpose(x),x)
      call inverse_matrix(xtx,inv,info)
      if (info /= 0) then
         res%status = info
         return
      end if
      rhs = matmul(transpose(x),y)
      allocate(res%ols(p),res%constrained(p))
      res%ols = matmul(inv,rhs)
      tmpm = matmul(rmat,inv)
      mid = matmul(tmpm,transpose(rmat))
      call inverse_matrix(mid,midinv,info)
      if (info /= 0) then
         res%status = info
         return
      end if
      delta = matmul(rmat,res%ols)-ca
      corr = matmul(inv,matmul(transpose(rmat),matmul(midinv,delta)))
      res%constrained = res%ols-corr
   end function cls_fit

   function robust_lm(y,x) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      type(regression_result) :: res
      type(regression_result) :: fit
      real(dp), allocatable :: xx(:,:),bread(:,:),meat(:,:),xr(:,:),cov(:,:)
      integer :: p,n,j,info

      xx = add_intercept(x)
      fit = lmfit(xx,y,.false.)
      if (fit%status /= 0) then
         res = fit
         return
      end if
      n = size(xx,1)
      p = size(xx,2)
      allocate(bread(p,p),meat(p,p),xr(n,p),cov(p,p))
      call inverse_matrix(matmul(transpose(xx),xx),bread,info)
      if (info /= 0) then
         res = fit
         res%status = info
         return
      end if
      do j = 1, p
         xr(:,j) = xx(:,j)*fit%residuals
      end do
      meat = matmul(transpose(xr),xr)
      cov = matmul(bread,matmul(meat,bread))
      res = fit
      res%covariance = cov
   end function robust_lm

   function heteroscedastic_lm(x,y) result(res)
      real(dp), intent(in) :: x(:,:),y(:)
      type(regression_result) :: res
      type(regression_result) :: first,aux
      real(dp), allocatable :: logu2(:),w(:)

      first = lmfit(x,y,.false.)
      if (first%status /= 0) then
         res = first
         return
      end if
      allocate(logu2(size(y)),w(size(y)))
      logu2 = log(max(first%residuals*first%residuals,tiny(1.0_dp)))
      aux = lmfit(x,logu2,.false.)
      if (aux%status /= 0) then
         res = aux
         return
      end if
      w = exp(-0.5_dp*aux%fitted)
      res = weighted_least_squares(x,y,w)
   end function heteroscedastic_lm

   function logistic_reg(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res
      real(dp), allocatable :: xx(:,:)

      xx = add_intercept(x)
      if (present(tol) .and. present(maxiter)) then
         res = glm_logistic(xx,y,tol,maxiter)
      else if (present(tol)) then
         res = glm_logistic(xx,y,tol=tol)
      else if (present(maxiter)) then
         res = glm_logistic(xx,y,maxiter=maxiter)
      else
         res = glm_logistic(xx,y)
      end if
   end function logistic_reg

   function poisson_reg(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res
      real(dp), allocatable :: xx(:,:)

      xx = add_intercept(x)
      if (present(tol) .and. present(maxiter)) then
         res = glm_poisson(xx,y,tol,maxiter)
      else if (present(tol)) then
         res = glm_poisson(xx,y,tol=tol)
      else if (present(maxiter)) then
         res = glm_poisson(xx,y,maxiter=maxiter)
      else
         res = glm_poisson(xx,y)
      end if
   end function poisson_reg

   function gamma_reg(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res

      if (present(tol) .and. present(maxiter)) then
         res = gamma_regression(y,x,tol,maxiter)
      else if (present(tol)) then
         res = gamma_regression(y,x,tol=tol)
      else if (present(maxiter)) then
         res = gamma_regression(y,x,maxiter=maxiter)
      else
         res = gamma_regression(y,x)
      end if
   end function gamma_reg

   function gamma_reg_scan(y,x,logged) result(out)
      real(dp), intent(in) :: y(:),x(:,:)
      logical, intent(in), optional :: logged
      real(dp) :: out(size(x,2),2)

      if (present(logged)) then
         out = gamma_regs(y,x,logged=logged)
      else
         out = gamma_regs(y,x)
      end if
   end function gamma_reg_scan

   function weibull_reg(y,x,tol,maxiter) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(weibull_regression_result) :: res

      if (present(tol) .and. present(maxiter)) then
         res = weibull_regression(y,x,tol,maxiter)
      else if (present(tol)) then
         res = weibull_regression(y,x,tol=tol)
      else if (present(maxiter)) then
         res = weibull_regression(y,x,maxiter=maxiter)
      else
         res = weibull_regression(y,x)
      end if
   end function weibull_reg

   function weibull_reg_scan(y,x,logged) result(out)
      real(dp), intent(in) :: y(:),x(:,:)
      logical, intent(in), optional :: logged
      real(dp) :: out(size(x,2),2)
      type(weibull_regression_result) :: nullfit,fit
      real(dp), allocatable :: one(:,:)
      real(dp) :: stat,pv
      integer :: j
      logical :: lg

      lg = .false.
      if (present(logged)) lg = logged
      allocate(one(size(y),1))
      one = 0.0_dp
      nullfit = weibull_regression(y,one)
      do j = 1, size(x,2)
         one(:,1) = x(:,j)
         fit = weibull_regression(y,one)
         stat = max(0.0_dp,2.0_dp*(fit%loglik-nullfit%loglik))
         pv = max(tiny(1.0_dp),1.0_dp-chi1_cdf(stat))
         out(j,1) = stat
         if (lg) then
            out(j,2) = log(pv)
         else
            out(j,2) = pv
         end if
      end do
   end function weibull_reg_scan

   real(dp) function chi1_cdf(x) result(p)
      use rfast_special, only : chisq_cdf
      real(dp), intent(in) :: x
      p = chisq_cdf(x,1.0_dp)
   end function chi1_cdf

   function multinom_reg(y,x,tol,maxiter) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(multinomial_result) :: res
      if (present(tol) .and. present(maxiter)) then
         res = multinomial_regression(y,x,tol,maxiter)
      else if (present(tol)) then
         res = multinomial_regression(y,x,tol=tol)
      else if (present(maxiter)) then
         res = multinomial_regression(y,x,maxiter=maxiter)
      else
         res = multinomial_regression(y,x)
      end if
   end function multinom_reg

   function multinom_reg_scan(y,x,logged) result(out)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: logged
      real(dp) :: out(size(x,2),2)

      if (present(logged)) then
         out = multinomial_regs(y,x,logged=logged)
      else
         out = multinomial_regs(y,x)
      end if
   end function multinom_reg_scan

   function normlog_reg_scan(y,x,logged) result(out)
      real(dp), intent(in) :: y(:),x(:,:)
      logical, intent(in), optional :: logged
      real(dp) :: out(size(x,2),2)

      if (present(logged)) then
         out = normlog_regs(y,x,logged=logged)
      else
         out = normlog_regs(y,x)
      end if
   end function normlog_reg_scan

   function binom_reg(y,ni,x,tol,maxiter) result(res)
      integer, intent(in) :: y(:),ni(:)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res
      real(dp), allocatable :: xx(:,:),b(:),bn(:),eta(:),p(:),w(:),score(:),h(:,:),step(:),wx(:,:)
      real(dp) :: eps,ll,llold
      integer :: n,d,it,mi,j,info

      xx = add_intercept(x)
      n = size(y)
      d = size(xx,2)
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mi = 100
      if (present(maxiter)) mi = maxiter
      allocate(b(d),bn(d),eta(n),p(n),w(n),score(d),h(d,d),step(d),wx(n,d))
      b = 0.0_dp
      p = real(sum(y),dp)/real(sum(ni),dp)
      b(1) = log(p(1)/(1.0_dp-p(1)))
      llold = -huge(1.0_dp)
      ll = llold
      do it = 1, mi
         eta = matmul(xx,b)
         p = 1.0_dp/(1.0_dp+exp(-max(-700.0_dp,min(700.0_dp,eta))))
         score = matmul(transpose(xx),real(y,dp)-real(ni,dp)*p)
         w = real(ni,dp)*p*(1.0_dp-p)
         do j = 1, d
            wx(:,j) = w*xx(:,j)
         end do
         h = matmul(transpose(xx),wx)
         call solve_linear(h,score,step,info)
         if (info /= 0) then
            res%status = info
            return
         end if
         bn = b+step
         eta = matmul(xx,bn)
         p = 1.0_dp/(1.0_dp+exp(-max(-700.0_dp,min(700.0_dp,eta))))
         ll = sum(real(y,dp)*log(max(p,tiny(1.0_dp)))+real(ni-y,dp)*log(max(1.0_dp-p,tiny(1.0_dp))))
         if (maxval(abs(bn-b)) <= eps .or. abs(ll-llold) <= eps) then
            b = bn
            exit
         end if
         b = bn
         llold = ll
      end do
      allocate(res%beta(d),res%fitted(n),res%residuals(n),res%covariance(d,d))
      res%beta = b
      res%fitted = p
      res%residuals = real(y,dp)-real(ni,dp)*p
      res%loglik = ll
      res%deviance = -2.0_dp*ll
      res%iterations = it
      call inverse_matrix(h,res%covariance,info)
      if (info /= 0) res%covariance = 0.0_dp
   end function binom_reg

   function ztp_reg(y,x,tol,maxiter) result(res)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(regression_result) :: res
      type(mle_result) :: ini
      real(dp), allocatable :: xx(:,:),b(:),bn(:),eta(:),mu(:),em(:),w(:),score(:),h(:,:),step(:),wx(:,:)
      real(dp) :: eps,ll,llold
      integer :: n,d,it,mi,j,info

      if (any(y <= 0)) then
         res%status = 1
         return
      end if
      xx = add_intercept(x)
      n = size(y)
      d = size(xx,2)
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mi = 100
      if (present(maxiter)) mi = maxiter
      ini = ztp_mle(y)
      allocate(b(d),bn(d),eta(n),mu(n),em(n),w(n),score(d),h(d,d),step(d),wx(n,d))
      b = 0.0_dp
      if (allocated(ini%param)) then
         b(1) = log(max(ini%param(1),tiny(1.0_dp)))
      else
         b(1) = log(max(real(sum(y),dp)/real(n,dp),tiny(1.0_dp)))
      end if
      llold = -huge(1.0_dp)
      ll = llold
      do it = 1, mi
         eta = matmul(xx,b)
         mu = exp(max(-50.0_dp,min(50.0_dp,eta)))
         em = exp(-mu)
         score = matmul(transpose(xx),real(y,dp)-mu/(1.0_dp-em))
         w = mu*(1.0_dp-em*(1.0_dp+mu))/(1.0_dp-em)**2
         do j = 1, d
            wx(:,j) = w*xx(:,j)
         end do
         h = matmul(transpose(xx),wx)
         call solve_linear(h,score,step,info)
         if (info /= 0) then
            res%status = info
            return
         end if
         bn = b+step
         eta = matmul(xx,bn)
         mu = exp(max(-50.0_dp,min(50.0_dp,eta)))
         ll = sum(real(y,dp)*eta-log_expm1_vec(mu))
         if (maxval(abs(bn-b)) <= eps .or. abs(ll-llold) <= eps) then
            b = bn
            exit
         end if
         b = bn
         llold = ll
      end do
      allocate(res%beta(d),res%fitted(n),res%residuals(n),res%covariance(d,d))
      res%beta = b
      res%fitted = mu/(1.0_dp-exp(-mu))
      res%residuals = real(y,dp)-res%fitted
      res%loglik = ll-sum(log_gamma(real(y+1,dp)))
      res%iterations = it
      call inverse_matrix(h,res%covariance,info)
      if (info /= 0) res%covariance = 0.0_dp
   end function ztp_reg

   function sp_logiregs(y,x,logged) result(out)
      real(dp), intent(in) :: y(:),x(:,:)
      logical, intent(in), optional :: logged
      real(dp) :: out(size(x,2),2),p,w,z(size(y)),zc(size(y)),den,stat,pv
      integer :: j,n
      logical :: lg

      n = size(y)
      p = sum(y)/real(n,dp)
      w = p*(1.0_dp-p)
      z = log(p/(1.0_dp-p))+(y-p)/w
      zc = z-sum(z)/real(n,dp)
      lg = .false.
      if (present(logged)) lg = logged
      do j = 1, size(x,2)
         den = sum(x(:,j)**2)-sum(x(:,j))**2/real(n,dp)
         stat = dot_product(zc,x(:,j))*sqrt(w)/sqrt(max(tiny(1.0_dp),den))
         pv = 2.0_dp*(1.0_dp-normal_cdf(abs(stat)))
         out(j,1) = stat
         if (lg) then
            out(j,2) = log(max(tiny(1.0_dp),pv))
         else
            out(j,2) = pv
         end if
      end do
   end function sp_logiregs

   function batch_logistic(y,x,k) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      integer, intent(in), optional :: k
      type(regression_result) :: res
      integer :: kk,n,d,i,lo,hi,rows,j
      type(regression_result) :: fit
      real(dp), allocatable :: xx(:,:),be(:,:),va(:,:),prec(:),b(:),p(:),fullx(:,:)

      n = size(y)
      d = size(x,2)+1
      kk = 10
      if (present(k)) kk = max(1,min(k,n))
      allocate(be(kk,d),va(kk,d))
      be = 0.0_dp
      va = huge(1.0_dp)
      do i = 1, kk
         lo = 1+(i-1)*n/kk
         hi = i*n/kk
         rows = hi-lo+1
         if (rows <= d) cycle
         xx = add_intercept(x(lo:hi,:))
         fit = glm_logistic(xx,y(lo:hi))
         if (fit%status == 0) then
            be(i,:) = fit%beta
            do j = 1, d
               va(i,j) = max(tiny(1.0_dp),fit%covariance(j,j))
            end do
         end if
      end do
      allocate(prec(d),b(d),fullx(n,d),p(n))
      do j = 1, d
         prec(j) = sum(1.0_dp/va(:,j),mask=va(:,j)<huge(1.0_dp)/2.0_dp)
         if (prec(j) > 0.0_dp) then
            b(j) = sum(be(:,j)/va(:,j),mask=va(:,j)<huge(1.0_dp)/2.0_dp)/prec(j)
         else
            b(j) = 0.0_dp
         end if
      end do
      fullx = add_intercept(x)
      p = 1.0_dp/(1.0_dp+exp(-matmul(fullx,b)))
      allocate(res%beta(d),res%fitted(n),res%residuals(n),res%covariance(d,d))
      res%beta = b
      res%fitted = p
      res%residuals = y-p
      res%covariance = 0.0_dp
      do j = 1, d
         if (prec(j) > 0.0_dp) res%covariance(j,j) = 1.0_dp/prec(j)
      end do
      res%deviance = -2.0_dp*sum(y*log(max(p,tiny(1.0_dp)))+(1.0_dp-y)*log(max(1.0_dp-p,tiny(1.0_dp))))
      res%iterations = kk
   end function batch_logistic

   function tobit_reg(y) result(res)
      real(dp), intent(in) :: y(:)
      type(tobit_result) :: res
      res = tobit_mle(y)
   end function tobit_reg

   elemental real(dp) function log_expm1_scalar(x) result(v)
      real(dp), intent(in) :: x
      if (x > 50.0_dp) then
         v = x + log(1.0_dp-exp(-x))
      else
         v = log(exp(x)-1.0_dp)
      end if
   end function log_expm1_scalar

   pure function log_expm1_vec(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v(size(x))
      integer :: i
      do i = 1, size(x)
         v(i) = log_expm1_scalar(x(i))
      end do
   end function log_expm1_vec

end module rfast2_regression
