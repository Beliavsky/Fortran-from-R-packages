module bzinb_fit
  use bzinb_kinds, only : dp
  use bzinb_distributions
  use bzinb_optimize, only : nelder_mead, golden_max
  use bzinb_linalg, only : invert_matrix, symmetrize
  use bzinb_special, only : inverse_digamma
  use bzinb_em, only : em_fit_result, bzinb_em_fit, bzinb_expectation_vec, em_expectation_result
  implicit none
  private
  public :: bp_fit_result, bzip_a_fit_result, bzip_b_fit_result
  public :: bnb_fit_result, bzinb_fit_result
  public :: fit_bp, fit_bzip_a, fit_bzip_b, fit_bnb, fit_bzinb
  public :: fit_bnb_direct, fit_bzinb_direct, fit_bnb_em, fit_bzinb_em
  public :: bzinb_standard_errors
  public :: loglik_bp, loglik_bzip_a, loglik_bzip_b, loglik_bnb, loglik_bzinb
  public :: weighted_pearson_correlation, pairwise_bzinb, inverse_digamma
  public :: pairwise_bzinb_result, pairwise_bzinb_full

  type :: bp_fit_result
    real(dp) :: param(3) = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    logical :: converged = .true.
  end type bp_fit_result

  type :: bzip_a_fit_result
    real(dp) :: param(4) = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
  end type bzip_a_fit_result

  type :: bzip_b_fit_result
    real(dp) :: param(7) = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
  end type bzip_b_fit_result

  type :: bnb_fit_result
    real(dp) :: param(5) = 0.0_dp
    real(dp) :: se(5) = 0.0_dp
    real(dp) :: covariance(5,5) = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp) :: rho_se = 0.0_dp
    real(dp) :: logit_rho = 0.0_dp
    real(dp) :: logit_rho_se = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    logical :: covariance_ok = .false.
    real(dp) :: information(5,5) = 0.0_dp
    real(dp), allocatable :: trajectory(:)
  end type bnb_fit_result

  type :: bzinb_fit_result
    real(dp) :: param(9) = 0.0_dp
    real(dp) :: se(9) = 0.0_dp
    real(dp) :: covariance(9,9) = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp) :: rho_se = 0.0_dp
    real(dp) :: logit_rho = 0.0_dp
    real(dp) :: logit_rho_se = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    logical :: covariance_ok = .false.
    real(dp) :: information(8,8) = 0.0_dp
    real(dp), allocatable :: trajectory(:)
  end type bzinb_fit_result

  type :: pairwise_bzinb_result
    integer :: npairs = 0
    integer, allocatable :: first(:), second(:)
    real(dp), allocatable :: rho(:), se_rho(:)
    real(dp), allocatable :: nonzero_first(:), nonzero_second(:), nonzero_min(:)
    real(dp), allocatable :: param(:,:), se_param(:,:), loglik(:)
    integer, allocatable :: iterations(:)
    logical, allocatable :: converged(:)
  end type pairwise_bzinb_result

  abstract interface
    function scalar_objective(z) result(v)
      import dp
      real(dp), intent(in) :: z(:)
      real(dp) :: v
    end function scalar_objective
  end interface
contains
  pure real(dp) function loglik_bp(x, y, p) result(v)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: p(3)
    integer :: i
    v = 0.0_dp
    do i = 1, size(x)
      v = v + bp_logpmf(x(i), y(i), p(1), p(2), p(3))
    end do
  end function loglik_bp

  pure real(dp) function loglik_bzip_a(x, y, p) result(v)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: p(4)
    integer :: i
    v = 0.0_dp
    do i = 1, size(x)
      v = v + bzip_a_logpmf(x(i), y(i), p(1), p(2), p(3), p(4))
    end do
  end function loglik_bzip_a

  pure real(dp) function loglik_bzip_b(x, y, p) result(v)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: p(7)
    integer :: i
    v = 0.0_dp
    do i = 1, size(x)
      v = v + bzip_b_logpmf(x(i), y(i), p(1), p(2), p(3), &
                            p(4), p(5), p(6), p(7))
    end do
  end function loglik_bzip_b

  pure real(dp) function loglik_bnb(x, y, p) result(v)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: p(5)
    integer :: i
    v = 0.0_dp
    do i = 1, size(x)
      v = v + bnb_logpmf(x(i), y(i), p(1), p(2), p(3), p(4), p(5))
    end do
  end function loglik_bnb

  pure real(dp) function loglik_bzinb(x, y, p) result(v)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: p(9)
    integer :: i
    v = 0.0_dp
    do i = 1, size(x)
      v = v + bzinb_logpmf(x(i), y(i), p(1), p(2), p(3), p(4), p(5), &
                           p(6), p(7), p(8), p(9))
    end do
  end function loglik_bzinb

  subroutine sample_stats(x, y, mx, my, vx, vy, r)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(out) :: mx, my, vx, vy, r
    integer :: n
    real(dp) :: c
    n = size(x)
    mx = sum(real(x,dp))/real(n,dp)
    my = sum(real(y,dp))/real(n,dp)
    if (n > 1) then
      vx = sum((real(x,dp)-mx)**2)/real(n-1,dp)
      vy = sum((real(y,dp)-my)**2)/real(n-1,dp)
      c = sum((real(x,dp)-mx)*(real(y,dp)-my))/real(n-1,dp)
      if (vx > 0.0_dp .and. vy > 0.0_dp) then
        r = c/sqrt(vx*vy)
      else
        r = 0.0_dp
      end if
    else
      vx = 1.0_dp
      vy = 1.0_dp
      r = 0.0_dp
    end if
  end subroutine sample_stats

  function fit_bp(x, y, tol) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    type(bp_fit_result) :: res
    real(dp) :: mx, my, m, f, t
    mx = sum(real(x,dp))/real(size(x),dp)
    my = sum(real(y,dp))/real(size(y),dp)
    m = min(mx,my)
    t = 1.0e-8_dp
    if (present(tol)) t = tol
    if (m <= 0.0_dp) then
      res%param = [0.0_dp, mx, my]
      res%loglik = loglik_bp(x,y,res%param)
      return
    end if
    call golden_max(obj, 0.0_dp, m, res%param(1), f, t, 300)
    res%param(2) = mx-res%param(1)
    res%param(3) = my-res%param(1)
    res%loglik = f
  contains
    function obj(a) result(v)
      real(dp), intent(in) :: a
      real(dp) :: v, p(3)
      p = [a, mx-a, my-a]
      v = loglik_bp(x,y,p)
    end function obj
  end function fit_bp

  function fit_bzip_a(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(4)
    type(bzip_a_fit_result) :: res
    real(dp) :: p(4), old(4), tt, q, den, eu, w
    integer :: i, s, m, it, mi
    tt = 1.0e-7_dp
    if (present(tol)) tt = tol
    mi = 5000
    if (present(maxiter)) mi = maxiter
    if (present(initial)) then
      p = initial
    else
      p = [1.0_dp,1.0_dp,1.0_dp,0.5_dp]
    end if
    do it = 1, mi
      old = p
      p = 0.0_dp
      do i = 1, size(x)
        if (x(i)+y(i) == 0) then
          q = old(4)/(old(4)+(1.0_dp-old(4))*exp(-sum(old(1:3))))
          p(1:3) = p(1:3) + old(1:3)*q
          p(4) = p(4) + q
        else
          m = min(x(i),y(i))
          den = 0.0_dp
          eu = 0.0_dp
          do s = 0, m
            w = exp(log_poisson_pmf(s,old(1)) + &
                    log_poisson_pmf(x(i)-s,old(2)) + &
                    log_poisson_pmf(y(i)-s,old(3)))
            den = den+w
            eu = eu+real(s,dp)*w
          end do
          eu = eu/max(den,tiny(1.0_dp))
          p(1) = p(1)+eu
          p(2) = p(2)+real(x(i),dp)-eu
          p(3) = p(3)+real(y(i),dp)-eu
        end if
      end do
      p = p/real(size(x),dp)
      if (maxval(abs(p-old)) <= tt) exit
    end do
    res%param = p
    res%iterations = min(it,mi)
    res%converged = it <= mi
    res%loglik = loglik_bzip_a(x,y,p)
  end function fit_bzip_a

  function fit_bzip_b(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(7)
    type(bzip_b_fit_result) :: res
    real(dp) :: p(7), old(7), post(4), sump, eu1, eu2, eu3
    real(dp) :: den1, num1, den2, num2, tt
    integer :: i, it, mi
    tt = 1.0e-7_dp
    if (present(tol)) tt = tol
    mi = 5000
    if (present(maxiter)) mi = maxiter
    if (present(initial)) then
      p = initial
    else
      p = [1.0_dp,1.0_dp,1.0_dp,0.25_dp,0.25_dp,0.25_dp,0.25_dp]
    end if
    do it = 1, mi
      old = p
      p = 0.0_dp
      den1 = 0.0_dp; num1 = 0.0_dp
      den2 = 0.0_dp; num2 = 0.0_dp
      do i = 1, size(x)
        post = 0.0_dp
        post(1) = old(4)*bp_pmf(x(i),y(i),old(1),old(2),old(3))
        if (y(i) == 0) then
          post(2) = old(5)*exp(log_poisson_pmf(x(i),old(1)+old(2)))
        end if
        if (x(i) == 0) then
          post(3) = old(6)*exp(log_poisson_pmf(y(i),old(1)+old(3)))
        end if
        if (x(i) == 0 .and. y(i) == 0) post(4) = old(7)
        sump = sum(post)
        if (sump <= 0.0_dp) cycle
        post = post/sump
        p(4:7) = p(4:7)+post
        eu1 = bp_common_mean(x(i),y(i),old(1),old(2),old(3))*post(1)
        eu2 = real(x(i),dp)*old(1)/max(old(1)+old(2),tiny(1.0_dp))*post(2)
        eu3 = real(y(i),dp)*old(1)/max(old(1)+old(3),tiny(1.0_dp))*post(3)
        p(1) = p(1)+eu1+eu2+eu3+old(1)*post(4)
        den1 = den1+post(1)+post(2)
        num1 = num1+(post(1)+post(2))*real(x(i),dp)-eu1-eu2
        den2 = den2+post(1)+post(3)
        num2 = num2+(post(1)+post(3))*real(y(i),dp)-eu1-eu3
      end do
      p(1) = p(1)/real(size(x),dp)
      p(2) = num1/max(den1,tiny(1.0_dp))
      p(3) = num2/max(den2,tiny(1.0_dp))
      p(4:7) = p(4:7)/real(size(x),dp)
      p(1:3) = max(p(1:3),1.0e-10_dp)
      p(4:7) = max(p(4:7),0.0_dp)
      p(4:7) = p(4:7)/sum(p(4:7))
      if (maxval(abs(p-old)) <= tt) exit
    end do
    res%param = p
    res%iterations = min(it,mi)
    res%converged = it <= mi
    res%loglik = loglik_bzip_b(x,y,p)
  end function fit_bzip_b

  pure real(dp) function bp_common_mean(xx, yy, m0, m1, m2) result(e)
    integer, intent(in) :: xx, yy
    real(dp), intent(in) :: m0, m1, m2
    integer :: s
    real(dp) :: w, d, n
    d = 0.0_dp; n = 0.0_dp
    do s = 0, min(xx,yy)
      w = exp(log_poisson_pmf(s,m0) + log_poisson_pmf(xx-s,m1) + &
              log_poisson_pmf(yy-s,m2))
      d = d+w
      n = n+real(s,dp)*w
    end do
    e = n/max(d,tiny(1.0_dp))
  end function bp_common_mean

  function fit_bnb_direct(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(5)
    type(bnb_fit_result) :: res
    real(dp) :: u(5), p(5), f, mx, my, vx, vy, r, tt
    integer :: mi, it
    call sample_stats(x,y,mx,my,vx,vy,r)
    if (present(initial)) then
      p = initial
    else
      p(4) = max(vx/max(mx,1.0e-4_dp),0.1_dp)
      p(5) = max(vy/max(my,1.0e-4_dp),0.1_dp)
      p(2) = max(mx/p(4),1.0e-5_dp)
      p(3) = max(my/p(5),1.0e-5_dp)
      p(1) = max(min(p(2),p(3))*abs(r),1.0e-5_dp)
      p(2) = max(p(2)-p(1),1.0e-5_dp)
      p(3) = max(p(3)-p(1),1.0e-5_dp)
    end if
    u = log(max(p,1.0e-8_dp))
    tt = 1.0e-8_dp
    if (present(tol)) tt = tol
    mi = 2500
    if (present(maxiter)) mi = maxiter
    call nelder_mead(obj,u,f,mi,tt,it,res%converged,0.12_dp)
    p = exp(u)
    res%param = p
    res%loglik = -f
    res%iterations = it
    res%rho = true_correlation(p(1),p(2),p(3),p(4),p(5))
    res%logit_rho = log(res%rho/(1.0_dp-res%rho))
    call covariance_bnb(x,y,u,res)
  contains
    function obj(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: v, pp(5)
      pp = exp(z)
      v = -loglik_bnb(x,y,pp)
      if (.not. (v < huge(v))) v = huge(v)/100.0_dp
    end function obj
  end function fit_bnb_direct

  function fit_bzinb_direct(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(9)
    type(bzinb_fit_result) :: res
    real(dp) :: u(8), p(9), f, mx, my, vx, vy, r, tt, prof(4)
    integer :: mi, it, i
    call sample_stats(x,y,mx,my,vx,vy,r)
    if (present(initial)) then
      p = initial
    else
      p(4) = max(vx/max(mx,1.0e-4_dp),0.1_dp)
      p(5) = max(vy/max(my,1.0e-4_dp),0.1_dp)
      p(2) = max(mx/max(p(4),0.1_dp),1.0e-5_dp)
      p(3) = max(my/max(p(5),0.1_dp),1.0e-5_dp)
      p(1) = max(min(p(2),p(3))*abs(r),1.0e-5_dp)
      p(2) = max(p(2)-p(1),1.0e-5_dp)
      p(3) = max(p(3)-p(1),1.0e-5_dp)
      prof = 0.0_dp
      do i = 1, size(x)
        if (x(i)>0 .and. y(i)>0) then
          prof(1)=prof(1)+1.0_dp
        else if (x(i)>0) then
          prof(2)=prof(2)+1.0_dp
        else if (y(i)>0) then
          prof(3)=prof(3)+1.0_dp
        else
          prof(4)=prof(4)+1.0_dp
        end if
      end do
      p(6:9) = max(prof/real(size(x),dp),1.0e-5_dp)
      p(6:9) = p(6:9)/sum(p(6:9))
    end if
    u(1:5) = log(max(p(1:5),1.0e-8_dp))
    u(6:8) = log(max(p(6:8),1.0e-12_dp)/max(p(9),1.0e-12_dp))
    tt = 1.0e-7_dp
    if (present(tol)) tt = tol
    mi = 3500
    if (present(maxiter)) mi = maxiter
    call nelder_mead(obj,u,f,mi,tt,it,res%converged,0.10_dp)
    call u_to_bzinb(u,p)
    res%param = p
    res%loglik = -f
    res%iterations = it
    res%rho = true_correlation(p(1),p(2),p(3),p(4),p(5))
    res%logit_rho = log(res%rho/(1.0_dp-res%rho))
    call covariance_bzinb(x,y,u,res)
  contains
    function obj(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: v, pp(9)
      call u_to_bzinb(z,pp)
      v = -loglik_bzinb(x,y,pp)
      if (.not. (v < huge(v))) v = huge(v)/100.0_dp
    end function obj
  end function fit_bzinb_direct

  function fit_bnb(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(5)
    type(bnb_fit_result) :: res
    res = fit_bnb_em(x,y,tol,maxiter,initial)
  end function fit_bnb

  function fit_bzinb(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(9)
    type(bzinb_fit_result) :: res
    res = fit_bzinb_em(x,y,tol,maxiter,initial)
  end function fit_bzinb

  function fit_bnb_em(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(5)
    type(bnb_fit_result) :: res
    real(dp) :: p0(9), mx,my,vx,vy,r,tt, cov5(5,5), dg(5),vr
    integer :: mi,i
    logical :: ok
    type(em_fit_result) :: emr
    if(maxval(x)==0 .and. maxval(y)==0) then
      res%param=1.0e-10_dp;res%loglik=0.0_dp;res%iterations=1;res%converged=.true.
      res%rho=true_correlation(res%param(1),res%param(2),res%param(3),res%param(4),res%param(5))
      res%logit_rho=log(res%rho/(1.0_dp-res%rho))
      allocate(res%trajectory(0:1));res%trajectory=0.0_dp
      return
    end if
    call sample_stats(x,y,mx,my,vx,vy,r)
    if(present(initial)) then
      p0(1:5)=initial
    else
      p0(4)=vx/max(mx,1.0e-4_dp)
      p0(5)=vy/max(my,1.0e-4_dp)
      p0(2)=mx/max(p0(4),0.1_dp)
      p0(3)=my/max(p0(5),0.1_dp)
      p0(1)=min(p0(2),p0(3))*abs(r)
      p0(2)=p0(2)-p0(1);p0(3)=p0(3)-p0(1)
      p0(1:5)=max(p0(1:5),1.0e-5_dp)
    end if
    p0(6:9)=[1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    tt=1.0e-8_dp;if(present(tol))tt=tol
    mi=50000;if(present(maxiter))mi=maxiter
    emr=bzinb_em_fit(x,y,p0,tol=tt,maxiter=mi,se=.true.,bnb=.true.)
    res%param=emr%param(1:5);res%loglik=emr%expt(1);res%iterations=emr%iterations
    res%converged=emr%converged;res%information=emr%information(1:5,1:5)
    allocate(res%trajectory(size(emr%trajectory)));res%trajectory=emr%trajectory
    call invert_matrix(res%information,cov5,ok)
    if(ok) then
      res%covariance=cov5
      do i=1,5;res%se(i)=sqrt(max(cov5(i,i),0.0_dp));end do
      res%covariance_ok=.true.
    end if
    res%rho=true_correlation(res%param(1),res%param(2),res%param(3),res%param(4),res%param(5))
    res%logit_rho=log(res%rho/(1.0_dp-res%rho))
    if(ok) then
      call rho_delta(res%param,res%covariance,dg,vr)
      res%rho_se=sqrt(max(vr,0.0_dp))
      res%logit_rho_se=res%rho_se/max(res%rho*(1.0_dp-res%rho),tiny(1.0_dp))
    end if
  end function fit_bnb_em

  function fit_bzinb_em(x, y, tol, maxiter, initial) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: initial(9)
    type(bzinb_fit_result) :: res
    real(dp) :: p0(9),mx,my,vx,vy,r,tt,prof(4),cov8(8,8),cov5tmp(5,5),dg(5),vr
    integer :: mi,i
    logical :: ok
    type(em_fit_result) :: emr
    if(maxval(x)==0 .and. maxval(y)==0) then
      res%param(1:5)=1.0e-10_dp;res%param(6:9)=[1.0_dp,0.0_dp,0.0_dp,0.0_dp]
      res%loglik=0.0_dp;res%iterations=1;res%converged=.true.
      res%rho=true_correlation(res%param(1),res%param(2),res%param(3),res%param(4),res%param(5))
      res%logit_rho=log(res%rho/(1.0_dp-res%rho))
      allocate(res%trajectory(0:1));res%trajectory=0.0_dp
      return
    end if
    call sample_stats(x,y,mx,my,vx,vy,r)
    if(present(initial)) then
      p0=initial
    else
      p0(4)=vx/max(mx,1.0e-4_dp)
      p0(5)=vy/max(my,1.0e-4_dp)
      p0(2)=mx/max(p0(4),0.1_dp);p0(3)=my/max(p0(5),0.1_dp)
      p0(1)=min(p0(2),p0(3))*abs(r);p0(2)=p0(2)-p0(1);p0(3)=p0(3)-p0(1)
      prof=0.0_dp
      do i=1,size(x)
        if(x(i)>0.and.y(i)>0)then;prof(1)=prof(1)+1
        else if(x(i)>0)then;prof(2)=prof(2)+1
        else if(y(i)>0)then;prof(3)=prof(3)+1
        else;prof(4)=prof(4)+1;end if
      end do
      p0(6:9)=prof/max(sum(prof),1.0_dp)
      p0=max(p0,1.0e-5_dp)
    end if
    tt=1.0e-8_dp;if(present(tol))tt=tol
    mi=50000;if(present(maxiter))mi=maxiter
    emr=bzinb_em_fit(x,y,p0,tol=tt,maxiter=mi,se=.true.,bnb=.false.)
    res%param=emr%param;res%loglik=emr%expt(1);res%iterations=emr%iterations
    res%converged=emr%converged;res%information=emr%information
    allocate(res%trajectory(size(emr%trajectory)));res%trajectory=emr%trajectory
    call invert_matrix(res%information,cov8,ok)
    if(ok) then
      res%covariance=0.0_dp
      res%covariance(1:8,1:8)=cov8
      do i=1,8
        res%covariance(i,9)=-sum(cov8(i,6:8));res%covariance(9,i)=res%covariance(i,9)
      end do
      res%covariance(9,9)=sum(cov8(6:8,6:8))
      call symmetrize(res%covariance)
      do i=1,9;res%se(i)=sqrt(max(res%covariance(i,i),0.0_dp));end do
      res%covariance_ok=.true.
    end if
    res%rho=true_correlation(res%param(1),res%param(2),res%param(3),res%param(4),res%param(5))
    res%logit_rho=log(res%rho/(1.0_dp-res%rho))
    if(ok) then
      cov5tmp=res%covariance(1:5,1:5)
      call rho_delta(res%param(1:5),cov5tmp,dg,vr)
      res%rho_se=sqrt(max(vr,0.0_dp))
      res%logit_rho_se=res%rho_se/max(res%rho*(1.0_dp-res%rho),tiny(1.0_dp))
    end if
  end function fit_bzinb_em

  pure subroutine u_to_bzinb(u,p)
    real(dp), intent(in) :: u(8)
    real(dp), intent(out) :: p(9)
    real(dp) :: q(4), m
    p(1:5) = exp(u(1:5))
    q(1:3) = u(6:8)
    q(4) = 0.0_dp
    m = maxval(q)
    q = exp(q-m)
    p(6:9) = q/sum(q)
  end subroutine u_to_bzinb

  subroutine numerical_hessian(fn,x,h)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(size(x),size(x))
    integer :: n, i, j
    real(dp) :: hi, hj, fpp, fpm, fmp, fmm, f0, fp, fm
    real(dp), allocatable :: z(:)
    n = size(x)
    allocate(z(n))
    h = 0.0_dp
    f0 = fn(x)
    do i = 1, n
      hi = 1.0e-4_dp*(1.0_dp+abs(x(i)))
      z=x; z(i)=x(i)+hi; fp=fn(z)
      z(i)=x(i)-hi; fm=fn(z)
      h(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
      do j = i+1, n
        hj = 1.0e-4_dp*(1.0_dp+abs(x(j)))
        z=x; z(i)=x(i)+hi; z(j)=x(j)+hj; fpp=fn(z)
        z(j)=x(j)-hj; fpm=fn(z)
        z=x; z(i)=x(i)-hi; z(j)=x(j)+hj; fmp=fn(z)
        z(j)=x(j)-hj; fmm=fn(z)
        h(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
        h(j,i)=h(i,j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine covariance_bnb(x,y,u,res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: u(5)
    type(bnb_fit_result), intent(inout) :: res
    real(dp) :: h(5,5), cu(5,5), jac(5,5), dg(5), vr
    logical :: ok
    integer :: i
    call numerical_hessian(obj,u,h)
    call invert_matrix(h,cu,ok)
    if (.not. ok) return
    jac = 0.0_dp
    do i=1,5
      jac(i,i)=res%param(i)
    end do
    res%covariance=matmul(jac,matmul(cu,transpose(jac)))
    call symmetrize(res%covariance)
    do i=1,5
      res%se(i)=sqrt(max(res%covariance(i,i),0.0_dp))
    end do
    call rho_delta(res%param,res%covariance,dg,vr)
    res%rho_se=sqrt(max(vr,0.0_dp))
    res%logit_rho_se=res%rho_se/max(res%rho*(1.0_dp-res%rho),tiny(1.0_dp))
    res%covariance_ok=.true.
  contains
    function obj(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: v, pp(5)
      pp=exp(z)
      v=-loglik_bnb(x,y,pp)
    end function obj
  end subroutine covariance_bnb

  subroutine covariance_bzinb(x,y,u,res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: u(8)
    type(bzinb_fit_result), intent(inout) :: res
    real(dp) :: h(8,8), cu(8,8), jac(9,8)
    real(dp) :: pp(9), pm(9), up(8), um(8), dh, dg(5), vr, cov5tmp(5,5)
    logical :: ok
    integer :: i, k
    call numerical_hessian(obj,u,h)
    call invert_matrix(h,cu,ok)
    if (.not. ok) return
    do k=1,8
      dh=1.0e-6_dp*(1.0_dp+abs(u(k)))
      up=u; um=u
      up(k)=up(k)+dh; um(k)=um(k)-dh
      call u_to_bzinb(up,pp)
      call u_to_bzinb(um,pm)
      jac(:,k)=(pp-pm)/(2.0_dp*dh)
    end do
    res%covariance=matmul(jac,matmul(cu,transpose(jac)))
    call symmetrize(res%covariance)
    do i=1,9
      res%se(i)=sqrt(max(res%covariance(i,i),0.0_dp))
    end do
    cov5tmp=res%covariance(1:5,1:5)
    call rho_delta(res%param(1:5),cov5tmp,dg,vr)
    res%rho_se=sqrt(max(vr,0.0_dp))
    res%logit_rho_se=res%rho_se/max(res%rho*(1.0_dp-res%rho),tiny(1.0_dp))
    res%covariance_ok=.true.
  contains
    function obj(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: v, p(9)
      call u_to_bzinb(z,p)
      v=-loglik_bzinb(x,y,p)
    end function obj
  end subroutine covariance_bzinb

  subroutine rho_delta(p,cov,dg,vr)
    real(dp), intent(in) :: p(5), cov(5,5)
    real(dp), intent(out) :: dg(5), vr
    real(dp) :: rho
    rho=true_correlation(p(1),p(2),p(3),p(4),p(5))
    dg(1)=rho*(1.0_dp/p(1)-1.0_dp/(2.0_dp*(p(1)+p(2))) &
              -1.0_dp/(2.0_dp*(p(1)+p(3))))
    dg(2)=rho*(-1.0_dp/(2.0_dp*(p(1)+p(2))))
    dg(3)=rho*(-1.0_dp/(2.0_dp*(p(1)+p(3))))
    dg(4)=rho/(2.0_dp*p(4)*(p(4)+1.0_dp))
    dg(5)=rho/(2.0_dp*p(5)*(p(5)+1.0_dp))
    vr=dot_product(dg,matmul(cov,dg))
  end subroutine rho_delta

  function bzinb_standard_errors(x,y,param) result(res)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: param(9)
    type(bzinb_fit_result) :: res
    integer, allocatable :: ux(:),uy(:),fr(:)
    integer :: i,j,nuniq
    type(em_expectation_result) :: er
    real(dp) :: cov8(8,8),cov5tmp(5,5),dg(5),vr
    logical :: ok
    res%param=param
    allocate(ux(size(x)),uy(size(x)),fr(size(x)));nuniq=0
    do i=1,size(x)
      j=0
      if(nuniq>0) then
        do j=1,nuniq
          if(ux(j)==x(i).and.uy(j)==y(i)) exit
        end do
        if(j<=nuniq .and. ux(j)==x(i).and.uy(j)==y(i)) then
          fr(j)=fr(j)+1;cycle
        end if
      end if
      nuniq=nuniq+1;ux(nuniq)=x(i);uy(nuniq)=y(i);fr(nuniq)=1
    end do
    call bzinb_expectation_vec(ux(1:nuniq),uy(1:nuniq),fr(1:nuniq),param,er,.true.,.false.)
    res%loglik=er%expt(1);res%information=er%information
    res%rho=true_correlation(param(1),param(2),param(3),param(4),param(5))
    res%logit_rho=log(res%rho/(1.0_dp-res%rho))
    call invert_matrix(res%information,cov8,ok)
    if(ok) then
      res%covariance=0.0_dp;res%covariance(1:8,1:8)=cov8
      do i=1,8
        res%covariance(i,9)=-sum(cov8(i,6:8));res%covariance(9,i)=res%covariance(i,9)
      end do
      res%covariance(9,9)=sum(cov8(6:8,6:8));call symmetrize(res%covariance)
      do i=1,9;res%se(i)=sqrt(max(res%covariance(i,i),0.0_dp));end do
      cov5tmp=res%covariance(1:5,1:5)
      call rho_delta(param(1:5),cov5tmp,dg,vr)
      res%rho_se=sqrt(max(vr,0.0_dp));res%logit_rho_se=res%rho_se/max(res%rho*(1-res%rho),tiny(1.0_dp))
      res%covariance_ok=.true.
    end if
  end function bzinb_standard_errors

  real(dp) function weighted_pearson_correlation(x,y,param) result(r)
    integer, intent(in) :: x(:), y(:)
    real(dp), intent(in) :: param(9)
    real(dp) :: w(size(x)), sw, ex, ey, exy, exx, eyy
    integer :: i
    do i=1,size(x)
      w(i)=nondropout_weight(x(i),y(i),param(1),param(2),param(3), &
                            param(4),param(5),param(6),param(7),param(8),param(9))
    end do
    sw=sum(w)
    if (sw<=0.0_dp) then
      r=0.0_dp
      return
    end if
    ex=sum(w*real(x,dp))/sw
    ey=sum(w*real(y,dp))/sw
    exy=sum(w*real(x*y,dp))/sw
    exx=sum(w*real(x*x,dp))/sw
    eyy=sum(w*real(y*y,dp))/sw
    r=(exy-ex*ey)/sqrt(max((exx-ex*ex)*(eyy-ey*ey),tiny(1.0_dp)))
  end function weighted_pearson_correlation

  subroutine pairwise_bzinb_full(data,result,nonzero_prop,full_param,nsample,maxiter)
    integer, intent(in) :: data(:,:)
    type(pairwise_bzinb_result), intent(out) :: result
    logical, intent(in), optional :: nonzero_prop, full_param
    integer, intent(in), optional :: nsample, maxiter
    integer :: p,n,total,k,i,j,mi,keep,t,idx
    integer, allocatable :: ai(:),aj(:),ord(:)
    real(dp) :: u
    logical :: do_nz,do_full
    type(bzinb_fit_result) :: f
    p=size(data,1);n=size(data,2);total=p*(p-1)/2
    do_nz=.true.;if(present(nonzero_prop))do_nz=nonzero_prop
    do_full=.false.;if(present(full_param))do_full=full_param
    keep=total;if(present(nsample))keep=min(max(nsample,1),total)
    mi=50000;if(present(maxiter))mi=maxiter
    allocate(ai(total),aj(total),ord(total));k=0
    do j=2,p
      do i=1,j-1
        k=k+1;ai(k)=i;aj(k)=j;ord(k)=k
      end do
    end do
    if(keep<total) then
      do i=1,keep
        call random_number(u)
        j=i+int(u*real(total-i+1,dp));j=min(j,total)
        t=ord(i);ord(i)=ord(j);ord(j)=t
      end do
    end if
    result%npairs=keep
    allocate(result%first(keep),result%second(keep),result%rho(keep),result%se_rho(keep))
    allocate(result%nonzero_first(keep),result%nonzero_second(keep),result%nonzero_min(keep))
    allocate(result%param(9,keep),result%se_param(9,keep),result%loglik(keep))
    allocate(result%iterations(keep),result%converged(keep))
    result%param=0.0_dp;result%se_param=0.0_dp
    do k=1,keep
      idx=ord(k);i=ai(idx);j=aj(idx)
      f=fit_bzinb(data(i,:),data(j,:),maxiter=mi)
      result%first(k)=i;result%second(k)=j;result%rho(k)=f%rho;result%se_rho(k)=f%rho_se
      result%loglik(k)=f%loglik;result%iterations(k)=f%iterations;result%converged(k)=f%converged
      if(do_nz) then
        result%nonzero_first(k)=real(count(data(i,:)/=0),dp)/real(n,dp)
        result%nonzero_second(k)=real(count(data(j,:)/=0),dp)/real(n,dp)
        result%nonzero_min(k)=min(result%nonzero_first(k),result%nonzero_second(k))
      else
        result%nonzero_first(k)=0.0_dp;result%nonzero_second(k)=0.0_dp;result%nonzero_min(k)=0.0_dp
      end if
      if(do_full) then
        result%param(:,k)=f%param;result%se_param(:,k)=f%se
      end if
    end do
  end subroutine pairwise_bzinb_full

  subroutine pairwise_bzinb(data,rho,se_rho,converged,maxiter)
    integer, intent(in) :: data(:,:)
    real(dp), intent(out) :: rho(:,:), se_rho(:,:)
    logical, intent(out) :: converged(:,:)
    integer, intent(in), optional :: maxiter
    integer :: i, j, mi
    type(bzinb_fit_result) :: f
    mi=1200
    if (present(maxiter)) mi=maxiter
    rho=0.0_dp; se_rho=0.0_dp; converged=.false.
    do i=1,size(data,1)
      rho(i,i)=1.0_dp
      converged(i,i)=.true.
      do j=i+1,size(data,1)
        f=fit_bzinb(data(i,:),data(j,:),maxiter=mi)
        rho(i,j)=f%rho; rho(j,i)=f%rho
        se_rho(i,j)=f%rho_se; se_rho(j,i)=f%rho_se
        converged(i,j)=f%converged; converged(j,i)=f%converged
      end do
    end do
  end subroutine pairwise_bzinb
end module bzinb_fit
