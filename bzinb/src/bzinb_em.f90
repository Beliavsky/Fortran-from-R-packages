module bzinb_em
  use, intrinsic :: iso_fortran_env, only : real32
  use bzinb_kinds, only : dp
  use bzinb_special, only : digamma_fn, trigamma_fn
  use bzinb_distributions, only : bnb_logpmf, bzinb_logpmf
  implicit none
  private
  public :: em_expectation_result, em_fit_result
  public :: bzinb_expectation, bzinb_expectation_vec, bzinb_em_fit
  public :: upstream_inverse_digamma, idigamma

  real(dp), parameter :: eps_inv_digamma = 1.0e-9_dp
  real(dp), parameter :: eps_mstep = 1.0e-9_dp
  integer, parameter :: iter_allowance = 100

  type :: em_expectation_result
    real(dp) :: expt(12) = 0.0_dp
    real(dp) :: score(8) = 0.0_dp
    real(dp) :: information(8,8) = 0.0_dp
  end type em_expectation_result

  type :: em_fit_result
    real(dp) :: param(9) = 0.0_dp
    real(dp) :: expt(12) = 0.0_dp
    real(dp) :: information(8,8) = 0.0_dp
    real(dp), allocatable :: trajectory(:)
    integer :: iterations = 0
    logical :: converged = .false.
  end type em_fit_result
contains
  pure real(dp) function upstream_digamma(x) result(v)
    real(dp), intent(in) :: x
    if (x < 600.0_dp) then
      v = digamma_fn(x)
    else
      v = log(x - 0.5_dp)
    end if
  end function upstream_digamma

  real(dp) function upstream_inverse_digamma(x0, y) result(x)
    real(dp), intent(in) :: x0, y
    real(dp) :: h
    integer :: it
    if (x0 >= 600.0_dp) then
      x = exp(y + 0.5_dp)
      return
    end if
    x = max(x0, tiny(1.0_dp))
    h = (upstream_digamma(x) - y)/trigamma_fn(x)
    it = 0
    do while (abs(h) >= eps_inv_digamma .and. it < 10000)
      h = (upstream_digamma(x) - y)/trigamma_fn(x)
      if (h > x) h = x/2.0_dp
      x = x - h
      it = it + 1
    end do
  end function upstream_inverse_digamma

  real(dp) function idigamma(y) result(x)
    real(dp), intent(in) :: y
    x = upstream_inverse_digamma(1.0_dp,y)
  end function idigamma

  pure real(dp) function l1_value(x,y,a0,a1,a2,k,m,adj) result(v)
    integer, intent(in) :: x,y,k,m
    real(dp), intent(in) :: a0,a1,a2,adj
    v = log_gamma(a1+real(k,dp))-log_gamma(real(k+1,dp))-log_gamma(a1) + &
      log_gamma(real(x+y-m-k,dp)+a0)-log_gamma(real(x-k+1,dp)) - &
      log_gamma(a0+real(y-m,dp)) + log_gamma(real(m,dp)+a2) - &
      log_gamma(real(m+1,dp))-log_gamma(a2) + log_gamma(real(y-m,dp)+a0) - &
      log_gamma(real(y-m+1,dp))-log_gamma(a0)-adj
    v = exp(v)
  end function l1_value

  pure real(dp) function l1c_value(t1,t2,k,m,adj) result(v)
    real(dp), intent(in) :: t1,t2,adj
    integer, intent(in) :: k,m
    v = exp(real(k,dp)*log(t1)+real(m,dp)*log(t2)-adj)
  end function l1c_value

  pure real(dp) function l1ac_value(t1,t2,x,y,a0,a1,a2,k,m,adj) result(v)
    real(dp), intent(in) :: t1,t2,a0,a1,a2,adj
    integer, intent(in) :: x,y,k,m
    v = log_gamma(a1+real(k,dp))-log_gamma(real(k+1,dp))-log_gamma(a1) + &
      log_gamma(real(x+y-m-k,dp)+a0)-log_gamma(real(x-k+1,dp)) - &
      log_gamma(a0+real(y-m,dp)) + log_gamma(real(m,dp)+a2) - &
      log_gamma(real(m+1,dp))-log_gamma(a2) + log_gamma(real(y-m,dp)+a0) - &
      log_gamma(real(y-m+1,dp))-log_gamma(a0) + real(k,dp)*log(t1) + &
      real(m,dp)*log(t2)-adj
    v = exp(v)
  end function l1ac_value

  pure real(dp) function l2a_value(x,a0,a1,k,adj) result(v)
    integer, intent(in) :: x,k
    real(dp), intent(in) :: a0,a1,adj
    v = exp(log_gamma(real(x,dp)+a0-real(k,dp)) + log_gamma(real(k,dp)+a1) - &
      log_gamma(a0)-log_gamma(real(x-k+1,dp))-log_gamma(a1)- &
      log_gamma(real(k+1,dp))-adj)
  end function l2a_value

  pure real(dp) function l3a_value(y,a0,a2,m,adj) result(v)
    integer, intent(in) :: y,m
    real(dp), intent(in) :: a0,a2,adj
    v = exp(log_gamma(real(y,dp)+a0-real(m,dp)) + log_gamma(real(m,dp)+a2) - &
      log_gamma(a0)-log_gamma(real(y-m+1,dp))-log_gamma(a2)- &
      log_gamma(real(m+1,dp))-adj)
  end function l3a_value

  subroutine bzinb_expectation(x,y,freq,param,res,se,bnb)
    integer, intent(in) :: x,y,freq
    real(dp), intent(in) :: param(9)
    type(em_expectation_result), intent(inout) :: res
    logical, intent(in) :: se, bnb
    real(dp) :: a0,a1,a2,b1,b2,p1,p2,p3,p4,t1,t2,zeta
    real(dp) :: adj_a,adj_b1,adj_c,adj_sum,l1b,l2b,l3b,l4b,sum_a,sum_c,sum_ac,l_sum
    real(dp) :: r0b(4),r1b(4),r2b(4),lr0e,lr1e,lr2e,r0e,r1e,r2e
    real(dp) :: ee(4), su, ve, dvec(9), fx,fy,density
    real(dp), allocatable :: la(:,:),lc(:,:),lac(:,:),r0m(:,:),r1m(:,:),r2m(:,:)
    real(dp), allocatable :: lr0m(:,:),lr1m(:,:),lr2m(:,:),l2a(:),l3a(:)
    real(dp), allocatable :: lr0m2(:),lr0m3(:),lr1m2(:),lr2m2(:),lr1m3(:),lr2m3(:)
    integer :: i,j,ii,jj

    a0=param(1);a1=param(2);a2=param(3);b1=param(4);b2=param(5)
    p1=param(6);p2=param(7);p3=param(8);p4=param(9)
    ! Preserve upstream's explicit C++ (float) cast in t1/t2.
    t1=real(real((b1+b2+1.0_dp)/(b1+1.0_dp),real32),dp)
    t2=real(real((b1+b2+1.0_dp)/(b2+1.0_dp),real32),dp)
    zeta=merge(1.0_dp,0.0_dp,x==0 .and. y==0)
    adj_a=0.0_dp;adj_b1=0.0_dp;adj_c=0.0_dp;adj_sum=0.0_dp
    allocate(la(0:x,0:y),lc(0:x,0:y),lac(0:x,0:y),r0m(0:x,0:y),r1m(0:x,0:y),r2m(0:x,0:y))
    allocate(lr0m(0:x,0:y),lr1m(0:x,0:y),lr2m(0:x,0:y),l2a(0:x),l3a(0:y))
    allocate(lr0m2(0:x),lr0m3(0:y),lr1m2(0:x),lr2m2(0:x),lr1m3(0:y),lr2m3(0:y))

    l1b=-(real(x+y,dp)+a0)*log(1.0_dp+b1+b2)+real(x,dp)*log(b1)+real(y,dp)*log(b2) - &
      a1*log(1.0_dp+b1)-a2*log(1.0_dp+b2)
    l2b=0.0_dp
    if(y==0) l2b=exp(-(real(x,dp)+a0+a1)*log(1.0_dp+b1)+real(x,dp)*log(b1))*p2
    l3b=0.0_dp
    if(x==0) l3b=exp(-(real(y,dp)+a0+a2)*log(1.0_dp+b2)+real(y,dp)*log(b2))*p3
    l4b=merge(p4,0.0_dp,x+y==0)
    sum_a=0.0_dp;sum_c=0.0_dp
    do i=0,x
      do j=0,y
        la(i,j)=l1_value(x,y,a0,a1,a2,i,j,adj_a)
        lc(i,j)=l1c_value(t1,t2,i,j,adj_c)
        sum_a=sum_a+la(i,j);sum_c=sum_c+lc(i,j)
      end do
    end do
    if(l1b < -200.0_dp) then
      if(l2b+l3b+l4b > 0.0_dp) then
        if(log(l2b+l3b+l4b) < 0.0_dp) then
          adj_b1=-l1b-200.0_dp
          l1b=l1b+adj_b1
        end if
      end if
    end if
    l1b=exp(l1b)*p1
    do while(sum_a>0.0_dp .and. log(sum_a)>250.0_dp)
      sum_a=0.0_dp;adj_a=adj_a+200.0_dp
      do i=0,x;do j=0,y;la(i,j)=l1_value(x,y,a0,a1,a2,i,j,adj_a);sum_a=sum_a+la(i,j);end do;end do
    end do
    do i=0,x;l2a(i)=l2a_value(x,a0,a1,i,adj_a);end do
    do j=0,y;l3a(j)=l3a_value(y,a0,a2,j,adj_a);end do
    do while(sum_c>0.0_dp .and. log(sum_c)>250.0_dp)
      sum_c=0.0_dp;adj_c=adj_c+200.0_dp
      do i=0,x;do j=0,y;lc(i,j)=l1c_value(t1,t2,i,j,adj_c);sum_c=sum_c+lc(i,j);end do;end do
    end do
    sum_a=0.0_dp;sum_ac=0.0_dp
    do i=0,x;do j=0,y;sum_ac=sum_ac+la(i,j)*lc(i,j);sum_a=sum_a+la(i,j);end do;end do
    if(sum_ac>0.0_dp .and. log(sum_ac)>200.0_dp) then
      adj_a=adj_a+100.0_dp;adj_c=adj_c+100.0_dp;sum_ac=0.0_dp;sum_a=0.0_dp
      do i=0,x;do j=0,y
        la(i,j)=l1_value(x,y,a0,a1,a2,i,j,adj_a);lc(i,j)=l1c_value(t1,t2,i,j,adj_c)
        sum_ac=sum_ac+la(i,j)*lc(i,j);sum_a=sum_a+la(i,j)
      end do;end do
    else if(sum_ac>0.0_dp .and. log(sum_ac)<-100.0_dp) then
      adj_a=adj_a-200.0_dp;adj_c=adj_c-200.0_dp;sum_ac=0.0_dp;sum_a=0.0_dp
      do i=0,x;do j=0,y
        la(i,j)=l1_value(x,y,a0,a1,a2,i,j,adj_a);lc(i,j)=l1c_value(t1,t2,i,j,adj_c)
        lac(i,j)=l1ac_value(t1,t2,x,y,a0,a1,a2,i,j,adj_a+adj_c)
        sum_ac=sum_ac+lac(i,j);sum_a=sum_a+la(i,j)
      end do;end do
    end if
    do i=0,x;l2a(i)=l2a_value(x,a0,a1,i,adj_a);end do
    do j=0,y;l3a(j)=l3a_value(y,a0,a2,j,adj_a);end do
    l_sum=sum_ac*l1b+sum_a*(l2b+l3b+l4b)*exp(-adj_c)
    if(l_sum==0.0_dp) then
      adj_sum=-floor(2.0_dp*log(max(sum_ac,tiny(1.0_dp)))/3.0_dp + &
        2.0_dp*log(max(l1b,tiny(1.0_dp)))/3.0_dp)
      l_sum=sum_ac*exp(adj_sum)*l1b+sum_a*exp(adj_sum)*(l2b+l3b+l4b)*exp(-adj_c)
    end if

    r0b=[b1/(1.0_dp+b1+b2), b1/(1.0_dp+b1), b1/(1.0_dp+b2), b1]
    r1b=[b1/(1.0_dp+b1), b1/(1.0_dp+b1), b1, b1]
    r2b=[b1/(1.0_dp+b2), b1, b1/(1.0_dp+b2), b1]
    lr0e=0.0_dp;lr1e=0.0_dp;lr2e=0.0_dp;r0e=0.0_dp;r1e=0.0_dp;r2e=0.0_dp
    do i=0,x
      lr0m2(i)=l2a(i)*(digamma_fn(real(x-i,dp)+a0)+log(r0b(2)))
      lr1m2(i)=l2a(i)*(digamma_fn(real(i,dp)+a1)+log(r1b(2)))
      lr2m2(i)=l2a(i)*(digamma_fn(a2)+log(r2b(2)))
      do j=0,y
        r0m(i,j)=(real(x-i+y-j,dp)+a0)*la(i,j)
        r1m(i,j)=(real(i,dp)+a1)*la(i,j)
        r2m(i,j)=(real(j,dp)+a2)*la(i,j)
        lr0m(i,j)=(digamma_fn(real(x-i+y-j,dp)+a0)+log(r0b(1)))*la(i,j)
        lr1m(i,j)=(digamma_fn(real(i,dp)+a1)+log(r1b(1)))*la(i,j)
        lr2m(i,j)=(digamma_fn(real(j,dp)+a2)+log(r2b(1)))*la(i,j)
        r0e=r0e+r0m(i,j)*lc(i,j)*exp(adj_sum)*l1b*r0b(1) + &
          r0m(i,j)*(l2b*r0b(2)+l3b*r0b(3)+l4b*r0b(4))*exp(-adj_c+adj_sum)
        r1e=r1e+r1m(i,j)*lc(i,j)*exp(adj_sum)*l1b*r1b(1) + &
          r1m(i,j)*(l2b*r1b(2)+l3b*r1b(3)+l4b*r1b(4))*exp(-adj_c+adj_sum)
        r2e=r2e+r2m(i,j)*lc(i,j)*exp(adj_sum)*l1b*r2b(1) + &
          r2m(i,j)*(l2b*r2b(2)+l3b*r2b(3)+l4b*r2b(4))*exp(-adj_c+adj_sum)
        lr0e=lr0e+lr0m(i,j)*lc(i,j)*exp(adj_sum)*l1b
        lr1e=lr1e+lr1m(i,j)*lc(i,j)*exp(adj_sum)*l1b
        lr2e=lr2e+lr2m(i,j)*lc(i,j)*exp(adj_sum)*l1b
      end do
      lr0e=lr0e+lr0m2(i)*l2b*exp(adj_sum-adj_c)
      lr1e=lr1e+lr1m2(i)*l2b*exp(adj_sum-adj_c)
      lr2e=lr2e+lr2m2(i)*l2b*exp(adj_sum-adj_c)
    end do
    do j=0,y
      lr0m3(j)=l3a(j)*(digamma_fn(real(y-j,dp)+a0)+log(r0b(3)))
      lr1m3(j)=l3a(j)*(digamma_fn(a1)+log(r1b(3)))
      lr2m3(j)=l3a(j)*(digamma_fn(real(j,dp)+a2)+log(r2b(3)))
      lr0e=lr0e+lr0m3(j)*l3b*exp(adj_sum-adj_c)
      lr1e=lr1e+lr1m3(j)*l3b*exp(adj_sum-adj_c)
      lr2e=lr2e+lr2m3(j)*l3b*exp(adj_sum-adj_c)
    end do
    r0e=r0e/l_sum;r1e=r1e/l_sum;r2e=r2e/l_sum
    lr0e=(lr0e+(digamma_fn(a0)+log(b1))*exp(adj_sum-adj_c)*l4b)/l_sum
    lr1e=(lr1e+(digamma_fn(a1)+log(b1))*exp(adj_sum-adj_c)*l4b)/l_sum
    lr2e=(lr2e+(digamma_fn(a2)+log(b1))*exp(adj_sum-adj_c)*l4b)/l_sum
    if(.not.bnb) then
      ee(1)=sum_ac*exp(adj_sum)*l1b;ee(2)=sum_a*l2b*exp(-adj_c+adj_sum)
      ee(3)=sum_a*l3b*exp(-adj_c+adj_sum);ee(4)=sum_a*l4b*exp(-adj_c+adj_sum)
      su=sum(ee);if(su>0.0_dp)ee=ee/su
    else
      ee=[1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    end if
    ve=merge((a0+a2)*b2*(ee(2)+ee(4)),real(y,dp),y==0)
    res%expt(1)=res%expt(1)+(log(l_sum)+adj_a-adj_b1+adj_c-adj_sum)*real(freq,dp)
    res%expt(2:12)=res%expt(2:12)+[r0e,r1e,r2e,lr0e,lr1e,lr2e,ee(1),ee(2),ee(3),ee(4),ve]*real(freq,dp)

    if(se) then
      dvec=0.0_dp
      do i=0,x
        do j=0,y
          dvec(1)=dvec(1)+la(i,j)*exp(adj_sum)*lc(i,j)*digamma_fn(real(x-i+y-j,dp)+a0)
          dvec(2)=dvec(2)+la(i,j)*exp(adj_sum)*lc(i,j)*digamma_fn(a1+real(i,dp))
          dvec(3)=dvec(3)+la(i,j)*exp(adj_sum)*lc(i,j)*digamma_fn(a2+real(j,dp))
          dvec(4)=dvec(4)+la(i,j)*exp(adj_sum)*lc(i,j)*(real(x,dp)/b1 - &
            (real(x+y-i-j,dp)+a0)/(b1+b2+1.0_dp)-(a1+real(i,dp))/(b1+1.0_dp))
          dvec(5)=dvec(5)+la(i,j)*exp(adj_sum)*lc(i,j)*(real(y,dp)/b2 - &
            (real(x+y-i-j,dp)+a0)/(b1+b2+1.0_dp)-(a2+real(j,dp))/(b2+1.0_dp))
        end do
      end do
      dvec(1)=dvec(1)*l1b-sum_ac*exp(adj_sum)*l1b*(digamma_fn(a0)+log(b1+b2+1.0_dp))
      dvec(2)=dvec(2)*l1b-sum_ac*exp(adj_sum)*l1b*(digamma_fn(a1)+log(b1+1.0_dp))
      dvec(3)=dvec(3)*l1b-sum_ac*exp(adj_sum)*l1b*(digamma_fn(a2)+log(b2+1.0_dp))
      dvec(4:5)=dvec(4:5)*l1b
      dvec(1:5)=dvec(1:5)*exp(adj_a+adj_c-adj_b1-adj_sum)
      fx=0.0_dp;fy=0.0_dp
      if(y==0) then
        fx=exp(log_gamma(real(x,dp)+a0+a1)-log_gamma(real(x+1,dp))-log_gamma(a0+a1) + &
          real(x,dp)*log(b1)-(real(x,dp)+a0+a1)*log(b1+1.0_dp))
        dvec(6)=fx*(digamma_fn(real(x,dp)+a0+a1)-digamma_fn(a0+a1)-log(b1+1.0_dp))*p2
        dvec(7)=fx*(real(x,dp)/b1-(real(x,dp)+a0+a1)/(b1+1.0_dp))*p2
      end if
      if(x==0) then
        fy=exp(log_gamma(real(y,dp)+a0+a2)-log_gamma(real(y+1,dp))-log_gamma(a0+a2) + &
          real(y,dp)*log(b2)-(real(y,dp)+a0+a2)*log(b2+1.0_dp))
        dvec(8)=fy*(digamma_fn(real(y,dp)+a0+a2)-digamma_fn(a0+a2)-log(b2+1.0_dp))*p3
        dvec(9)=fy*(real(y,dp)/b2-(real(y,dp)+a0+a2)/(b2+1.0_dp))*p3
      end if
      res%score(1)=dvec(1)+dvec(6)+dvec(8)
      res%score(2)=dvec(2)+dvec(6);res%score(3)=dvec(3)+dvec(8)
      res%score(4)=dvec(4)+dvec(7);res%score(5)=dvec(5)+dvec(9)
      res%score(6)=sum_ac*l1b*exp(adj_a+adj_c-adj_b1)-zeta
      res%score(7)=fy-zeta;res%score(8)=fx-zeta
      density=l_sum*exp(adj_a-adj_b1+adj_c-adj_sum)
      if(density>0.0_dp) res%score=res%score/density
      do ii=1,8;do jj=1,8
        res%information(ii,jj)=res%information(ii,jj)+res%score(ii)*res%score(jj)*real(freq,dp)
      end do;end do
    end if
  end subroutine bzinb_expectation

  subroutine bzinb_expectation_vec(x,y,freq,param,res,se,bnb)
    integer, intent(in) :: x(:),y(:),freq(:)
    real(dp), intent(in) :: param(9)
    type(em_expectation_result), intent(out) :: res
    logical, intent(in) :: se,bnb
    integer :: i,nfreq
    res%expt=0.0_dp;res%score=0.0_dp;res%information=0.0_dp
    nfreq=sum(freq)
    do i=1,size(x)
      call bzinb_expectation(x(i),y(i),freq(i),param,res,se,bnb)
    end do
    if(nfreq>0) res%expt(2:12)=res%expt(2:12)/real(nfreq,dp)
  end subroutine bzinb_expectation_vec

  subroutine mstep_shapes(expt,param)
    real(dp), intent(in) :: expt(12)
    real(dp), intent(inout) :: param(9)
    real(dp) :: lb,h,obj,denom,idg(3),maxe
    integer :: it,k
    lb=log(param(4))
    do it=1,10
      if(maxval(expt(5:7)-lb)<600.0_dp) then
        do k=1,3
          idg(k)=upstream_inverse_digamma(param(k),expt(4+k)-lb)
        end do
        obj=log(sum(idg))+lb-log(expt(2)+expt(3)+expt(4))
        denom=1.0_dp
        do k=1,3
          denom=denom-1.0_dp/(trigamma_fn(idg(k))*sum(idg))
        end do
        h=obj/denom
      else
        maxe=maxval(expt(5:7))-0.5_dp
        obj=maxe+log(sum(exp(expt(5:7)-maxe)))-log(expt(2)+expt(3)+expt(4))
        h=obj*(1.0_dp+2.0_dp*exp(maxe+0.5_dp)-lb+0.5_dp)
      end if
      if(abs(h)<eps_mstep) exit
      lb=lb-h
    end do
    ! Upstream calls hfunc once before entering the loop, so idg should correspond
    ! to the current lb whenever the non-large branch is active.
    if(maxval(expt(5:7)-lb)<600.0_dp) then
      do k=1,3;idg(k)=upstream_inverse_digamma(param(k),expt(4+k)-lb);end do
    else
      do k=1,3;idg(k)=exp(expt(4+k)-lb+0.5_dp);end do
    end if
    param(1:3)=idg
    param(4)=exp(lb)
    param(5)=param(4)*expt(12)/(expt(2)+expt(4))
  end subroutine mstep_shapes

  function bzinb_em_fit(x,y,initial,tol,maxiter,se,bnb) result(fit)
    integer, intent(in) :: x(:),y(:)
    real(dp), intent(in) :: initial(9)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    logical, intent(in), optional :: se,bnb
    type(em_fit_result) :: fit
    integer, allocatable :: ux(:),uy(:),fr(:)
    integer :: i,j,nuniq,mi,it,iter_max,last_traj
    real(dp) :: p(9),old(9),pmax(9),emax(12),pdiff,tt
    logical :: wantse,is_bnb
    type(em_expectation_result) :: er
    real(dp), allocatable :: traj(:)

    tt=1.0e-8_dp;if(present(tol))tt=tol
    mi=50000;if(present(maxiter))mi=maxiter
    wantse=.true.;if(present(se))wantse=se
    is_bnb=.false.;if(present(bnb))is_bnb=bnb
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
    allocate(traj(0:mi));traj=0.0_dp
    p=initial;pdiff=1.0_dp;iter_max=0;emax=0.0_dp;pmax=p;er%expt=0.0_dp;last_traj=0
    do it=1,mi
      if(it<3) emax(1)=er%expt(1)
      old=p
      call bzinb_expectation_vec(ux(1:nuniq),uy(1:nuniq),fr(1:nuniq),p,er,.false.,is_bnb)
      p(6:9)=er%expt(8:11)
      call mstep_shapes(er%expt,p)
      pdiff=maxval(abs(p-old))
      if(er%expt(1)>emax(1)) then
        iter_max=it;pmax=old;emax=er%expt
      end if
      if(it>=iter_max+iter_allowance) exit
      traj(it)=er%expt(1);last_traj=it
      if(pdiff<=tt) exit
    end do
    fit%converged=.not.(it>=mi .and. pdiff>tt)
    fit%iterations=iter_max
    fit%param=pmax;fit%expt=emax
    if(wantse) then
      call bzinb_expectation_vec(ux(1:nuniq),uy(1:nuniq),fr(1:nuniq),fit%param,er,.true.,is_bnb)
      fit%expt=er%expt;fit%information=er%information
    end if
    allocate(fit%trajectory(0:max(last_traj,1)))
    fit%trajectory=traj(0:max(last_traj,1))
  end function bzinb_em_fit
end module bzinb_em
