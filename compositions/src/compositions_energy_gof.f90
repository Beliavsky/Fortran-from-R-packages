! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_energy_gof
  !! Native equivalents of the energy::eqdist.etest and energy::mvnorm.etest
  !! calls used by compositions.  The equal-distribution wrapper defaults to
  !! the source package's behavior of passing the raw composition matrix to
  !! energy (the locally computed ILR list in gsi.AcompGOFEtest is unused).
  use compositions_kinds, only: dp
  use compositions_geometry, only: ilr_rows
  use compositions_linalg, only: covariance_matrix, symmetric_eigen
  use bayesm_rng, only: rng_seed, randn, rand_uniform
  implicit none
  private

  type, public :: energy_test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp), allocatable :: replicates(:)
  end type energy_test_result

  public :: energy_ksample_statistic, acomp_energy_test
  public :: mvnorm_energy_statistic, acomp_normal_energy_test

contains

  real(dp) function energy_ksample_statistic(x,sizes) result(stat)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: sizes(:)
    integer :: g,h,i,j,ig0,ih0,ng,nh
    real(dp) :: mij,mii,mjj,dist
    if(sum(sizes)/=size(x,1).or.any(sizes<1)) error stop 'energy_ksample_statistic: size mismatch'
    stat=0.0_dp; ig0=0
    do g=1,size(sizes)-1
      ng=sizes(g); ih0=ig0+ng
      do h=g+1,size(sizes)
        nh=sizes(h); mij=0.0_dp; mii=0.0_dp; mjj=0.0_dp
        do i=1,ng
          do j=1,nh
            mij=mij+euclid(x(ig0+i,:),x(ih0+j,:))
          end do
        end do
        mij=mij/real(ng*nh,dp)
        do i=1,ng; do j=1,ng
          mii=mii+euclid(x(ig0+i,:),x(ig0+j,:))
        end do; end do
        mii=mii/real(ng*ng,dp)
        do i=1,nh; do j=1,nh
          mjj=mjj+euclid(x(ih0+i,:),x(ih0+j,:))
        end do; end do
        mjj=mjj/real(nh*nh,dp)
        dist=real(ng*nh,dp)/real(ng+nh,dp)*(2.0_dp*mij-mii-mjj)
        stat=stat+dist
        ih0=ih0+nh
      end do
      ig0=ig0+ng
    end do
  end function energy_ksample_statistic

  function acomp_energy_test(x,sizes,reps,seed,use_ilr) result(res)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: sizes(:)
    integer, intent(in), optional :: reps,seed
    logical, intent(in), optional :: use_ilr
    type(energy_test_result) :: res
    real(dp), allocatable :: z(:,:),zp(:,:)
    integer, allocatable :: idx(:)
    integer :: r,nr,n,i,j,tmp,sd
    logical :: ilr_mode
    ilr_mode=.false.; if(present(use_ilr)) ilr_mode=use_ilr
    if(ilr_mode) then; z=ilr_rows(x); else; z=x; end if
    nr=999; if(present(reps)) nr=max(0,reps)
    sd=202603; if(present(seed)) sd=seed
    res%statistic=energy_ksample_statistic(z,sizes)
    allocate(res%replicates(nr)); if(nr==0) return
    call rng_seed(sd); n=size(z,1); allocate(idx(n),zp(n,size(z,2)))
    do r=1,nr
      idx=[(i,i=1,n)]
      do i=1,n-1
        j=i+int(rand_uniform()*real(n-i+1,dp)); j=min(n,max(i,j))
        tmp=idx(i); idx(i)=idx(j); idx(j)=tmp
      end do
      zp=z(idx,:); res%replicates(r)=energy_ksample_statistic(zp,sizes)
    end do
    res%p_value=real(1+count(res%replicates>=res%statistic),dp)/real(nr+1,dp)
  end function acomp_energy_test

  real(dp) function mvnorm_energy_statistic(x) result(stat)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cov(:,:),mu(:),vals(:),vecs(:,:),isqrt(:,:),y(:,:)
    real(dp) :: mean1,mean2,mean3,const,r2
    integer :: n,d,i,j,info
    n=size(x,1); d=size(x,2)
    if(n<2.or.d<1) error stop 'mvnorm_energy_statistic: invalid dimensions'
    call covariance_matrix(x,cov,mu); call symmetric_eigen(cov,vals,vecs,info)
    if(info/=0.or.any(vals<=0.0_dp)) error stop 'mvnorm_energy_statistic: singular covariance'
    allocate(isqrt(d,d)); isqrt=matmul(vecs,matmul(diagonal(1.0_dp/sqrt(vals)),transpose(vecs)))
    y=matmul(x-spread(mu,1,n),isqrt)
    const=exp(log_gamma(0.5_dp*real(d+1,dp))-log_gamma(0.5_dp*real(d,dp)))
    mean2=2.0_dp*const; mean1=0.0_dp
    do i=1,n
      r2=sum(y(i,:)**2)
      mean1=mean1+noncentral_chi_mean(real(d,dp),r2)
    end do
    mean1=mean1/real(n,dp); mean3=0.0_dp
    do i=1,n-1; do j=i+1,n
      mean3=mean3+euclid(y(i,:),y(j,:))
    end do; end do
    mean3=2.0_dp*mean3/real(n*n,dp)
    stat=real(n,dp)*(2.0_dp*mean1-mean2-mean3)
  end function mvnorm_energy_statistic

  function acomp_normal_energy_test(x,reps,seed) result(res)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: reps,seed
    type(energy_test_result) :: res
    real(dp), allocatable :: z(:,:),sim(:,:)
    integer :: nr,sd,r,i,j
    z=ilr_rows(x); nr=999; if(present(reps)) nr=max(0,reps)
    sd=97531; if(present(seed)) sd=seed
    res%statistic=mvnorm_energy_statistic(z); allocate(res%replicates(nr))
    if(nr==0) return
    call rng_seed(sd); allocate(sim(size(z,1),size(z,2)))
    do r=1,nr
      do j=1,size(z,2); do i=1,size(z,1); sim(i,j)=randn(); end do; end do
      res%replicates(r)=mvnorm_energy_statistic(sim)
    end do
    ! energy::mvnorm.etest uses 1 - mean(t < t0), without the +1 correction.
    res%p_value=1.0_dp-real(count(res%replicates<res%statistic),dp)/real(nr,dp)
  end function acomp_normal_energy_test

  pure real(dp) function euclid(a,b) result(d)
    real(dp), intent(in) :: a(:),b(:)
    d=sqrt(sum((a-b)**2))
  end function euclid

  pure function diagonal(x) result(a)
    real(dp), intent(in) :: x(:)
    real(dp) :: a(size(x),size(x))
    integer :: i
    a=0.0_dp; do i=1,size(x); a(i,i)=x(i); end do
  end function diagonal

  real(dp) function noncentral_chi_mean(df,lambda) result(mu)
    !! Mean of a noncentral chi distribution using its Poisson-mixture
    !! representation.  This is algebraically equivalent to the 1F1 formula
    !! used by energy::mvnorm.e, but avoids a GSL dependency.
    real(dp), intent(in) :: df,lambda
    integer :: k,kmax
    real(dp) :: a,logw,term,maxlog,sumw,summ
    if(lambda<0.0_dp.or.df<=0.0_dp) error stop 'noncentral_chi_mean: invalid parameter'
    a=0.5_dp*lambda; kmax=max(40,int(a+12.0_dp*sqrt(a+1.0_dp)+40.0_dp))
    maxlog=-huge(1.0_dp)
    do k=0,kmax
      logw=-a+real(k,dp)*log(max(a,tiny(1.0_dp)))-log_gamma(real(k+1,dp))
      if(a==0.0_dp.and.k==0) logw=0.0_dp
      if(a==0.0_dp.and.k>0) logw=-huge(1.0_dp)
      maxlog=max(maxlog,logw)
    end do
    sumw=0.0_dp; summ=0.0_dp
    do k=0,kmax
      logw=-a+real(k,dp)*log(max(a,tiny(1.0_dp)))-log_gamma(real(k+1,dp))
      if(a==0.0_dp.and.k==0) logw=0.0_dp
      if(a==0.0_dp.and.k>0) cycle
      term=exp(logw-maxlog)
      sumw=sumw+term
      summ=summ+term*sqrt(2.0_dp)*exp(log_gamma(0.5_dp*(df+2.0_dp*real(k,dp)+1.0_dp)) &
        -log_gamma(0.5_dp*(df+2.0_dp*real(k,dp))))
    end do
    mu=summ/sumw
  end function noncentral_chi_mean

end module compositions_energy_gof
