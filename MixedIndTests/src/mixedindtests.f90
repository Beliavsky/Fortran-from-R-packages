! SPDX-License-Identifier: GPL-3.0-only
module mixedindtests
  use iso_fortran_env, only : int64
  use mixedind_kinds, only : dp
  use mixedind_types
  use mixedind_rng, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_poisson
  use mixedind_special, only : normal_cdf, normal_quantile, chi_square_survival, &
    student_t_cdf, student_t_quantile, poisson_quantile, negative_binomial_quantile
  use mixedind_core, only : preparedata_core, stat_dep_core, stat_dep_serial_core, &
    sn_nonserial_core, sn_serial_core, sn_serial_vector_core, bootstrap_core, &
    moebius_nonserial_core, moebius_serial_core
  implicit none
  private

  public :: dp
  public :: prepared_data_result, pair_dependence_result, dependence_result
  public :: serial_dependence_result, sn_result, bootstrap_result
  public :: copula_test_result, moebius_result
  public :: mixedind_success, mixedind_invalid_argument, mixedind_allocation_error
  public :: mixedind_numerical_error

  public :: preparedata, stat_dep, stat_dep_ser
  public :: Sn_A, Sn_Aserial, Sn_AserialVec, Sn_serial, bootstrap
  public :: EstDep, EstDepSerial, TestIndCopula, TestIndSerCopula
  public :: TestIndSerCopulaMulti, EstDepMoebius, EstDepSerialMoebius
  public :: select_p, SimAR1Poisson, SimCopulaSeries, Finv

contains

  function preparedata(x) result(out)
    real(dp), intent(in) :: x(:)
    type(prepared_data_result) :: out
    out = preparedata_core(x)
  end function preparedata

  function stat_dep(x,y) result(out)
    real(dp), intent(in) :: x(:),y(:)
    type(pair_dependence_result) :: out
    out = stat_dep_core(x,y)
  end function stat_dep

  function stat_dep_ser(x,lag) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    type(pair_dependence_result) :: out
    out = stat_dep_serial_core(x,lag)
  end function stat_dep_ser

  function Sn_A(x,trunc_level) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: trunc_level
    type(sn_result) :: out
    out = sn_nonserial_core(x,trunc_level)
  end function Sn_A

  function Sn_Aserial(x,p,trunc_level) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p,trunc_level
    type(sn_result) :: out
    out = sn_serial_core(x,p,trunc_level)
  end function Sn_Aserial

  function Sn_AserialVec(x,p,trunc_level) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: p,trunc_level
    type(sn_result) :: out
    out = sn_serial_vector_core(x,p,trunc_level)
  end function Sn_AserialVec

  function Sn_serial(x,p) result(value)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p
    real(dp) :: value
    type(sn_result) :: out
    out = sn_serial_core(x,p,p,.false.)
    if (out%status == mixedind_success) then
      value = out%sn
    else
      value = 0.0_dp
    end if
  end function Sn_serial

  function bootstrap(multiplier,sn_multiplier,xi) result(out)
    real(dp), intent(in) :: multiplier(:,:,:),sn_multiplier(:,:),xi(:)
    type(bootstrap_result) :: out
    out = bootstrap_core(multiplier,sn_multiplier,xi)
  end function bootstrap

  function EstDep(x) result(out)
    real(dp), intent(in) :: x(:,:)
    type(dependence_result) :: out
    type(pair_dependence_result) :: pair
    integer :: n,d,i,j,df
    real(dp) :: zrho,ztau,z

    n=size(x,1); d=size(x,2)
    if(n<2 .or. d<2) then
      out%status=mixedind_invalid_argument
      return
    end if
    allocate(out%tau(d,d),out%rho(d,d),out%p_tau(d,d),out%p_rho(d,d))
    out%tau=0.0_dp; out%rho=0.0_dp
    out%p_tau=0.0_dp; out%p_rho=0.0_dp
    do i=1,d
      out%tau(i,i)=1.0_dp; out%rho(i,i)=1.0_dp
      out%p_tau(i,i)=100.0_dp; out%p_rho(i,i)=100.0_dp
    end do
    do j=1,d-1
      do i=j+1,d
        pair=stat_dep_core(x(:,i),x(:,j))
        out%tau(i,j)=pair%tau; out%tau(j,i)=pair%tau
        out%rho(i,j)=pair%rho; out%rho(j,i)=pair%rho
        zrho=sqrt(real(n,dp))*pair%rho
        z=sqrt(real(n,dp))*pair%tau
        if(pair%scale>tiny(1.0_dp)) then
          ztau=0.5_dp*z/pair%scale
        else
          ztau=0.0_dp
        end if
        out%p_tau(i,j)=200.0_dp*normal_cdf(-abs(ztau))
        out%p_tau(j,i)=out%p_tau(i,j)
        out%p_rho(i,j)=200.0_dp*normal_cdf(-abs(zrho))
        out%p_rho(j,i)=out%p_rho(i,j)
        out%lb_tau=out%lb_tau+ztau*ztau
        out%lb_rho=out%lb_rho+zrho*zrho
      end do
    end do
    df=d*(d-1)/2
    out%p_lb_tau=100.0_dp*chi_square_survival(out%lb_tau,real(df,dp))
    out%p_lb_rho=100.0_dp*chi_square_survival(out%lb_rho,real(df,dp))
  end function EstDep

  function EstDepSerial(x,lag) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    type(serial_dependence_result) :: out
    type(pair_dependence_result) :: pair
    integer :: n,j
    real(dp) :: z

    n=size(x)
    if(n<2 .or. lag<1 .or. lag>=n) then
      out%status=mixedind_invalid_argument
      return
    end if
    allocate(out%tau(lag),out%rho(lag),out%p_tau(lag),out%p_rho(lag))
    do j=1,lag
      pair=stat_dep_serial_core(x,j)
      out%tau(j)=pair%tau
      out%rho(j)=pair%rho
      out%scale=pair%scale
    end do
    do j=1,lag
      z=sqrt(real(n,dp))*out%rho(j)
      out%p_rho(j)=200.0_dp*normal_cdf(-abs(z))
      out%lb_rho=out%lb_rho+z*z
      if(out%scale>tiny(1.0_dp)) then
        z=0.5_dp*sqrt(real(n,dp))*out%tau(j)/out%scale
      else
        z=0.0_dp
      end if
      out%p_tau(j)=200.0_dp*normal_cdf(-abs(z))
      out%lb_tau=out%lb_tau+z*z
    end do
    out%p_lb_tau=100.0_dp*chi_square_survival(out%lb_tau,real(lag,dp))
    out%p_lb_rho=100.0_dp*chi_square_survival(out%lb_rho,real(lag,dp))
  end function EstDepSerial

  function TestIndCopula(x,trunc_level,b,seed) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: trunc_level,b
    integer(int64), intent(in), optional :: seed
    type(copula_test_result) :: out
    type(sn_result) :: raw
    integer :: tr,nb,n,k,it
    integer(int64) :: sd
    type(rng_state) :: rng
    real(dp), allocatable :: xi(:)
    type(bootstrap_result) :: boot

    tr=2; if(present(trunc_level)) tr=trunc_level
    nb=1000; if(present(b)) nb=b
    sd=1357911_int64; if(present(seed)) sd=seed
    raw=sn_nonserial_core(x,tr)
    if(raw%status/=mixedind_success .or. nb<1) then
      out%status=merge(raw%status,mixedind_invalid_argument,raw%status/=mixedind_success)
      return
    end if
    n=size(x,1)
    allocate(out%cvm(size(raw%stats)),out%p_cvm(size(raw%stats)), &
      out%cardinality(size(raw%cardinality)),out%subsets(size(raw%subsets,1),size(raw%subsets,2)),xi(n))
    out%cvm=raw%stats; out%cardinality=raw%cardinality; out%subsets=raw%subsets
    out%sn=raw%sn; out%p_cvm=0.0_dp; out%p_sn=0.0_dp
    call rng_seed(rng,sd)
    do it=1,nb
      do k=1,n
        xi(k)=rng_normal(rng)
      end do
      boot=bootstrap_core(raw%multiplier,raw%sn_multiplier,xi)
      where(boot%cvm>=raw%stats) out%p_cvm=out%p_cvm+1.0_dp
      if(boot%sn>=raw%sn) out%p_sn=out%p_sn+1.0_dp
    end do
    out%p_cvm=100.0_dp*out%p_cvm/real(nb,dp)
    out%p_sn=100.0_dp*out%p_sn/real(nb,dp)
    call finalize_copula_test(out,size(x,2)*(size(x,2)-1))
  end function TestIndCopula

  function TestIndSerCopula(x,p,trunc_level,b,seed) result(out)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p
    integer, intent(in), optional :: trunc_level,b
    integer(int64), intent(in), optional :: seed
    type(copula_test_result) :: out
    type(sn_result) :: raw
    integer :: tr,nb,n,k,it
    integer(int64) :: sd
    type(rng_state) :: rng
    real(dp), allocatable :: xi(:)
    type(bootstrap_result) :: boot

    tr=2; if(present(trunc_level)) tr=trunc_level
    nb=1000; if(present(b)) nb=b
    sd=2468021_int64; if(present(seed)) sd=seed
    raw=sn_serial_core(x,p,tr)
    if(raw%status/=mixedind_success .or. nb<1) then
      out%status=merge(raw%status,mixedind_invalid_argument,raw%status/=mixedind_success)
      return
    end if
    n=size(x)
    allocate(out%cvm(size(raw%stats)),out%p_cvm(size(raw%stats)), &
      out%cardinality(size(raw%cardinality)),out%subsets(size(raw%subsets,1),size(raw%subsets,2)),xi(n))
    out%cvm=raw%stats; out%cardinality=raw%cardinality; out%subsets=raw%subsets
    out%sn=raw%sn; out%p_cvm=0.0_dp; out%p_sn=0.0_dp
    call rng_seed(rng,sd)
    do it=1,nb
      do k=1,n
        xi(k)=rng_normal(rng)
      end do
      boot=bootstrap_core(raw%multiplier,raw%sn_multiplier,xi)
      where(boot%cvm>raw%stats) out%p_cvm=out%p_cvm+1.0_dp
      if(boot%sn>raw%sn) out%p_sn=out%p_sn+1.0_dp
    end do
    out%p_cvm=100.0_dp*out%p_cvm/real(nb,dp)
    out%p_sn=100.0_dp*out%p_sn/real(nb,dp)
    call finalize_copula_test(out,2*(p-1))
  end function TestIndSerCopula

  function TestIndSerCopulaMulti(x,p,trunc_level,b,seed) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: p
    integer, intent(in), optional :: trunc_level,b
    integer(int64), intent(in), optional :: seed
    type(copula_test_result) :: out
    type(sn_result) :: one,pair
    real(dp), allocatable :: msum(:,:,:),jsum(:,:),stats(:),xi(:)
    integer :: tr,nb,n,d,i,j,k,it
    integer(int64) :: sd
    type(rng_state) :: rng
    type(bootstrap_result) :: boot

    tr=2; if(present(trunc_level)) tr=trunc_level
    nb=1000; if(present(b)) nb=b
    sd=9753186_int64; if(present(seed)) sd=seed
    n=size(x,1); d=size(x,2)
    if(n<2 .or. d<1 .or. nb<1) then
      out%status=mixedind_invalid_argument
      return
    end if
    one=sn_serial_core(x(:,1),p,tr)
    if(one%status/=mixedind_success) then
      out%status=one%status
      return
    end if
    allocate(msum(size(one%multiplier,1),n,n),jsum(n,n),stats(size(one%stats)))
    msum=0.0_dp; jsum=0.0_dp; stats=0.0_dp
    out%sn=0.0_dp
    do k=1,d
      one=sn_serial_core(x(:,k),p,tr)
      msum=msum+one%multiplier; jsum=jsum+one%sn_multiplier
      stats=stats+one%stats; out%sn=out%sn+one%sn
    end do
    do i=1,d-1
      do j=i+1,d
        pair=sn_serial_vector_core(x(:,[i,j]),p,tr)
        msum=msum+pair%multiplier; jsum=jsum+pair%sn_multiplier
        stats=stats+pair%stats; out%sn=out%sn+pair%sn
      end do
    end do
    allocate(out%cvm(size(stats)),out%p_cvm(size(stats)), &
      out%cardinality(size(one%cardinality)),out%subsets(size(one%subsets,1),size(one%subsets,2)),xi(n))
    out%cvm=stats; out%p_cvm=0.0_dp
    out%cardinality=one%cardinality; out%subsets=one%subsets; out%p_sn=0.0_dp
    call rng_seed(rng,sd)
    do it=1,nb
      do k=1,n
        xi(k)=rng_normal(rng)
      end do
      boot=bootstrap_core(msum,jsum,xi)
      where(boot%cvm>stats) out%p_cvm=out%p_cvm+1.0_dp
      if(boot%sn>out%sn) out%p_sn=out%p_sn+1.0_dp
    end do
    out%p_cvm=100.0_dp*out%p_cvm/real(nb,dp)
    out%p_sn=100.0_dp*out%p_sn/real(nb,dp)
    call finalize_copula_test(out,2*(p-1))
  end function TestIndSerCopulaMulti

  subroutine finalize_copula_test(out,df_pairs)
    type(copula_test_result), intent(inout) :: out
    integer, intent(in) :: df_pairs
    real(dp), allocatable :: prob(:)
    integer :: k,npairs

    call sort_copula_by_cardinality(out)
    allocate(prob(size(out%p_cvm)))
    prob=max(out%p_cvm/100.0_dp,1.0e-20_dp)
    out%tn=-2.0_dp*sum(log(prob))
    npairs=count(out%cardinality==2)
    if(npairs>0) then
      out%tn2=-2.0_dp*sum(log(pack(prob,out%cardinality==2)))
    else
      out%tn2=0.0_dp
    end if
    out%p_tn=100.0_dp*chi_square_survival(out%tn,real(2*size(prob),dp))
    out%p_tn2=100.0_dp*chi_square_survival(out%tn2,real(df_pairs,dp))
    k=npairs
  end subroutine finalize_copula_test

  subroutine sort_copula_by_cardinality(out)
    type(copula_test_result), intent(inout) :: out
    integer :: i,j,tmpc
    real(dp) :: tmpr
    integer, allocatable :: tmps(:)
    allocate(tmps(size(out%subsets,2)))
    do i=2,size(out%cardinality)
      j=i
      do while(j>1 .and. out%cardinality(j)<out%cardinality(j-1))
        tmpc=out%cardinality(j); out%cardinality(j)=out%cardinality(j-1); out%cardinality(j-1)=tmpc
        tmpr=out%cvm(j); out%cvm(j)=out%cvm(j-1); out%cvm(j-1)=tmpr
        tmpr=out%p_cvm(j); out%p_cvm(j)=out%p_cvm(j-1); out%p_cvm(j-1)=tmpr
        tmps=out%subsets(j,:); out%subsets(j,:)=out%subsets(j-1,:); out%subsets(j-1,:)=tmps
        j=j-1
      end do
    end do
  end subroutine sort_copula_by_cardinality

  function EstDepMoebius(x,trunc_level) result(out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: trunc_level
    type(moebius_result) :: out
    integer :: tr,status,n,d,dfpairs

    tr=2; if(present(trunc_level)) tr=trunc_level
    n=size(x,1); d=size(x,2)
    call moebius_nonserial_core(x,tr,out%spearman,out%vdw,out%savage, &
      out%cardinality,out%subsets,status)
    out%status=status
    if(status/=mixedind_success) return
    call sort_moebius(out)
    out%ln_spearman=real(n,dp)*sum(out%spearman**2)
    out%ln_vdw=real(n,dp)*sum(out%vdw**2)
    out%ln_savage=real(n,dp)*sum(out%savage**2)
    out%ln2_spearman=real(n,dp)*sum(pack(out%spearman**2,out%cardinality==2))
    out%ln2_vdw=real(n,dp)*sum(pack(out%vdw**2,out%cardinality==2))
    out%ln2_savage=real(n,dp)*sum(pack(out%savage**2,out%cardinality==2))
    dfpairs=d*(d-1)/2
    call set_moebius_pvalues(out,size(out%cardinality),dfpairs)
  end function EstDepMoebius

  function EstDepSerialMoebius(y,p,trunc_level) result(out)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: p
    integer, intent(in), optional :: trunc_level
    type(moebius_result) :: out
    integer :: tr,status,n

    tr=2; if(present(trunc_level)) tr=trunc_level
    n=size(y)
    call moebius_serial_core(y,p,tr,out%spearman,out%vdw,out%savage, &
      out%cardinality,out%subsets,status)
    out%status=status
    if(status/=mixedind_success) return
    call sort_moebius(out)
    out%ln_spearman=real(n,dp)*sum(out%spearman**2)
    out%ln_vdw=real(n,dp)*sum(out%vdw**2)
    out%ln_savage=real(n,dp)*sum(out%savage**2)
    out%ln2_spearman=real(n,dp)*sum(pack(out%spearman**2,out%cardinality==2))
    out%ln2_vdw=real(n,dp)*sum(pack(out%vdw**2,out%cardinality==2))
    out%ln2_savage=real(n,dp)*sum(pack(out%savage**2,out%cardinality==2))
    call set_moebius_pvalues(out,size(out%cardinality),p-1)
  end function EstDepSerialMoebius

  subroutine sort_moebius(out)
    type(moebius_result), intent(inout) :: out
    integer :: i,j,tmpc
    real(dp) :: tmpr
    integer, allocatable :: tmps(:)
    allocate(tmps(size(out%subsets,2)))
    do i=2,size(out%cardinality)
      j=i
      do while(j>1 .and. out%cardinality(j)<out%cardinality(j-1))
        tmpc=out%cardinality(j); out%cardinality(j)=out%cardinality(j-1); out%cardinality(j-1)=tmpc
        tmpr=out%spearman(j); out%spearman(j)=out%spearman(j-1); out%spearman(j-1)=tmpr
        tmpr=out%vdw(j); out%vdw(j)=out%vdw(j-1); out%vdw(j-1)=tmpr
        tmpr=out%savage(j); out%savage(j)=out%savage(j-1); out%savage(j-1)=tmpr
        tmps=out%subsets(j,:); out%subsets(j,:)=out%subsets(j-1,:); out%subsets(j-1,:)=tmps
        j=j-1
      end do
    end do
  end subroutine sort_moebius

  subroutine set_moebius_pvalues(out,df_all,df_pairs)
    type(moebius_result), intent(inout) :: out
    integer, intent(in) :: df_all,df_pairs
    out%p_ln_spearman=100.0_dp*chi_square_survival(out%ln_spearman,real(df_all,dp))
    out%p_ln_vdw=100.0_dp*chi_square_survival(out%ln_vdw,real(df_all,dp))
    out%p_ln_savage=100.0_dp*chi_square_survival(out%ln_savage,real(df_all,dp))
    out%p_ln2_spearman=100.0_dp*chi_square_survival(out%ln2_spearman,real(df_pairs,dp))
    out%p_ln2_vdw=100.0_dp*chi_square_survival(out%ln2_vdw,real(df_pairs,dp))
    out%p_ln2_savage=100.0_dp*chi_square_survival(out%ln2_savage,real(df_pairs,dp))
  end subroutine set_moebius_pvalues

  function select_p(x,p0,d,q,lambda) result(p_selected)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: p0,d
    real(dp), intent(in), optional :: q,lambda
    integer :: p_selected
    integer :: pmin,pmax,p,k,best
    real(dp) :: q0,lam,m,cte,maxc,score,bestscore
    real(dp), allocatable :: criteria(:),scaled(:)
    type(prepared_data_result) :: pd
    type(sn_result) :: raw

    pmin=2; if(present(p0)) pmin=p0
    pmax=5; if(present(d)) pmax=d
    q0=2.4_dp; if(present(q)) q0=q
    lam=0.25_dp; if(present(lambda)) lam=lambda
    if(pmin<2 .or. pmax<pmin .or. pmax>size(x)) then
      p_selected=-1
      return
    end if
    pd=preparedata_core(x)
    m=(1.0_dp-sum(pd%pdf**2))/6.0_dp
    if(m<=0.0_dp) then
      p_selected=pmin
      return
    end if
    raw=sn_serial_core(x,pmax,pmax,.false.)
    allocate(scaled(size(raw%stats)))
    scaled=raw%stats/(m**real(raw%cardinality,dp))
    maxc=maxval(scaled)
    if(maxc<=q0*log(real(size(x),dp))) then
      cte=lam*log(real(size(x),dp))
    else
      cte=2.0_dp*lam
    end if
    allocate(criteria(pmax-pmin+1))
    do p=pmin,pmax
      raw=sn_serial_core(x,p,p,.false.)
      score=sum(raw%stats/(m**real(raw%cardinality,dp)))-real(size(raw%stats),dp)*cte
      criteria(p-pmin+1)=score
    end do
    best=1; bestscore=criteria(1)
    do k=2,size(criteria)
      if(criteria(k)>=bestscore) then
        best=k; bestscore=criteria(k)
      end if
    end do
    p_selected=pmin+best-1
  end function select_p

  function SimAR1Poisson(param,n,seed) result(x)
    real(dp), intent(in) :: param(2)
    integer, intent(in) :: n
    integer(int64), intent(in), optional :: seed
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: work(:)
    type(rng_state) :: rng
    integer(int64) :: sd
    integer :: i,n0
    real(dp) :: lambda

    if(n<1 .or. param(1)<=0.0_dp .or. param(2)<0.0_dp .or. param(2)>=1.0_dp) then
      allocate(x(0)); return
    end if
    sd=123456789_int64; if(present(seed)) sd=seed
    call rng_seed(rng,sd)
    n0=n+100
    allocate(work(n0),x(n))
    work(1)=real(rng_poisson(rng,param(1)),dp)
    do i=2,n0
      lambda=param(1)+param(2)*work(i-1)
      work(i)=real(rng_poisson(rng,lambda),dp)
    end do
    x=work(101:n0)
  end function SimAR1Poisson

  function Finv(u,k) result(x)
    real(dp), intent(in) :: u(:)
    integer, intent(in) :: k
    real(dp), allocatable :: x(:)
    real(dp), parameter :: w=0.1_dp
    real(dp) :: c0,p0,c1,p,pa,pb,uu
    integer :: i
    allocate(x(size(u)))
    c0=exp(-6.67_dp); p0=w+(1.0_dp-w)*c0; c1=0.5_dp*(1.0_dp-w)
    do i=1,size(u)
      uu=min(max(u(i),epsilon(1.0_dp)),1.0_dp-epsilon(1.0_dp))
      select case(k)
      case(1)
        x(i)=merge(1.0_dp,0.0_dp,uu>0.2_dp)
      case(2)
        x(i)=real(poisson_quantile(uu,6.0_dp),dp)
      case(3)
        x(i)=real(negative_binomial_quantile(uu,1.5_dp,0.2_dp),dp)
      case(4)
        if(uu<=p0) then
          x(i)=0.0_dp
        else
          p=max(c0,(uu-w)/(1.0_dp-w))
          x(i)=real(poisson_quantile(p,6.67_dp),dp)
        end if
      case(5)
        if(uu<c1) then
          pa=min(0.5_dp,uu/(1.0_dp-w))
          x(i)=normal_quantile(pa)
        else if(uu>w+c1) then
          pb=max(w,(uu-w)/(1.0_dp-w))
          x(i)=normal_quantile(pb)
        else
          x(i)=0.0_dp
        end if
      case(6)
        x(i)=floor(200.0_dp*normal_quantile(uu))
      case(7)
        x(i)=-1.0_dp+ceiling(1.0_dp/(1.0_dp-uu))
      case default
        x(i)=0.0_dp
      end select
    end do
  end function Finv

  function SimCopulaSeries(family,n,tau,param,seed) result(u)
    character(len=*), intent(in) :: family
    integer, intent(in) :: n
    real(dp), intent(in), optional :: tau,param
    integer(int64), intent(in), optional :: seed
    real(dp), allocatable :: u(:)
    type(rng_state) :: rng
    integer(int64) :: sd
    real(dp) :: tau0,par,theta,rho,z,cte,xq,yq,w,a
    integer :: i,pord
    character(len=:), allocatable :: fam

    if(n<1) then
      allocate(u(0)); return
    end if
    tau0=0.0_dp; if(present(tau)) tau0=tau
    par=0.0_dp; if(present(param)) par=param
    sd=11235813_int64; if(present(seed)) sd=seed
    call rng_seed(rng,sd)
    allocate(u(n))
    fam=lowercase(trim(family))
    select case(fam)
    case('ind','independence')
      do i=1,n; u(i)=rng_uniform(rng); end do
    case('tent')
      u(1)=rng_uniform(rng)
      do i=2,n; u(i)=1.999_dp*min(u(i-1),1.0_dp-u(i-1)); end do
    case('gaussian','normal')
      rho=sin(0.5_dp*acos(-1.0_dp)*tau0)
      cte=sqrt(max(0.0_dp,1.0_dp-rho*rho))
      z=rng_normal(rng); u(1)=normal_cdf(z)
      do i=2,n
        z=rho*z+cte*rng_normal(rng)
        u(i)=normal_cdf(z)
      end do
    case('t','student')
      if(par<=0.0_dp) par=5.0_dp
      rho=sin(0.5_dp*acos(-1.0_dp)*tau0)
      cte=sqrt(max(0.0_dp,1.0_dp-rho*rho))
      u(1)=rng_uniform(rng); xq=student_t_quantile(u(1),par)
      do i=2,n
        yq=student_t_quantile(rng_uniform(rng),par+1.0_dp)
        z=rho*xq+yq*cte*sqrt((par+xq*xq)/(par+1.0_dp))
        u(i)=student_t_cdf(z,par)
        xq=student_t_quantile(u(i),par)
      end do
    case('fgm')
      pord=max(2,nint(par)); theta=max(-1.0_dp,min(1.0_dp,4.5_dp*tau0))
      do i=1,min(pord-1,n); u(i)=rng_uniform(rng); end do
      do i=pord,n
        w=rng_uniform(rng)
        a=theta*product(1.0_dp-2.0_dp*u(i-pord+1:i-1))
        if(abs(a)<1.0e-14_dp) then
          u(i)=w
        else
          u(i)=2.0_dp*w/(1.0_dp+a+sqrt(max(0.0_dp,(1.0_dp+a)**2-4.0_dp*a*w)))
        end if
      end do
    case('clayton','frank','gumbel','joe','plackett')
      theta=tau_to_parameter(fam,tau0)
      u(1)=rng_uniform(rng)
      do i=2,n
        u(i)=conditional_copula_quantile(fam,u(i-1),rng_uniform(rng),theta)
      end do
    case default
      do i=1,n; u(i)=rng_uniform(rng); end do
    end select
    u=min(max(u,epsilon(1.0_dp)),1.0_dp-epsilon(1.0_dp))
  end function SimCopulaSeries

  function lowercase(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i,c
    do i=1,len(s)
      c=iachar(s(i:i))
      if(c>=iachar('A') .and. c<=iachar('Z')) then
        t(i:i)=achar(c+32)
      else
        t(i:i)=s(i:i)
      end if
    end do
  end function lowercase

  function tau_to_parameter(family,tau) result(theta)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: tau
    real(dp) :: theta,lo,hi,mid,tval
    integer :: iter
    select case(family)
    case('clayton')
      theta=2.0_dp*tau/max(1.0e-12_dp,1.0_dp-tau)
      theta=max(-0.99_dp,theta)
    case('gumbel')
      theta=1.0_dp/max(1.0e-12_dp,1.0_dp-max(0.0_dp,tau))
    case('frank')
      if(abs(tau)<1.0e-10_dp) then; theta=0.0_dp; return; end if
      if(tau>0.0_dp) then; lo=1.0e-6_dp; hi=50.0_dp; else; lo=-50.0_dp; hi=-1.0e-6_dp; end if
      do iter=1,100
        mid=0.5_dp*(lo+hi); tval=frank_tau(mid)
        if(tval<tau) then; lo=mid; else; hi=mid; end if
      end do
      theta=0.5_dp*(lo+hi)
    case('plackett')
      if(abs(tau)<1.0e-10_dp) then; theta=1.0_dp; return; end if
      lo=1.0e-4_dp; hi=1.0e4_dp
      do iter=1,120
        mid=sqrt(lo*hi); tval=plackett_tau(mid)
        if(tval<tau) then; lo=mid; else; hi=mid; end if
      end do
      theta=sqrt(lo*hi)
    case('joe')
      if(tau<=0.0_dp) then; theta=1.0_dp; return; end if
      lo=1.0_dp; hi=100.0_dp
      do iter=1,100
        mid=0.5_dp*(lo+hi); tval=joe_tau(mid)
        if(tval<tau) then; lo=mid; else; hi=mid; end if
      end do
      theta=0.5_dp*(lo+hi)
    case default
      theta=1.0_dp
    end select
  end function tau_to_parameter

  function frank_tau(theta) result(tau)
    real(dp), intent(in) :: theta
    real(dp) :: tau,integ,h,x,sumv,fx
    integer, parameter :: m=400
    integer :: i
    if(abs(theta)<1.0e-8_dp) then; tau=theta/9.0_dp; return; end if
    h=theta/real(m,dp); sumv=0.0_dp
    do i=0,m
      x=real(i,dp)*h
      if(abs(x)<1.0e-10_dp) then; fx=1.0_dp; else; fx=x/safe_expm1(x); end if
      if(i==0 .or. i==m) then; sumv=sumv+fx
      else if(mod(i,2)==0) then; sumv=sumv+2.0_dp*fx
      else; sumv=sumv+4.0_dp*fx; end if
    end do
    integ=h*sumv/3.0_dp
    tau=1.0_dp-4.0_dp/theta+4.0_dp*integ/(theta*theta)
  end function frank_tau

  function plackett_tau(theta) result(tau)
    real(dp), intent(in) :: theta
    real(dp) :: tau
    if(abs(theta-1.0_dp)<1.0e-7_dp) then
      tau=(theta-1.0_dp)/3.0_dp
    else
      tau=(theta+1.0_dp)/(theta-1.0_dp)-2.0_dp*theta*log(theta)/(theta-1.0_dp)**2
    end if
  end function plackett_tau

  function joe_tau(theta) result(tau)
    real(dp), intent(in) :: theta
    real(dp) :: tau,h,t,sumv,phi,phip,integrand
    integer, parameter :: m=1000
    integer :: i
    if(theta<=1.0_dp+1.0e-10_dp) then; tau=0.0_dp; return; end if
    h=1.0_dp/real(m,dp); sumv=0.0_dp
    do i=0,m-1
      t=(real(i,dp)+0.5_dp)*h
      phi=(1.0_dp-t)**theta
      if (phi < 1.0e-100_dp) then
        integrand=-(1.0_dp-t)/theta
      else
        phip=1.0_dp-phi
        integrand=log(phip)*phip/(theta*(1.0_dp-t)**(theta-1.0_dp))
      end if
      sumv=sumv+integrand
    end do
    tau=1.0_dp+4.0_dp*h*sumv
  end function joe_tau

  function conditional_copula_quantile(family,u,w,theta) result(v)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: u,w,theta
    real(dp) :: v,lo,hi,mid,cdf
    integer :: iter
    if((family=='frank' .and. abs(theta)<1.0e-8_dp) .or. &
       (family=='plackett' .and. abs(theta-1.0_dp)<1.0e-8_dp) .or. &
       ((family=='gumbel' .or. family=='joe') .and. abs(theta-1.0_dp)<1.0e-8_dp)) then
      v=w; return
    end if
    lo=epsilon(1.0_dp); hi=1.0_dp-epsilon(1.0_dp)
    do iter=1,80
      mid=0.5_dp*(lo+hi)
      cdf=copula_conditional_cdf(family,u,mid,theta)
      if(cdf<w) then; lo=mid; else; hi=mid; end if
    end do
    v=0.5_dp*(lo+hi)
  end function conditional_copula_quantile

  function copula_conditional_cdf(family,u,v,theta) result(p)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: u,v,theta
    real(dp) :: p,h,up,um
    h=min(1.0e-5_dp,0.25_dp*min(u,1.0_dp-u))
    h=max(h,1.0e-8_dp)
    up=min(1.0_dp-epsilon(1.0_dp),u+h)
    um=max(epsilon(1.0_dp),u-h)
    p=(copula_cdf(family,up,v,theta)-copula_cdf(family,um,v,theta))/(up-um)
    p=min(max(p,0.0_dp),1.0_dp)
  end function copula_conditional_cdf

  function copula_cdf(family,u,v,theta) result(c)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: u,v,theta
    real(dp) :: c,a,b,disc
    select case(family)
    case('clayton')
      if(abs(theta)<1.0e-10_dp) then
        c=u*v
      else
        a=u**(-theta)+v**(-theta)-1.0_dp
        if(a<=0.0_dp) then; c=0.0_dp; else; c=a**(-1.0_dp/theta); end if
      end if
    case('frank')
      if(abs(theta)<1.0e-10_dp) then
        c=u*v
      else
        c=-log(1.0_dp+(exp(-theta*u)-1.0_dp)*(exp(-theta*v)-1.0_dp)/safe_expm1(-theta))/theta
      end if
    case('gumbel')
      c=exp(-(((-log(u))**theta+(-log(v))**theta)**(1.0_dp/theta)))
    case('joe')
      a=(1.0_dp-u)**theta; b=(1.0_dp-v)**theta
      c=1.0_dp-(a+b-a*b)**(1.0_dp/theta)
    case('plackett')
      if(abs(theta-1.0_dp)<1.0e-10_dp) then
        c=u*v
      else
        a=1.0_dp+(theta-1.0_dp)*(u+v)
        disc=max(0.0_dp,a*a-4.0_dp*theta*(theta-1.0_dp)*u*v)
        c=(a-sqrt(disc))/(2.0_dp*(theta-1.0_dp))
      end if
    case default
      c=u*v
    end select
  end function copula_cdf

  elemental function safe_expm1(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v
    if (abs(x) < 1.0e-5_dp) then
      v = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x/24.0_dp)))
    else
      v = exp(x)-1.0_dp
    end if
  end function safe_expm1

end module mixedindtests
