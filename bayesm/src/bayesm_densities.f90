module bayesm_densities
  use bayesm_kinds, only: dp, pi, log2pi
  use bayesm_linalg, only: chol_upper, inverse_upper, trace_matrix
  implicit none
  private
  public :: lnd_mvn, lnd_mvst, lnd_iwishart, lnd_ichisq
  interface lnd_ichisq
    module procedure lnd_ichisq_scalar
    module procedure lnd_ichisq_vector
    module procedure lnd_ichisq_matrix
  end interface
contains
  pure real(dp) function lnd_mvn(x,mu,rooti) result(v)
    real(dp), intent(in) :: x(:),mu(:),rooti(:,:)
    real(dp) :: z(size(x))
    integer :: i
    z=matmul(transpose(rooti),x-mu)
    v=-0.5_dp*real(size(x),dp)*log2pi-0.5_dp*dot_product(z,z)
    do i=1,size(x)
      v=v+log(abs(rooti(i,i)))
    end do
  end function lnd_mvn

  pure real(dp) function lnd_mvst(x,nu,mu,rooti,normc) result(v)
    real(dp), intent(in) :: x(:),nu,mu(:),rooti(:,:)
    logical, intent(in), optional :: normc
    real(dp) :: z(size(x)),constant
    integer :: i,d
    logical :: nc
    d=size(x); nc=.false.; if (present(normc)) nc=normc
    constant=0.0_dp
    if (nc) then
      constant=0.5_dp*nu*log(nu)+log_gamma(0.5_dp*(nu+real(d,dp))) &
        -0.5_dp*real(d,dp)*log(pi)-log_gamma(0.5_dp*nu)
    end if
    z=matmul(transpose(rooti),x-mu)
    v=constant-0.5_dp*(real(d,dp)+nu)*log(nu+dot_product(z,z))
    do i=1,d
      v=v+log(abs(rooti(i,i)))
    end do
  end function lnd_mvst

  real(dp) function lnd_iwishart(nu,vmat,iw) result(val)
    real(dp), intent(in) :: nu,vmat(:,:),iw(:,:)
    integer :: k,i,info
    real(dp) :: uiw(size(iw,1),size(iw,2)),uiwi(size(iw,1),size(iw,2))
    real(dp) :: viwi(size(iw,1),size(iw,2)),iw_inv(size(iw,1),size(iw,2))
    real(dp) :: cholv(size(vmat,1),size(vmat,2)),cnst,lndetvd2,lndetiwd2,arg
    k=size(vmat,1)
    call chol_upper(iw,uiw,info)
    if (info/=0) then; val=-huge(1.0_dp); return; end if
    call inverse_upper(uiw,uiwi,info)
    if (info/=0) then; val=-huge(1.0_dp); return; end if
    iw_inv=matmul(uiwi,transpose(uiwi))
    call chol_upper(vmat,cholv,info)
    if (info/=0) then; val=-huge(1.0_dp); return; end if
    lndetvd2=0.0_dp; lndetiwd2=0.0_dp
    do i=1,k
      lndetvd2=lndetvd2+log(cholv(i,i)); lndetiwd2=lndetiwd2+log(uiw(i,i))
    end do
    cnst=0.5_dp*nu*real(k,dp)*log(2.0_dp)+0.25_dp*real(k*(k-1),dp)*log(pi)
    do i=1,k
      arg=0.5_dp*(nu+1.0_dp-real(i,dp)); cnst=cnst+log_gamma(arg)
    end do
    viwi=matmul(vmat,iw_inv)
    val=-cnst+nu*lndetvd2-(nu+real(k+1,dp))*lndetiwd2-0.5_dp*trace_matrix(viwi)
  end function lnd_iwishart

  pure elemental real(dp) function lnd_ichisq_scalar(nu,ssq,x) result(v)
    real(dp), intent(in) :: nu,ssq,x
    if (x<=0.0_dp) then
      v=-huge(1.0_dp)
    else
      v=-log_gamma(0.5_dp*nu)+0.5_dp*nu*log(0.5_dp*nu*ssq) &
        -(0.5_dp*nu+1.0_dp)*log(x)-0.5_dp*nu*ssq/x
    end if
  end function lnd_ichisq_scalar

  pure function lnd_ichisq_vector(nu,ssq,x) result(v)
    real(dp), intent(in) :: nu,ssq,x(:)
    real(dp) :: v(size(x))
    v=lnd_ichisq_scalar(nu,ssq,x)
  end function lnd_ichisq_vector

  pure function lnd_ichisq_matrix(nu,ssq,x) result(v)
    real(dp), intent(in) :: nu,ssq,x(:,:)
    real(dp) :: v(size(x,1),size(x,2))
    v=lnd_ichisq_scalar(nu,ssq,x)
  end function lnd_ichisq_matrix
end module bayesm_densities
