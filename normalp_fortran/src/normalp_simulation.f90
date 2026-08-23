module normalp_simulation
  use normalp_special, only: dp
  use normalp_distribution, only: rnormp
  use normalp_estimation, only: normalp_params, paramp_fit
  use normalp_regression, only: lmp_result, lmp_fit
  implicit none
  private
  public :: simul_mp, simul_lmp, simulation_summary
  type :: simulation_summary
    real(dp), allocatable :: draws(:,:)
    real(dp), allocatable :: means(:), variances(:)
    integer :: convergence_failures = 0
  end type
contains
  subroutine simul_mp(n, m, out, mu, sigmap, p)
    integer, intent(in) :: n, m
    type(simulation_summary), intent(out) :: out
    real(dp), intent(in), optional :: mu, sigmap, p
    real(dp) :: muv, sv, pv
    real(dp), allocatable :: x(:)
    type(normalp_params) :: fit
    integer :: i
    muv=0.0_dp; if(present(mu)) muv=mu
    sv=1.0_dp; if(present(sigmap)) sv=sigmap
    pv=2.0_dp; if(present(p)) pv=p
    allocate(x(n),out%draws(m,5),out%means(5),out%variances(5))
    out%convergence_failures=0
    do i=1,m
      call rnormp(x,muv,sv,pv)
      call paramp_fit(x,fit)
      out%draws(i,:)=[fit%mean,fit%mp,fit%sd,fit%sp,fit%p]
      out%convergence_failures=out%convergence_failures+fit%no_conv
    end do
    call summarize(out%draws,out%means,out%variances)
  end subroutine simul_mp

  subroutine simul_lmp(n, m, beta, out, sigmap, p, estimate_shape)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: beta(:)
    type(simulation_summary), intent(out) :: out
    real(dp), intent(in), optional :: sigmap, p
    logical, intent(in), optional :: estimate_shape
    real(dp) :: sv,pv,u
    real(dp), allocatable :: x(:,:),y(:),e(:)
    type(lmp_result) :: fit
    integer :: i,j,k
    logical :: est
    sv=1.0_dp; if(present(sigmap)) sv=sigmap
    pv=2.0_dp; if(present(p)) pv=p
    est=.true.; if(present(estimate_shape)) est=estimate_shape
    k=size(beta)
    allocate(x(n,k),y(n),e(n),out%draws(m,k+2),out%means(k+2),out%variances(k+2))
    out%convergence_failures=0
    do i=1,m
      x(:,1)=1.0_dp
      do j=2,k
        call fill_uniform(x(:,j))
      end do
      call rnormp(e,0.0_dp,sv,pv)
      y=matmul(x,beta)+e
      if(est) then
        call lmp_fit(x,y,fit)
      else
        call lmp_fit(x,y,fit,p_fixed=pv)
      end if
      out%draws(i,1:k)=fit%coef
      out%draws(i,k+1)=fit%sigma
      out%draws(i,k+2)=fit%p
      out%convergence_failures=out%convergence_failures+fit%no_conv
    end do
    call summarize(out%draws,out%means,out%variances)
  contains
    subroutine fill_uniform(v)
      real(dp),intent(out)::v(:); integer::ii
      do ii=1,size(v); call random_number(u); v(ii)=u; end do
    end subroutine
  end subroutine simul_lmp

  subroutine summarize(a,means,variances)
    real(dp),intent(in)::a(:,:); real(dp),intent(out)::means(:),variances(:)
    integer::j,n
    n=size(a,1)
    do j=1,size(a,2)
      means(j)=sum(a(:,j))/real(n,dp)
      if(n>1) then
        variances(j)=sum((a(:,j)-means(j))**2)/real(n,dp)
      else
        variances(j)=0.0_dp
      end if
    end do
  end subroutine
end module normalp_simulation
