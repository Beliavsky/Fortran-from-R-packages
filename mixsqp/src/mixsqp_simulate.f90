module mixsqp_simulate
  use mixsqp_kinds, only : dp
  use mixsqp_utils, only : logspace, normalize_likelihoods
  implicit none
  private
  public :: simulate_mix_data
contains
  subroutine simulate_mix_data(n,m,x,s,L,simtype,log_output,normalize_rows)
    integer, intent(in) :: n,m
    real(dp), allocatable, intent(out) :: x(:),s(:),L(:,:)
    character(len=*), intent(in), optional :: simtype
    logical, intent(in), optional :: log_output,normalize_rows
    character(len=2) :: st
    logical :: lg,norm
    integer :: k,i,j
    real(dp) :: smax,sd,pi,u1,u2,z
    real(dp), allocatable :: ls(:)
    if (n<1 .or. m<2) error stop "n>=1 and m>=2 required"
    st='n '; if (present(simtype)) st=simtype
    lg=.false.; if (present(log_output)) lg=log_output
    norm=.not.lg; if (present(normalize_rows)) norm=normalize_rows
    if (lg .and. norm) error stop "log output and row normalization are incompatible"
    allocate(x(n),s(m),L(n,m)); pi=acos(-1.0_dp)
    k=n/4
    do i=1,n
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp)); z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
      if (st=='nt' .and. i>n-k) then
        call random_number(u1); u1=max(u1,tiny(1.0_dp))
        x(i)=z/sqrt(-log(u1))
      else if (st=='nt' .and. i>n-2*k) then
        x(i)=5.0_dp*z
      else if (i>n-k) then
        x(i)=6.0_dp*z
      else if (i>n-2*k) then
        x(i)=3.0_dp*z
      else if (st=='nt' .and. i>n-3*k) then
        x(i)=3.0_dp*z
      else
        x(i)=z
      end if
    end do
    if (all(x*x<1.0_dp)) then
      smax=1.0_dp
    else
      smax=2.0_dp*sqrt(maxval(max(x*x-1.0_dp,0.0_dp)))
      if (smax<0.1_dp) smax=0.1_dp
    end if
    s(1)=0.0_dp; ls=logspace(0.1_dp,smax,m-1); s(2:)=ls
    do j=1,m
      sd=sqrt(1.0_dp+s(j)*s(j))
      if (lg) then
        L(:,j)=-0.5_dp*log(2.0_dp*pi)-log(sd)-0.5_dp*(x/sd)**2
      else
        L(:,j)=exp(-0.5_dp*(x/sd)**2)/(sqrt(2.0_dp*pi)*sd)
      end if
    end do
    if (norm) call normalize_likelihoods(L)
  end subroutine simulate_mix_data
end module mixsqp_simulate
