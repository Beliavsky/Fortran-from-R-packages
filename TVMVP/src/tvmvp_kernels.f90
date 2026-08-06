! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_kernels
  use tvmvp_kinds, only : dp
  implicit none
  private
  public :: kernel_function, epanechnikov_kernel, boundary_kernel
  public :: two_fold_convolution_kernel, integrate_kernel

  abstract interface
    pure real(dp) function kernel_function(u)
      import dp
      real(dp), intent(in) :: u
    end function kernel_function
  end interface
contains
  pure real(dp) function epanechnikov_kernel(u)
    real(dp), intent(in) :: u
    if (abs(u) <= 1.0_dp) then
      epanechnikov_kernel = 0.75_dp*(1.0_dp-u*u)
    else
      epanechnikov_kernel = 0.0_dp
    end if
  end function epanechnikov_kernel

  pure real(dp) function epanechnikov_integral(a,b)
    real(dp), intent(in) :: a,b
    real(dp) :: lo,hi
    lo=max(-1.0_dp,min(1.0_dp,a))
    hi=max(-1.0_dp,min(1.0_dp,b))
    if (hi<=lo) then
      epanechnikov_integral=0.0_dp
    else
      epanechnikov_integral=0.75_dp*((hi-hi**3/3.0_dp)-(lo-lo**3/3.0_dp))
    end if
  end function epanechnikov_integral

  real(dp) function integrate_kernel(kernel,a,b,n_intervals)
    procedure(kernel_function) :: kernel
    real(dp), intent(in) :: a,b
    integer, intent(in), optional :: n_intervals
    integer :: n,i
    real(dp) :: h,s,x
    n=400
    if (present(n_intervals)) n=max(20,n_intervals)
    if (mod(n,2)/=0) n=n+1
    if (b<=a) then
      integrate_kernel=0.0_dp
      return
    end if
    h=(b-a)/real(n,dp)
    s=kernel(a)+kernel(b)
    do i=1,n-1
      x=a+h*real(i,dp)
      if (mod(i,2)==0) then
        s=s+2.0_dp*kernel(x)
      else
        s=s+4.0_dp*kernel(x)
      end if
    end do
    integrate_kernel=s*h/3.0_dp
  end function integrate_kernel

  real(dp) function boundary_kernel(target, observation, n_time, bandwidth, kernel, source_compatible)
    integer, intent(in) :: target, observation, n_time
    real(dp), intent(in) :: bandwidth
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible
    logical :: compat
    integer :: boundary_index, th_floor
    real(dp) :: scaled_diff,kval,denom,lo,hi
    compat=.true.
    if (present(source_compatible)) compat=source_compatible
    if (bandwidth<=0.0_dp .or. n_time<1) then
      boundary_kernel=0.0_dp
      return
    end if
    scaled_diff=real(target-observation,dp)/(real(n_time,dp)*bandwidth)
    if (present(kernel)) then
      kval=kernel(scaled_diff)/bandwidth
    else
      kval=epanechnikov_kernel(scaled_diff)/bandwidth
    end if
    if (compat) then
      boundary_index=observation
    else
      boundary_index=target
    end if
    th_floor=floor(real(n_time,dp)*bandwidth)
    if (boundary_index<th_floor) then
      lo=-real(boundary_index,dp)/(real(n_time,dp)*bandwidth); hi=1.0_dp
      if (present(kernel)) then
        denom=integrate_kernel(kernel,lo,hi)
      else
        denom=epanechnikov_integral(lo,hi)
      end if
    else if (boundary_index>(n_time-th_floor)) then
      lo=-1.0_dp; hi=(1.0_dp-real(boundary_index,dp)/real(n_time,dp))/bandwidth
      if (present(kernel)) then
        denom=integrate_kernel(kernel,lo,hi)
      else
        denom=epanechnikov_integral(lo,hi)
      end if
    else
      denom=1.0_dp
    end if
    if (abs(denom)<=tiny(1.0_dp)) then
      boundary_kernel=0.0_dp
    else
      boundary_kernel=kval/denom
    end if
  end function boundary_kernel

  real(dp) function two_fold_convolution_kernel(u,kernel)
    real(dp), intent(in) :: u
    procedure(kernel_function), optional :: kernel
    real(dp) :: a
    integer :: n,i
    real(dp) :: h,s,v
    if (abs(u)>2.0_dp) then
      two_fold_convolution_kernel=0.0_dp
      return
    end if
    if (.not.present(kernel)) then
      a=abs(u)
      two_fold_convolution_kernel=0.6_dp-0.75_dp*a*a+0.375_dp*a**3-0.01875_dp*a**5
      return
    end if
    n=400; h=2.0_dp/real(n,dp)
    s=kernel(-1.0_dp)*kernel(u+1.0_dp)+kernel(1.0_dp)*kernel(u-1.0_dp)
    do i=1,n-1
      v=-1.0_dp+h*real(i,dp)
      if (mod(i,2)==0) then
        s=s+2.0_dp*kernel(v)*kernel(u-v)
      else
        s=s+4.0_dp*kernel(v)*kernel(u-v)
      end if
    end do
    two_fold_convolution_kernel=s*h/3.0_dp
  end function two_fold_convolution_kernel
end module tvmvp_kernels
