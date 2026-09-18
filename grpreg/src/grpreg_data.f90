module grpreg_data
   use grpreg_kinds, only : dp
   implicit none
   private
   public :: gen_nonlinear_data

contains

   subroutine gen_nonlinear_data(n, p, x, y, mu, seed)
      integer, intent(in) :: n !! Positive sample size for the nonlinear sparse-additive example.
      integer, intent(in) :: p !! Number of predictors; at least six are required by the upstream generator.
      real(dp), allocatable, intent(out) :: x(:,:) !! Uniform[0,1] predictor matrix with n rows and p columns.
      real(dp), allocatable, intent(out) :: y(:) !! Gaussian response with standard deviation 0.25 about mu.
      real(dp), allocatable, intent(out) :: mu(:) !! Conditional response mean formed from the first six nonlinear signals.
      integer, optional, intent(in) :: seed !! Optional deterministic Fortran RNG seed; streams need not match R.
      integer, allocatable :: old_seed(:), new_seed(:)
      real(dp), allocatable :: eta(:,:)
      real(dp) :: u1, u2, z, den
      integer :: i, j, ns
      if (n < 1) error stop 'gen_nonlinear_data: n must be positive'
      if (p < 6) error stop 'gen_nonlinear_data: p must be at least six'
      call random_seed(size=ns)
      allocate(old_seed(ns))
      call random_seed(get=old_seed)
      if (present(seed)) then
         allocate(new_seed(ns))
         do i = 1, ns
            new_seed(i) = modulo(abs(seed)+104729*i,huge(1)-1)+1
         end do
         call random_seed(put=new_seed)
      end if
      allocate(x(n,p),y(n),mu(n),eta(n,6))
      call random_number(x)
      den = 1.0_dp-exp(-10.0_dp)
      eta(:,1) = 2.0_dp*(exp(-10.0_dp*x(:,1))-exp(-10.0_dp))/den-1.0_dp
      eta(:,2) = -2.0_dp*(exp(-10.0_dp*x(:,2))-exp(-10.0_dp))/den+1.0_dp
      eta(:,3) = 2.0_dp*x(:,3)-1.0_dp
      eta(:,4) = -2.0_dp*x(:,4)+1.0_dp
      eta(:,5) = 8.0_dp*(x(:,5)-0.5_dp)**2-1.0_dp
      eta(:,6) = -8.0_dp*(x(:,6)-0.5_dp)**2+1.0_dp
      mu = 5.0_dp
      do j = 1, 6
         mu = mu+eta(:,j)
      end do
      i = 1
      do while (i <= n)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1,tiny(1.0_dp))
         z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
         y(i) = mu(i)+0.25_dp*z
         i = i+1
         if (i <= n) then
            z = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2)
            y(i) = mu(i)+0.25_dp*z
            i = i+1
         end if
      end do
      if (present(seed)) call random_seed(put=old_seed)
   end subroutine gen_nonlinear_data

end module grpreg_data
