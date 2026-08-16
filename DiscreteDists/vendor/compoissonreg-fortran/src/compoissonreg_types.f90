module compoissonreg_types
   use compoissonreg_kinds, only : dp
   implicit none
   private

   type, public :: cmp_control_t
      integer :: ymax = 1000000
      integer :: max_iter = 150
      real(dp) :: hybrid_tol = 1.0e-2_dp
      real(dp) :: truncate_tol = 1.0e-6_dp
      real(dp) :: optim_tol = 1.0e-5_dp
      real(dp) :: fd_step = 1.0e-5_dp
   end type cmp_control_t

   type, public :: cmp_init_t
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: zeta(:)
   end type cmp_init_t

   type, public :: cmp_fixed_t
      logical, allocatable :: beta(:)
      logical, allocatable :: gamma(:)
      logical, allocatable :: zeta(:)
   end type cmp_fixed_t

   type, public :: cmp_offset_t
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: w(:)
   end type cmp_offset_t

   type, public :: cmp_fit_t
      integer, allocatable :: y(:)
      real(dp), allocatable :: xmat(:,:), smat(:,:)
      real(dp), allocatable :: beta(:), gamma(:)
      real(dp), allocatable :: hessian(:,:), covariance(:,:)
      type(cmp_offset_t) :: offset
      type(cmp_fixed_t) :: fixed
      type(cmp_control_t) :: control
      real(dp) :: loglik = -huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0
   end type cmp_fit_t

   type, public :: zicmp_fit_t
      integer, allocatable :: y(:)
      real(dp), allocatable :: xmat(:,:), smat(:,:), wmat(:,:)
      real(dp), allocatable :: beta(:), gamma(:), zeta(:)
      real(dp), allocatable :: hessian(:,:), covariance(:,:)
      type(cmp_offset_t) :: offset
      type(cmp_fixed_t) :: fixed
      type(cmp_control_t) :: control
      real(dp) :: loglik = -huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0
   end type zicmp_fit_t

   type, public :: equitest_t
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: df = 0
   end type equitest_t

   public :: default_init, default_fixed, default_offset

contains

   function default_init(d1, d2, d3) result(out)
      integer, intent(in) :: d1, d2
      integer, intent(in), optional :: d3
      type(cmp_init_t) :: out
      integer :: k3
      k3 = 0
      if (present(d3)) k3 = d3
      allocate(out%beta(d1), out%gamma(d2), out%zeta(k3))
      out%beta = 0.0_dp
      out%gamma = 0.0_dp
      if (k3 > 0) out%zeta = 0.0_dp
   end function default_init

   function default_fixed(d1, d2, d3) result(out)
      integer, intent(in) :: d1, d2
      integer, intent(in), optional :: d3
      type(cmp_fixed_t) :: out
      integer :: k3
      k3 = 0
      if (present(d3)) k3 = d3
      allocate(out%beta(d1), out%gamma(d2), out%zeta(k3))
      out%beta = .false.
      out%gamma = .false.
      if (k3 > 0) out%zeta = .false.
   end function default_fixed

   function default_offset(n) result(out)
      integer, intent(in) :: n
      type(cmp_offset_t) :: out
      allocate(out%x(n), out%s(n), out%w(n))
      out%x = 0.0_dp
      out%s = 0.0_dp
      out%w = 0.0_dp
   end function default_offset

end module compoissonreg_types
