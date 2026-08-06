module lme4_quadrature
   use lme4_kinds, only : dp, pi
   use lme4_types, only : gh_rule_t
   use lme4_linalg, only : jacobi_eigen
   implicit none
   private
   public :: gh_rule, gh_integrate

   abstract interface
      function scalar_function(x) result(y)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   subroutine gh_rule(order, rule, info)
      integer, intent(in) :: order
      type(gh_rule_t), intent(out) :: rule
      integer, intent(out) :: info
      real(dp), allocatable :: jmat(:,:), roots(:), vectors(:,:)
      integer :: i

      if (order < 0 .or. order > 200) then
         info = -1
         allocate(rule%nodes(0),rule%weights(0),rule%log_density(0))
         return
      end if
      if (order == 0) then
         info = 0
         allocate(rule%nodes(0),rule%weights(0),rule%log_density(0))
         return
      end if
      allocate(jmat(order,order))
      jmat = 0.0_dp
      do i = 1, order - 1
         jmat(i,i+1) = sqrt(real(i,dp)/2.0_dp)
         jmat(i+1,i) = jmat(i,i+1)
      end do
      call jacobi_eigen(jmat, roots, vectors, info, tolerance=1.0e-14_dp, max_sweeps=200*order*order)
      if (info /= 0) return
      allocate(rule%nodes(order),rule%weights(order),rule%log_density(order))
      do i = 1, order
         rule%nodes(i) = sqrt(2.0_dp)*roots(order-i+1)
         rule%weights(i) = vectors(1,order-i+1)**2
         rule%log_density(i) = -0.5_dp*(log(2.0_dp*pi)+rule%nodes(i)**2)
      end do
      rule%weights = rule%weights/sum(rule%weights)
   end subroutine gh_rule

   real(dp) function gh_integrate(fun, rule) result(value)
      procedure(scalar_function) :: fun
      type(gh_rule_t), intent(in) :: rule
      integer :: i
      value = 0.0_dp
      do i = 1, size(rule%nodes)
         value = value + rule%weights(i)*fun(rule%nodes(i))
      end do
   end function gh_integrate

end module lme4_quadrature
