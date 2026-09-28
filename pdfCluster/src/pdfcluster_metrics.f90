module pdfcluster_metrics
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use pdfcluster_kinds, only : dp
   implicit none
   private

   public :: adjusted_rand_index

contains

   pure real(dp) function adjusted_rand_index(cl1, cl2) result(ari)
      integer, intent(in) :: cl1(:) !! First clustering labels; values need not be consecutive but must be positive.
      integer, intent(in) :: cl2(:) !! Second clustering labels with the same number of observations as cl1.
      integer, allocatable :: tab(:, :)
      integer, allocatable :: marg1(:)
      integer, allocatable :: marg2(:)
      real(dp) :: den
      real(dp) :: expected
      real(dp) :: num
      real(dp) :: s1
      real(dp) :: s2
      real(dp) :: sij
      integer :: i
      integer :: j
      integer :: n

      ari = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(cl1) /= size(cl2) .or. size(cl1) < 2) return
      if (any(cl1 <= 0) .or. any(cl2 <= 0)) return
      n = size(cl1)
      allocate(tab(maxval(cl1), maxval(cl2)), marg1(maxval(cl1)), marg2(maxval(cl2)))
      tab = 0
      do i = 1, n
         tab(cl1(i), cl2(i)) = tab(cl1(i), cl2(i)) + 1
      end do
      marg1 = sum(tab, dim=2)
      marg2 = sum(tab, dim=1)
      sij = 0.0_dp
      do i = 1, size(tab, 1)
         do j = 1, size(tab, 2)
            sij = sij + choose2(tab(i, j))
         end do
      end do
      s1 = sum([(choose2(marg1(i)), i=1, size(marg1))])
      s2 = sum([(choose2(marg2(j)), j=1, size(marg2))])
      expected = s1*s2/choose2(n)
      num = sij - expected
      den = 0.5_dp*(s1 + s2) - expected
      if (abs(den) > tiny(1.0_dp)) ari = num/den
   end function adjusted_rand_index

   pure real(dp) function choose2(n) result(v)
      integer, intent(in) :: n !! Nonnegative count converted to the number of unordered pairs.
      v = real(n*(n - 1), dp)/2.0_dp
   end function choose2

end module pdfcluster_metrics
