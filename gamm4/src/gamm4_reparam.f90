module gamm4_reparam
   use gamm4_kinds, only : dp
   use gamm4_types, only : gamm4_smooth_t
   use r_linalg, only : symmetric_eigen
   implicit none
   private

   type, public :: smooth_reparam_t
      real(dp), allocatable :: transform(:, :)
      real(dp), allocatable :: fixed_design(:, :)
      real(dp), allocatable :: random_design(:, :)
      real(dp), allocatable :: eigenvalues(:)
      integer :: nullity = 0
      integer :: rank = 0
   end type smooth_reparam_t

   public :: reparameterize_smooth

contains

   subroutine reparameterize_smooth(term, rep, tolerance, status)
      type(gamm4_smooth_t), intent(in) :: term !! Smooth basis and one quadratic penalty to reparameterize.
      type(smooth_reparam_t), intent(out) :: rep !! Demmler-Reinsch-style fixed/random representation and back-transform.
      real(dp), intent(in), optional :: tolerance !! Relative eigenvalue threshold used to identify the penalty null space.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or unsupported multi-penalty smooths.
      real(dp), allocatable :: values(:), vectors(:, :)
      real(dp) :: threshold, tol, vmax
      integer, allocatable :: null_idx(:), pen_idx(:)
      integer :: i, info, k, n, p

      status = 0
      if (.not. allocated(term%basis) .or. .not. allocated(term%spec%penalties)) then
         status = 1
         return
      end if
      n = size(term%basis, 1)
      p = size(term%basis, 2)
      if (n < 1 .or. p < 1) then
         status = 1
         return
      end if
      if (size(term%spec%penalties, 1) /= p .or. size(term%spec%penalties, 2) /= p) then
         status = 2
         return
      end if
      if (size(term%spec%penalties, 3) /= 1) then
         status = 3
         return
      end if

      call symmetric_eigen(term%spec%penalties(:, :, 1), values, vectors, info)
      if (info /= 0) then
         status = 4
         return
      end if
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      vmax = max(1.0_dp, maxval(abs(values)))
      threshold = tol * vmax
      rep%nullity = count(values <= threshold)
      rep%rank = p - rep%nullity
      allocate(null_idx(rep%nullity), pen_idx(rep%rank))
      k = 0
      do i = 1, p
         if (values(i) <= threshold) then
            k = k + 1
            null_idx(k) = i
         end if
      end do
      k = 0
      do i = 1, p
         if (values(i) > threshold) then
            k = k + 1
            pen_idx(k) = i
         end if
      end do

      allocate(rep%transform(p, p), rep%fixed_design(n, rep%nullity), &
         rep%random_design(n, rep%rank), rep%eigenvalues(rep%rank))
      rep%transform = 0.0_dp
      if (rep%nullity > 0) then
         rep%transform(:, 1:rep%nullity) = vectors(:, null_idx)
         rep%fixed_design = matmul(term%basis, rep%transform(:, 1:rep%nullity))
      end if
      if (rep%rank > 0) then
         rep%eigenvalues = values(pen_idx)
         do i = 1, rep%rank
            rep%transform(:, rep%nullity + i) = vectors(:, pen_idx(i)) / sqrt(rep%eigenvalues(i))
         end do
         rep%random_design = matmul(term%basis, rep%transform(:, rep%nullity + 1:p))
      end if
   end subroutine reparameterize_smooth

end module gamm4_reparam
