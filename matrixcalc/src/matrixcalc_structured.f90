module matrixcalc_structured
   use r_kinds, only : dp
   implicit none
   private

   type, public :: u_vectors_result
      integer, allocatable :: k(:, :)
      real(dp), allocatable :: ident(:, :)
   end type u_vectors_result

   public :: e_matrices, h_matrices, t_matrices, u_vectors
   public :: k_matrix, commutation_matrix, d_matrix, duplication_matrix
   public :: l_matrix, elimination_matrix, n_matrix

contains

   pure function e_matrices(n) result(e)
      integer, intent(in) :: n !! Matrix order; upstream requires n >= 2.
      real(dp), allocatable :: e(:, :, :, :)
      integer :: i
      integer :: j

      allocate(e(n, n, n, n))
      e = 0.0_dp
      do j = 1, n
         do i = 1, n
            e(i, j, i, j) = 1.0_dp
         end do
      end do
   end function e_matrices

   pure function h_matrices(r, c) result(h)
      integer, intent(in) :: r !! Number of rows in each basis matrix.
      integer, intent(in), optional :: c !! Number of columns; defaults to r.
      real(dp), allocatable :: h(:, :, :, :)
      integer :: cc
      integer :: i
      integer :: j

      cc = r
      if (present(c)) cc = c
      allocate(h(r, cc, r, cc))
      h = 0.0_dp
      do j = 1, cc
         do i = 1, r
            h(i, j, i, j) = 1.0_dp
         end do
      end do
   end function h_matrices

   pure function t_matrices(n) result(tmat)
      integer, intent(in) :: n !! Matrix order; upstream requires n >= 2.
      real(dp), allocatable :: tmat(:, :, :, :)
      integer :: i
      integer :: j

      allocate(tmat(n, n, n, n))
      tmat = 0.0_dp
      do j = 1, n
         do i = 1, n
            tmat(i, j, i, j) = 1.0_dp
            if (i /= j) tmat(j, i, i, j) = 1.0_dp
         end do
      end do
   end function t_matrices

   pure function u_vectors(n) result(u)
      integer, intent(in) :: n !! Matrix order used to enumerate lower-triangular coordinates.
      type(u_vectors_result) :: u
      integer :: p
      integer :: i
      integer :: j
      integer :: idx

      p = n * (n + 1) / 2
      allocate(u%k(n, n), u%ident(p, p))
      u%k = 0
      u%ident = 0.0_dp
      do i = 1, p
         u%ident(i, i) = 1.0_dp
      end do
      do j = 1, n
         do i = j, n
            idx = (j - 1) * n + i - j * (j - 1) / 2
            u%k(i, j) = idx
         end do
      end do
   end function u_vectors

   pure function k_matrix(r, c) result(kmat)
      integer, intent(in) :: r !! Number of rows in the matrix before transposition.
      integer, intent(in), optional :: c !! Number of columns; defaults to r.
      real(dp), allocatable :: kmat(:, :)
      integer :: cc
      integer :: i
      integer :: j
      integer :: input_index
      integer :: output_index

      cc = r
      if (present(c)) cc = c
      allocate(kmat(r * cc, r * cc))
      kmat = 0.0_dp
      do j = 1, cc
         do i = 1, r
            input_index = i + (j - 1) * r
            output_index = j + (i - 1) * cc
            kmat(output_index, input_index) = 1.0_dp
         end do
      end do
   end function k_matrix

   pure function commutation_matrix(r, c) result(kmat)
      integer, intent(in) :: r !! Number of rows in the matrix before transposition.
      integer, intent(in), optional :: c !! Number of columns; defaults to r.
      real(dp), allocatable :: kmat(:, :)

      if (present(c)) then
         kmat = k_matrix(r, c)
      else
         kmat = k_matrix(r)
      end if
   end function commutation_matrix

   pure function d_matrix(n) result(dmat)
      integer, intent(in) :: n !! Order of the symmetric matrix represented by vech.
      real(dp), allocatable :: dmat(:, :)
      integer :: p
      integer :: i
      integer :: j
      integer :: col_index
      integer :: row_ij
      integer :: row_ji

      p = n * (n + 1) / 2
      allocate(dmat(n * n, p))
      dmat = 0.0_dp
      col_index = 0
      do j = 1, n
         do i = j, n
            col_index = col_index + 1
            row_ij = i + (j - 1) * n
            dmat(row_ij, col_index) = 1.0_dp
            if (i /= j) then
               row_ji = j + (i - 1) * n
               dmat(row_ji, col_index) = 1.0_dp
            end if
         end do
      end do
   end function d_matrix

   pure function duplication_matrix(n) result(dmat)
      integer, intent(in) :: n !! Order of the symmetric matrix represented by vech.
      real(dp), allocatable :: dmat(:, :)

      dmat = d_matrix(n)
   end function duplication_matrix

   pure function l_matrix(n) result(lmat)
      integer, intent(in) :: n !! Order of the matrix represented by vec.
      real(dp), allocatable :: lmat(:, :)
      integer :: p
      integer :: i
      integer :: j
      integer :: row_index
      integer :: col_index

      p = n * (n + 1) / 2
      allocate(lmat(p, n * n))
      lmat = 0.0_dp
      row_index = 0
      do j = 1, n
         do i = j, n
            row_index = row_index + 1
            col_index = i + (j - 1) * n
            lmat(row_index, col_index) = 1.0_dp
         end do
      end do
   end function l_matrix

   pure function elimination_matrix(n) result(lmat)
      integer, intent(in) :: n !! Order of the matrix represented by vec.
      real(dp), allocatable :: lmat(:, :)

      lmat = l_matrix(n)
   end function elimination_matrix

   pure function n_matrix(n) result(nmat)
      integer, intent(in) :: n !! Order of the underlying square matrix.
      real(dp), allocatable :: nmat(:, :)
      integer :: i

      nmat = 0.5_dp * k_matrix(n)
      do i = 1, n * n
         nmat(i, i) = nmat(i, i) + 0.5_dp
      end do
   end function n_matrix

end module matrixcalc_structured
