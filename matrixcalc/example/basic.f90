program basic
   use matrixcalc, only : dp, hilbert_matrix, k_matrix, pascal_matrix, vec
   implicit none

   real(dp) :: a(2, 3)
   real(dp), allocatable :: h(:, :)
   real(dp), allocatable :: k(:, :)
   real(dp), allocatable :: p(:, :)
   real(dp), allocatable :: transposed_vec(:)

   h = hilbert_matrix(4)
   p = pascal_matrix(4)
   a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2, 3])
   k = k_matrix(2, 3)
   transposed_vec = matmul(k, vec(a))

   print '(a,1x,f10.7)', 'Hilbert(4,4):', h(4, 4)
   print '(a,1x,f6.1)', 'Pascal(4,2):', p(4, 2)
   print '(a,6(1x,f4.1))', 'K*vec(A):', transposed_vec
end program basic
