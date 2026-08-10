program mixed_sign
   use nnls, only : dp, nnls_result, nnnpls_fit
   implicit none
   real(dp) :: a(4,3), xtrue(3), b(4), con(3)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp ], [4,3])
   xtrue = [0.5_dp, -1.25_dp, 2.0_dp]
   b = matmul(a, xtrue)
   con = [1.0_dp, -1.0_dp, 1.0_dp]
   call nnnpls_fit(a, b, con, result)
   print '(a,*(f12.6,1x))', 'x = ', result%x
   print '(a,es14.6)', 'RSS = ', result%deviance
end program mixed_sign
