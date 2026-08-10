program basic_nnls
   use nnls, only : dp, nnls_result, nnls_fit
   implicit none
   real(dp) :: a(4,3), b(4)
   type(nnls_result) :: result

   a = reshape([ &
      1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp ], [4,3])
   b = [1.5_dp, 2.0_dp, 1.0_dp, 2.5_dp]
   call nnls_fit(a, b, result)
   print '(a,*(f12.6,1x))', 'x = ', result%x
   print '(a,es14.6)', 'RSS = ', result%deviance
end program basic_nnls
