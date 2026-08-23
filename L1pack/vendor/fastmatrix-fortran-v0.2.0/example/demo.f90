program demo
  use fastmatrix
  implicit none
  real(dp) :: a(3,3), s(3,3), eval(3), evec(3,3)
  integer :: info
  a = cor_ar1(0.6_dp,3)
  call matrix_sqrt(a,s,info)
  call jacobi_eigen(a,eval,evec,info=info)
  print '(a)', 'AR(1) correlation matrix:'
  print '(3f10.5)', transpose(a)
  print '(a,3f10.5)', 'eigenvalues: ',eval
end program
