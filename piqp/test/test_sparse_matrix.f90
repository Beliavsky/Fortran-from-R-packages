program test_sparse_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use piqp
   use matrix_sparse, only : csr_matrix, csc_matrix, csr_from_dense, csc_from_csr
   implicit none
   real(dp) :: pd(2,2), ad(1,2), gd(2,2), c(2), b(1), hu(2), xl(2), xu(2), inf
   type(csr_matrix) :: pr, ar, gr
   type(csc_matrix) :: p, a, g
   type(piqp_result_type) :: r
   inf=ieee_value(0.0_dp,ieee_positive_inf)
   pd=reshape([6.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
   ad=reshape([1.0_dp,-2.0_dp],[1,2])
   gd=reshape([1.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,2])
   call csr_from_dense(pd,pr); call csc_from_csr(pr,p)
   call csr_from_dense(ad,ar); call csc_from_csr(ar,a)
   call csr_from_dense(gd,gr); call csc_from_csr(gr,g)
   c=[-1.0_dp,-4.0_dp]; b=0.0_dp; hu=[1.0_dp,1.0_dp]; xl=[-inf,-1.0_dp]; xu=[inf,1.0_dp]
   call solve_piqp(p,c,r,a,b,g,h_u=hu,x_l=xl,x_u=xu)
   if(r%info%status/=PIQP_SOLVED) error stop 'sparse status'
   if(maxval(abs(r%x-[0.4285714_dp,0.2142857_dp]))>1.0e-6_dp) error stop 'sparse x'
   print *, 'test_sparse_matrix: PASS'
end program test_sparse_matrix
