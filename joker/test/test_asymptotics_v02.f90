program test_asymptotics_v02
  use joker
  implicit none
  integer :: fails
  real(dp) :: d2(2,2), d3(3,3)
  real(dp), allocatable :: d(:,:)
  real(dp) :: alpha(3), shape(3), prob(3)
  fails = 0

  call check_close(avar_bern_mle(0.3_dp), 0.21_dp, 1e-12_dp, "Bern MLE")
  call check_close(avar_bern_me(0.3_dp), 0.21_dp, 1e-12_dp, "Bern ME")
  call check_close(avar_binom_mle(5,0.3_dp), 0.042_dp, 1e-12_dp, "Binom MLE")
  call check_close(avar_geom_mle(0.4_dp), 0.096_dp, 1e-12_dp, "Geom MLE")
  call check_close(avar_exp_mle(2.0_dp), 4.0_dp, 1e-12_dp, "Exp MLE")
  call check_close(avar_pois_mle(3.2_dp), 3.2_dp, 1e-12_dp, "Pois MLE")
  call check_close(avar_nbinom_mle(4.0_dp,0.6_dp), 0.036_dp, 1e-12_dp, "Nbinom MLE")
  call check_close(avar_chisq_me(7.0_dp), 14.0_dp, 1e-12_dp, "Chisq ME")

  d2 = avar_beta_mle(2.3_dp,4.1_dp)
  call check_mat2(d2, reshape([9.42439720_dp,14.90048477_dp, &
    14.90048477_dp,32.90886240_dp],[2,2]), 2e-5_dp, "Beta MLE")
  d2 = avar_beta_me(2.3_dp,4.1_dp)
  call check_mat2(d2, reshape([10.45657431_dp,16.46546291_dp, &
    16.46546291_dp,35.34214380_dp],[2,2]), 2e-7_dp, "Beta ME")
  d2 = avar_beta_same(2.3_dp,4.1_dp)
  call check_mat2(d2, reshape([9.53383144_dp,15.00590164_dp, &
    15.00590164_dp,33.07068482_dp],[2,2]), 2e-5_dp, "Beta SAME")

  d2 = avar_gamma_mle(2.7_dp,1.4_dp)
  call check_mat2(d2, reshape([13.01376751_dp,-6.74787945_dp, &
    -6.74787945_dp,4.22482638_dp],[2,2]), 2e-5_dp, "Gamma MLE")
  d2 = avar_gamma_me(2.7_dp,1.4_dp)
  call check_mat2(d2, reshape([19.98_dp,-10.36_dp,-10.36_dp, &
    6.09777777777778_dp],[2,2]), 1e-10_dp, "Gamma ME")
  d2 = avar_gamma_same(2.7_dp,1.4_dp)
  call check_mat2(d2, reshape([13.39247515_dp,-6.94424638_dp, &
    -6.94424638_dp,4.32664627_dp],[2,2]), 2e-5_dp, "Gamma SAME")

  ! These two use corrected scale powers; upstream tests used scale = 1 and
  ! therefore did not expose the source-level Fisher-information typo.
  d2 = avar_laplace_mle(2.0_dp)
  call check_mat2(d2, reshape([4.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2]), &
    1e-12_dp, "Laplace corrected scale")
  d2 = avar_lnorm_mle(2.0_dp)
  call check_mat2(d2, reshape([4.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]), &
    1e-12_dp, "Lnorm corrected scale")

  prob = [0.7_dp,0.2_dp,0.1_dp]
  d = avar_cat_mle(prob)
  d2 = d
  call check_mat2(d2, reshape([0.21_dp,-0.14_dp,-0.14_dp,0.16_dp], &
    [2,2]), 1e-12_dp, "Categorical covariance")
  d = avar_multinom_mle(1000,prob)
  d2 = d
  call check_mat2(d2, reshape([0.00021_dp,-0.00014_dp,-0.00014_dp, &
    0.00016_dp],[2,2]), 1e-14_dp, "Multinomial covariance")

  alpha = [1.7_dp,2.4_dp,3.2_dp]
  d = avar_dir_mle(alpha)
  d3 = d
  call check_mat3(d3, reshape([3.06562268_dp,2.77927902_dp,3.90846228_dp, &
    2.77927902_dp,6.22071193_dp,6.01825779_dp, &
    3.90846228_dp,6.01825779_dp,11.19323889_dp],[3,3]), &
    3e-5_dp, "Dirichlet MLE")
  d = avar_dir_me(alpha)
  d3 = d
  call check_mat3(d3, reshape([3.97545748_dp,3.63190476_dp,5.04114934_dp, &
    3.63190476_dp,7.42224434_dp,7.36225805_dp, &
    5.04114934_dp,7.36225805_dp,13.00465536_dp],[3,3]), &
    3e-7_dp, "Dirichlet ME")
  d = avar_dir_same(alpha)
  d3 = d
  call check_mat3(d3, reshape([3.22216755_dp,2.74593038_dp,3.91043729_dp, &
    2.74593038_dp,6.42203637_dp,6.10006456_dp, &
    3.91043729_dp,6.10006456_dp,11.41695354_dp],[3,3]), &
    3e-5_dp, "Dirichlet SAME")

  shape = [1.4_dp,2.2_dp,3.1_dp]
  d = avar_multigam_mle(shape,1.3_dp)
  call check_sym(d,1e-11_dp,"Multigamma MLE symmetry")
  call check_close(d(4,4),1.25333147_dp,3e-5_dp,"Multigamma MLE scale")
  d = avar_multigam_me(shape,1.3_dp)
  call check_sym(d,1e-11_dp,"Multigamma ME symmetry")
  call check_close(d(1,1),2.74772239_dp,3e-6_dp,"Multigamma ME 11")
  call check_close(d(4,4),1.96682866_dp,3e-6_dp,"Multigamma ME scale")
  d = avar_multigam_same(shape,1.3_dp)
  call check_sym(d,1e-11_dp,"Multigamma SAME symmetry")
  call check_close(d(1,1),1.96372498_dp,3e-5_dp,"Multigamma SAME 11")
  call check_close(d(4,4),1.29083089_dp,3e-5_dp,"Multigamma SAME scale")

  ! Fisher-information inversion checks for the matrix-valued MLE families.
  d2 = avar_beta_mle(2.3_dp,4.1_dp)
  call check_identity(matmul(finf_beta(2.3_dp,4.1_dp),d2), 3e-11_dp, "Beta information")
  d = avar_cat_mle(prob)
  call check_identity(matmul(finf_cat(prob),d), 3e-12_dp, "Cat information")
  d = avar_dir_mle(alpha)
  call check_identity(matmul(finf_dir(alpha),d), 5e-10_dp, "Dir information")
  d2 = avar_gamma_mle(2.7_dp,1.4_dp)
  call check_identity(matmul(finf_gamma(2.7_dp,1.4_dp),d2), 3e-11_dp, "Gamma information")
  d = avar_multigam_mle(shape,1.3_dp)
  call check_identity(matmul(finf_multigam(shape,1.3_dp),d), 5e-10_dp, "Multigamma information")
  d = avar_multinom_mle(1000,prob)
  call check_identity(matmul(finf_multinom(1000,prob),d), 3e-12_dp, "Multinom information")

  d2 = avar_cauchy_mle(1.5_dp)
  call check_mat2(d2, reshape([4.5_dp,0.0_dp,0.0_dp,4.5_dp],[2,2]), &
    1e-12_dp, "Cauchy MLE")
  d2 = avar_norm_mle(2.0_dp)
  call check_mat2(d2, reshape([4.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]), &
    1e-12_dp, "Normal MLE")
  d2 = avar_norm_me(2.0_dp)
  call check_mat2(d2, avar_norm_mle(2.0_dp), 1e-12_dp, "Normal ME")

  if (fails == 0) then
    print '(a)', "test_asymptotics_v02: PASS"
  else
    print '(a,i0)', "test_asymptotics_v02: FAIL ", fails
    error stop 1
  end if
contains
  subroutine check_close(a,b,tol,name)
    real(dp),intent(in)::a,b,tol
    character(*),intent(in)::name
    if(abs(a-b)>tol)then
      print '(a,2es18.8)',trim(name)//" mismatch: ",a,b
      fails=fails+1
    end if
  end subroutine check_close

  subroutine check_mat2(a,b,tol,name)
    real(dp),intent(in)::a(2,2),b(2,2),tol
    character(*),intent(in)::name
    if(maxval(abs(a-b))>tol)then
      print '(a,es18.8)',trim(name)//" max error: ",maxval(abs(a-b))
      fails=fails+1
    end if
  end subroutine check_mat2

  subroutine check_mat3(a,b,tol,name)
    real(dp),intent(in)::a(3,3),b(3,3),tol
    character(*),intent(in)::name
    if(maxval(abs(a-b))>tol)then
      print '(a,es18.8)',trim(name)//" max error: ",maxval(abs(a-b))
      fails=fails+1
    end if
  end subroutine check_mat3

  subroutine check_identity(a,tol,name)
    real(dp),intent(in)::a(:,:),tol
    character(*),intent(in)::name
    real(dp),allocatable::eye(:,:)
    integer::i,n
    n=size(a,1)
    allocate(eye(n,n))
    eye=0.0_dp
    do i=1,n
      eye(i,i)=1.0_dp
    end do
    if(maxval(abs(a-eye))>tol)then
      print '(a,es18.8)',trim(name)//" max identity error: ",maxval(abs(a-eye))
      fails=fails+1
    end if
  end subroutine check_identity

  subroutine check_sym(a,tol,name)
    real(dp),intent(in)::a(:,:),tol
    character(*),intent(in)::name
    if(maxval(abs(a-transpose(a)))>tol)then
      print '(a,es18.8)',trim(name)//" max asymmetry: ", &
        maxval(abs(a-transpose(a)))
      fails=fails+1
    end if
  end subroutine check_sym
end program test_asymptotics_v02
