program test_linalg
  use corpcor, only : dp, svd_result, rank_condition_result, fast_svd, pseudoinverse, &
    mpower, rank_condition, is_positive_definite, make_positive_definite
  implicit none
  real(dp), parameter :: tol = 2.0e-9_dp
  real(dp) :: a(3, 2), sym(2, 2), bad(2, 2)
  real(dp), allocatable :: pinv(:, :), recon(:, :), root(:, :), fixed(:, :)
  type(svd_result) :: s
  type(rank_condition_result) :: rc

  a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], [3, 2])
  s = fast_svd(a)
  if (s%rank /= 1) error stop 'rank-deficient SVD rank mismatch'
  pinv = pseudoinverse(a)
  recon = matmul(a, matmul(pinv, a))
  if (maxval(abs(recon - a)) > tol) error stop 'pseudoinverse identity failed'
  rc = rank_condition(a)
  if (rc%rank /= 1) error stop 'rank_condition rank mismatch'

  sym = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
  root = mpower(sym, 0.5_dp)
  if (maxval(abs(matmul(root, root) - sym)) > tol) error stop 'matrix square root failed'
  if (.not. is_positive_definite(sym)) error stop 'SPD test failed'

  bad = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
  if (is_positive_definite(bad)) error stop 'indefinite matrix accepted'
  fixed = make_positive_definite(bad)
  if (.not. is_positive_definite(fixed)) error stop 'positive-definite repair failed'
  print '(a)', 'test_linalg: PASS'
end program test_linalg
