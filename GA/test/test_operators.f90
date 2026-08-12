program test_operators
  use ga
  implicit none
  integer :: p1(10), p2(10), c1(10), c2(10), m(10), s, i
  integer :: b1(12), b2(12), bc1(12), bc2(12)
  real(dp) :: r1(4), r2(4), rc1(4), rc2(4), rm(4), lo(4), up(4), f(8)
  integer :: sel(8)
  p1=[(i,i=1,10)];p2=p1(10:1:-1)
  call ga_seed(1234)
  do s=CROSS_PERM_CYCLE,CROSS_PERM_PBX
    call crossover_int(p1,p2,s,c1,c2)
    call check_perm(c1);call check_perm(c2)
  end do
  do s=MUT_PERM_INVERSION,MUT_PERM_SCRAMBLE
    call mutation_permutation(p1,s,m);call check_perm(m)
  end do
  b1=0;b2=1
  call crossover_int(b1,b2,CROSS_BINARY_UNIFORM,bc1,bc2)
  if(any((bc1/=0).and.(bc1/=1))) error stop "binary crossover"
  call mutation_binary(b1,bc1)
  if(count(b1/=bc1)/=1) error stop "binary mutation"
  lo=-2.0_dp;up=2.0_dp;r1=[-1.5_dp,-0.5_dp,0.5_dp,1.5_dp];r2=-r1
  do s=CROSS_REAL_WEIGHTED,CROSS_REAL_LAPLACE
    call crossover_real(r1,r2,lo,up,s,rc1,rc2)
    if(any(rc1<lo).or.any(rc1>up).or.any(rc2<lo).or.any(rc2>up)) error stop "real crossover bounds"
  end do
  do s=MUT_REAL_RANDOM,MUT_REAL_POWER
    call mutation_real(r1,lo,up,s,10,100,rm)
    if(any(rm<lo).or.any(rm>up)) error stop "real mutation bounds"
  end do
  f=[1.0_dp,4.0_dp,2.0_dp,8.0_dp,3.0_dp,7.0_dp,5.0_dp,6.0_dp]
  do s=SEL_LINEAR_RANK,SEL_REAL_SIGMA
    call select_indices(f,s,sel)
    if(any(sel<1).or.any(sel>8)) error stop "selection indices"
  end do
  print *, "test_operators: PASS"
contains
  subroutine check_perm(x)
    integer,intent(in)::x(:)
    integer::j
    do j=1,size(x)
      if(count(x==j)/=1) error stop "invalid permutation"
    end do
  end subroutine check_perm
end program test_operators
