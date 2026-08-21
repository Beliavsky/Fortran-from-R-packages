program test_principal_balances
  use compositions
  implicit none
  real(dp) :: x(8,4)
  real(dp), allocatable :: vh(:,:),va(:,:),vm(:,:)
  x=reshape([ &
   0.55_dp,0.25_dp,0.12_dp,0.08_dp, 0.50_dp,0.28_dp,0.14_dp,0.08_dp, &
   0.58_dp,0.22_dp,0.12_dp,0.08_dp, 0.48_dp,0.30_dp,0.14_dp,0.08_dp, &
   0.15_dp,0.15_dp,0.35_dp,0.35_dp, 0.14_dp,0.16_dp,0.34_dp,0.36_dp, &
   0.16_dp,0.14_dp,0.36_dp,0.34_dp, 0.15_dp,0.15_dp,0.34_dp,0.36_dp],[8,4],order=[2,1])
  vh=principal_balance_hclust(x); va=principal_balance_angprox(x); vm=principal_balance_maxvar(x)
  call check_basis(vh,'PBhclust'); call check_basis(va,'PBangprox'); call check_basis(vm,'PBmaxvar')
  print *, 'test_principal_balances: PASS'
contains
  subroutine check_basis(v,name)
    real(dp), intent(in) :: v(:,:)
    character(len=*), intent(in) :: name
    real(dp) :: g(size(v,2),size(v,2))
    integer :: q
    if(any(shape(v)/=[4,3])) then; print *,name,' shape'; error stop 1; end if
    g=matmul(transpose(v),v)
    do q=1,3; g(q,q)=g(q,q)-1.0_dp; end do
    if(maxval(abs(g))>1.0e-9_dp) then; print *,name,' orthogonality'; error stop 1; end if
    if(maxval(abs(sum(v,dim=1)))>1.0e-9_dp) then; print *,name,' clr subspace'; error stop 1; end if
  end subroutine
end program
