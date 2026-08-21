! SPDX-License-Identifier: GPL-2.0-or-later
program test_stodes
  use rootsolve, only : dp, stodes, steady_1d, sparse_options, steady_options, steady_result
  implicit none
  integer,parameter::n=80
  real(dp)::y0(n)
  type(sparse_options)::sp
  type(steady_options)::op
  type(steady_result)::a,b
  integer::i
  do i=1,n;y0(i)=real(mod(7*i,19),dp)/10.0_dp;end do
  sp%base%rtol=1e-9_dp;sp%base%atol=1e-11_dp
  sp%sparsetype='1D';sp%nspec=1;sp%dims=[n,0,0]
  a=stodes(rhs,y0,options=sp)
  if(.not.a%steady.or.maxval(abs(a%y-1.0_dp))>2e-8_dp)error stop 1
  if(size(a%colind) /= 3*n-2) error stop 2
  op=sp%base
  b=steady_1d(rhs,y0,1,dimens=n,method='stodes',options=op)
  if(.not.b%steady.or.maxval(abs(b%y-1.0_dp))>2e-8_dp)error stop 3
  print *, 'test_stodes: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    integer::k,nn
    nn=size(y)
    dy(1)=2*y(1)-y(2)-1
    do k=2,nn-1
      dy(k)=-y(k-1)+2*y(k)-y(k+1)
    end do
    dy(nn)=-y(nn-1)+2*y(nn)-1
    if(t < -huge(1.0_dp))error stop 99
  end subroutine rhs
end program test_stodes
