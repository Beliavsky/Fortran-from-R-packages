program zdt1_example
  use rmoo
  implicit none
  type(rmoo_real_result) :: result
  real(dp),parameter :: lower(10)=0.0_dp, upper(10)=1.0_dp
  integer :: i

  call rmoo_optimize_real(zdt1,lower,upper,2,60,50,result, &
    algorithm=ALG_NSGA2,seed=1234)

  print '(a,i0)', 'Nondominated solutions: ',count(result%rank==1)
  print '(a)', 'First few nondominated objective pairs:'
  do i=1,size(result%fitness,1)
    if(result%rank(i)==1)then
      print '(2f12.6)',result%fitness(i,:)
    end if
  end do
contains
  subroutine zdt1(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    real(dp)::g
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1)
    f(2)=g*(1.0_dp-sqrt(x(1)/g))
  end subroutine zdt1
end program zdt1_example
