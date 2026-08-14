program test_optimizers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rmoo
  implicit none
  type(rmoo_real_result) :: r2, r3, rr
  type(rmoo_integer_result) :: rb, ri, rp
  real(dp) :: lo(8), hi(8)
  real(dp), allocatable :: ref3(:,:), refr(:,:)
  integer :: ilo(5), ihi(5), i, j

  lo = 0.0_dp
  hi = 1.0_dp
  call rmoo_optimize_real(zdt1,lo,hi,2,50,35,r2, &
    algorithm=ALG_NSGA2,seed=1101)
  call check_real_result(r2,lo,hi,2,"NSGA-II")
  call assert_true(count(r2%rank==1)>=5,"NSGA-II front size")
  call assert_true(mean_zdt1_gap(r2)<0.45_dp,"NSGA-II ZDT1 convergence")

  call generate_reference_points(3,6,ref3)
  call rmoo_optimize_real(dtlz2,spread(0.0_dp,1,10), &
    spread(1.0_dp,1,10),3,50,20,r3,algorithm=ALG_NSGA3, &
    reference_dirs=ref3,seed=2202)
  call check_real_result(r3,spread(0.0_dp,1,10), &
    spread(1.0_dp,1,10),3,"NSGA-III")
  call assert_true(count(r3%rank==1)>=5,"NSGA-III front size")

  refr = reshape([1.0_dp,0.0_dp,0.75_dp,0.25_dp,0.5_dp,0.5_dp, &
    0.25_dp,0.75_dp,0.0_dp,1.0_dp],[5,2],order=[2,1])
  call rmoo_optimize_real(zdt2,lo,hi,2,40,20,rr, &
    algorithm=ALG_RNSGA2,reference_dirs=refr,epsilon=0.01_dp,seed=3303)
  call check_real_result(rr,lo,hi,2,"R-NSGA-II")

  call rmoo_optimize_binary(binary_obj,12,2,30,8,rb, &
    algorithm=ALG_NSGA2,seed=4404)
  call assert_true(all((rb%population==0).or.(rb%population==1)), &
    "binary population validity")
  call assert_true(all(ieee_is_finite(rb%fitness)),"binary finite fitness")

  ilo = 0
  ihi = 5
  call rmoo_optimize_integer(integer_obj,ilo,ihi,2,30,8,ri, &
    algorithm=ALG_NSGA2,seed=5505)
  call assert_true(all(ri%population>=0.and.ri%population<=5), &
    "integer population validity")

  call rmoo_optimize_permutation(perm_obj,1,7,2,30,8,rp, &
    algorithm=ALG_NSGA2,seed=6606)
  do i=1,size(rp%population,1)
    do j=1,7
      call assert_true(count(rp%population(i,:)==j)==1, &
        "permutation population validity")
    end do
  end do

  print *, "test_optimizers: PASS"
contains
  subroutine zdt1(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    real(dp)::g
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1)
    f(2)=g*(1.0_dp-sqrt(x(1)/g))
  end subroutine zdt1

  subroutine zdt2(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    real(dp)::g
    g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
    f(1)=x(1)
    f(2)=g*(1.0_dp-(x(1)/g)**2)
  end subroutine zdt2

  subroutine dtlz2(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    real(dp)::g,pi
    pi=acos(-1.0_dp)
    g=sum((x(3:)-0.5_dp)**2)
    f(1)=(1.0_dp+g)*cos(0.5_dp*pi*x(1))*cos(0.5_dp*pi*x(2))
    f(2)=(1.0_dp+g)*cos(0.5_dp*pi*x(1))*sin(0.5_dp*pi*x(2))
    f(3)=(1.0_dp+g)*sin(0.5_dp*pi*x(1))
  end subroutine dtlz2

  subroutine binary_obj(x,f)
    integer,intent(in)::x(:)
    real(dp),intent(out)::f(:)
    f(1)=real(sum(x),dp)
    f(2)=real(size(x)-sum(x),dp)
  end subroutine binary_obj

  subroutine integer_obj(x,f)
    integer,intent(in)::x(:)
    real(dp),intent(out)::f(:)
    f(1)=sum(real(x,dp)**2)
    f(2)=sum((real(x,dp)-3.0_dp)**2)
  end subroutine integer_obj

  subroutine perm_obj(x,f)
    integer,intent(in)::x(:)
    real(dp),intent(out)::f(:)
    integer::k,n
    n=size(x)
    f=0.0_dp
    do k=1,n
      f(1)=f(1)+abs(real(x(k)-k,dp))
      f(2)=f(2)+abs(real(x(k)-(n-k+1),dp))
    end do
  end subroutine perm_obj

  subroutine check_real_result(res,l,u,m,name)
    type(rmoo_real_result),intent(in)::res
    real(dp),intent(in)::l(:),u(:)
    integer,intent(in)::m
    character(*),intent(in)::name
    integer::k
    call assert_true(size(res%fitness,2)==m,trim(name)//" objective count")
    call assert_true(all(ieee_is_finite(res%fitness)),trim(name)//" finite fitness")
    do k=1,size(l)
      call assert_true(all(res%population(:,k)>=l(k)-1.0e-12_dp), &
        trim(name)//" lower bounds")
      call assert_true(all(res%population(:,k)<=u(k)+1.0e-12_dp), &
        trim(name)//" upper bounds")
    end do
    call assert_true(all(res%rank>=1),trim(name)//" valid ranks")
  end subroutine check_real_result

  real(dp) function mean_zdt1_gap(res) result(v)
    type(rmoo_real_result),intent(in)::res
    integer::k,n
    v=0.0_dp;n=0
    do k=1,size(res%fitness,1)
      if(res%rank(k)==1)then
        v=v+abs(res%fitness(k,2)-(1.0_dp-sqrt(max(0.0_dp,res%fitness(k,1)))))
        n=n+1
      end if
    end do
    if(n>0)v=v/real(n,dp)
  end function mean_zdt1_gap

  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(*),intent(in)::msg
    if(.not.cond)then
      write(*,*)"FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine assert_true
end program test_optimizers
