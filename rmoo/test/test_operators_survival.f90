program test_operators_survival
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rmoo
  implicit none
  real(dp) :: p1(3),p2(3),lo(3),hi(3),c1(3),c2(3),mut(3)
  integer :: b1(8),b2(8),bc1(8),bc2(8)
  integer :: q1(6),q2(6),qc1(6),qc2(6),qm(6)
  real(dp) :: f(8,2),ref(3,2),ideal(2),worst(2),nadir(2),smin(2),extreme(2,2)
  integer :: keep(4)

  call ga_seed(12345)
  p1=[0.1_dp,0.3_dp,0.9_dp];p2=[0.8_dp,0.7_dp,0.2_dp]
  lo=0.0_dp;hi=1.0_dp
  call sbx_crossover(p1,p2,lo,hi,c1,c2,20.0_dp,1.0_dp)
  call assert_true(all(c1>=lo.and.c1<=hi),"SBX child 1 bounds")
  call assert_true(all(c2>=lo.and.c2<=hi),"SBX child 2 bounds")
  call polynomial_mutation(c1,lo,hi,mut,20.0_dp,1.0_dp)
  call assert_true(all(mut>=lo.and.mut<=hi),"polynomial mutation bounds")

  b1=[0,0,0,0,1,1,1,1];b2=1-b1
  call hux_crossover(b1,b2,bc1,bc2)
  call assert_true(all((bc1==0).or.(bc1==1)),"HUX binary")
  call assert_true(count(bc1/=b1)==4,"HUX swaps half")

  q1=[1,2,3,4,5,6];q2=[6,5,4,3,2,1]
  call ox_crossover(q1,q2,qc1,qc2)
  call assert_perm(qc1,"OX child 1")
  call assert_perm(qc2,"OX child 2")
  call inversion_mutation(q1,qm)
  call assert_perm(qm,"inversion mutation")

  f = reshape([ &
    0.0_dp,1.0_dp, 0.15_dp,0.75_dp, 0.30_dp,0.58_dp, 0.45_dp,0.42_dp, &
    0.60_dp,0.30_dp, 0.75_dp,0.22_dp, 0.90_dp,0.18_dp, 0.55_dp,0.65_dp], &
    shape(f),order=[2,1])
  call nsga2_survivors(f,4,keep)
  call assert_true(size(keep)==4,"NSGA2 survivor count")
  call assert_true(all(keep>=1.and.keep<=8),"NSGA2 survivor indices")

  ref=reshape([1.0_dp,0.0_dp,0.5_dp,0.5_dp,0.0_dp,1.0_dp],shape(ref),order=[2,1])
  ideal=huge(1.0_dp);worst=-huge(1.0_dp);nadir=0.0_dp
  smin=huge(1.0_dp);extreme=0.0_dp
  call nsga3_survivors(f,4,ref,ideal,worst,smin,extreme,nadir,keep)
  call assert_true(all(keep>=1.and.keep<=8),"NSGA3 survivor indices")
  call assert_true(all(ieee_is_finite(nadir)),"NSGA3 nadir finite")

  call rnsga2_survivors(f,4,ref,0.05_dp,keep)
  call assert_true(all(keep>=1.and.keep<=8),"RNSGA2 survivor indices")

  print *, "test_operators_survival: PASS"
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(*),intent(in)::msg
    if(.not.cond)then
      write(*,*)"FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine
  subroutine assert_perm(x,msg)
    integer,intent(in)::x(:)
    character(*),intent(in)::msg
    integer::j
    do j=1,size(x)
      if(count(x==j)/=1)then
        write(*,*)"FAIL: ",trim(msg),x
        error stop 1
      end if
    end do
  end subroutine
end program test_operators_survival
