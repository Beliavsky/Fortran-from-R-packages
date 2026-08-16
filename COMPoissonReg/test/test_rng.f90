program test_rng
   use compoissonreg
   implicit none
   integer,parameter :: n=20000
   integer :: x(n),fails
   real(dp) :: m,zr
   fails=0
   call set_seed(24681357)
   call rcmp(n,4.0_dp,1.0_dp,x)
   m=sum(real(x,dp))/real(n,dp)
   if(abs(m-4.0_dp)>0.08_dp)then;print *,'FAIL CMP RNG mean',m;fails=fails+1;end if
   call rzicmp(n,3.0_dp,1.0_dp,0.2_dp,x)
   m=sum(real(x,dp))/real(n,dp)
   zr=real(count(x==0),dp)/real(n,dp)
   if(abs(m-2.4_dp)>0.08_dp)then;print *,'FAIL ZICMP RNG mean',m;fails=fails+1;end if
   if(abs(zr-(0.2_dp+0.8_dp*exp(-3.0_dp)))>0.02_dp)then;print *,'FAIL ZICMP zeros',zr;fails=fails+1;end if
   if(fails==0)then;print *,'test_rng: PASS';else;error stop 1;end if
contains
   subroutine set_seed(v)
      integer,intent(in)::v;integer::nseed,i;integer,allocatable::seed(:)
      call random_seed(size=nseed);allocate(seed(nseed));do i=1,nseed;seed(i)=v+37*i;end do;call random_seed(put=seed)
   end subroutine set_seed
end program test_rng
