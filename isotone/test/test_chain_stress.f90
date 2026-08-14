program test_chain_stress
   use isotone
   use iso_fortran_env, only : int64
   implicit none
   integer, parameter :: n=20, ncase=100
   integer :: a(n-1,2), i, ic
   integer(int64) :: state
   real(dp) :: y(n),w(n),u
   type(gpava_result) :: p
   type(active_set_options) :: o
   type(active_set_result) :: r
   do i=1,n-1; a(i,:)=[i,i+1]; end do
   allocate(o%weights(n));o%solver=ISO_LS
   state=13579_int64
   do ic=1,ncase
      do i=1,n
         call rng_u(state,u); y(i)=6.0_dp*(u-0.5_dp)
         call rng_u(state,u); w(i)=0.25_dp+2.0_dp*u
      end do
      o%weights=w
      call gpava_fit(y,p,weights=w)
      call active_set(a,y,r,o,maxiter=200)
      if(r%status/=ISO_SUCCESS) error stop 'chain stress active-set failure'
      if(maxval(abs(p%x-r%x))>2.0e-10_dp) error stop 'chain stress mismatch'
   end do
   print *, 'test_chain_stress: PASS'
contains
   subroutine rng_u(s,u)
      integer(int64),intent(inout)::s
      real(dp),intent(out)::u
      s=modulo(16807_int64*s,2147483647_int64)
      u=real(s,dp)/2147483647.0_dp
   end subroutine rng_u
end program
