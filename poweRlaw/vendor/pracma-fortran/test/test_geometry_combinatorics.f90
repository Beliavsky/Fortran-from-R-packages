! SPDX-License-Identifier: GPL-3.0-or-later
program test_geometry_combinatorics
   use pracma
   implicit none
   type(circle_result)::cr
   type(qpspecial_result)::qr
   real(dp),allocatable::h(:,:),q(:,:),tx(:)
   integer,allocatable::ps(:),cb(:,:),bv(:)
   logical,allocatable::inside(:)
   integer(i8),allocatable::fac(:)
   real(dp)::pint(2)
   logical::hit
   cr=circlefit([1.0_dp,0.0_dp,-1.0_dp,0.0_dp],[0.0_dp,1.0_dp,0.0_dp,-1.0_dp])
   call check(cr%status==pracma_ok.and.abs(cr%radius-1.0_dp)<1e-12_dp,'circlefit')
   inside=inpolygon([0.5_dp,2.0_dp],[0.5_dp,0.5_dp],[0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
                    [0.0_dp,0.0_dp,1.0_dp,1.0_dp])
   call check(all(inside.eqv.[.true.,.false.]),'inpolygon')
   pint=segment_intersection([0.0_dp,0.0_dp],[1.0_dp,1.0_dp],[0.0_dp,1.0_dp],[1.0_dp,0.0_dp],hit)
   call check(hit.and.maxval(abs(pint-[0.5_dp,0.5_dp]))<1e-14_dp,'segment intersection')
   call check(isprime(97_i8).and..not.isprime(91_i8),'isprime')
   ps=primes(20); call check(all(ps==[2,3,5,7,11,13,17,19]),'primes')
   fac=factors(84_i8); call check(all(fac==[2_i8,2_i8,3_i8,7_i8]),'factors')
   bv=bits(5_i8); call check(all(bv==[1,0,1]),'bits compatibility')
   tx=trisolve([2.0_dp,2.0_dp],[-1.0_dp],[-1.0_dp],[1.0_dp,0.0_dp])
   call check(maxval(abs(tx-[2.0_dp/3.0_dp,1.0_dp/3.0_dp]))<1e-13_dp,'trisolve compatibility')
   cb=combs([1,2,3,4],2); call check(size(cb,1)==6.and.size(cb,2)==2,'combs')
   h=hadamard(4); call check(maxval(abs(matmul(h,transpose(h))-4.0_dp*eye(4)))<1e-14_dp,'hadamard')
   q=randortho(4,123_i8); call check(maxval(abs(matmul(transpose(q),q)-eye(4)))<1e-10_dp,'randortho')
   qr=qpspecial(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))
   call check(qr%converged.and.abs(sum(qr%x)-1.0_dp)<1e-10_dp.and.minval(qr%x)>=-1e-12_dp,'qpspecial')
   print '(a)','test_geometry_combinatorics: PASS'
contains
   subroutine check(ok,name)
      logical,intent(in)::ok; character(len=*),intent(in)::name
      if(.not.ok)then; write(*,'(a,1x,a)')'FAIL:',trim(name); error stop 1; end if
   end subroutine check
end program test_geometry_combinatorics
