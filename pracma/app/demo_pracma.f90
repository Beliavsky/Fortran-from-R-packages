! SPDX-License-Identifier: GPL-3.0-or-later
module demo_pracma_functions
   use pracma_kinds, only : dp
   implicit none
contains
   function froot(x) result(y)
      real(dp),intent(in)::x
      real(dp)::y
      y=cos(x)-x
   end function froot
end module demo_pracma_functions

program demo_pracma
   use pracma
   use demo_pracma_functions
   implicit none
   type(root_result) :: root
   type(circle_result) :: circle
   real(dp), allocatable :: xi(:), yi(:), h(:,:)

   root=brentDekker(froot,0.0_dp,1.0_dp)
   xi=linspace(0.0_dp,2.0_dp*pi_dp,9)
   yi=interp1_pchip(xi,sin(xi),[pi_dp/6.0_dp,pi_dp/3.0_dp])
   circle=circlefit(cos(xi(:8)),sin(xi(:8)))
   h=hadamard(4)

   print '(a,f14.10)','root of cos(x)-x: ',root%root
   print '(a,2f12.6)','PCHIP sine values: ',yi
   print '(a,f12.6)','fitted unit-circle radius: ',circle%radius
   print '(a,f12.6)','Hadamard orthogonality error: ',maxval(abs(matmul(h,transpose(h))-4.0_dp*eye(4)))
end program demo_pracma
