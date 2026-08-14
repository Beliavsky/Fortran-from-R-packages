! SPDX-License-Identifier: GPL-2.0-only
module mco_test_functions
   use mco_kinds, only : dp, pi
   implicit none
   private
   public :: belegundu, belegundu_constraints, binh1, binh2, binh2_constraints, binh3
   public :: deb3, fonseca1, fonseca2, gianna
   public :: hanne1, hanne1_constraints, hanne2, hanne2_constraints
   public :: hanne3, hanne3_constraints, hanne4, hanne4_constraints
   public :: hanne5, hanne5_constraints, jimenez, jimenez_constraints
   public :: vnt, zdt1, zdt2, zdt3
contains
   subroutine belegundu(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[-2.0_dp*x(1)+x(2),2.0_dp*x(1)+x(2)]
   end subroutine
   subroutine belegundu_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
      g=[x(1)-x(2)+1.0_dp,-x(1)-x(2)+7.0_dp]
   end subroutine
   subroutine binh1(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[sum(x*x),sum((x-5.0_dp)**2)]
   end subroutine
   subroutine binh2(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[sum((4.0_dp*x)**2),sum((x-5.0_dp)**2)]
   end subroutine
   subroutine binh2_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
      g(1)=-((x(1)-5.0_dp)**2+x(2)**2-25.0_dp)
      g(2)=-(-(x(1)-8.0_dp)**2-(x(2)+3.0_dp)**2+7.7_dp)
   end subroutine
   subroutine binh3(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[x(1)-1.0e7_dp,x(2)-2.0e-6_dp,x(1)*x(2)-2.0_dp]
   end subroutine
   subroutine deb3(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      real(dp)::ff,gg,hh
      ff=4.0_dp*x(1)
      if(x(2)<=0.4_dp) then
         gg=4.0_dp-3.0_dp*exp(-((x(2)-0.2_dp)/0.02_dp)**2)
      else
         gg=4.0_dp-2.0_dp*exp(-((x(2)-0.7_dp)/0.2_dp)**2)
      end if
      if(ff<=gg) then
         hh=gg*(1.0_dp-(ff/gg)**(0.25_dp+3.75_dp*(gg-1.0_dp)))
      else
         hh=0.0_dp
      end if
      f=[ff,hh]
   end subroutine
   subroutine fonseca1(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[1.0_dp-exp(-(x(1)-1.0_dp)**2-(x(2)+1.0_dp)**2), &
         1.0_dp-exp(-(x(1)+1.0_dp)**2-(x(2)-1.0_dp)**2)]
   end subroutine
   subroutine fonseca2(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      real(dp)::c
      c=1.0_dp/sqrt(real(size(x),dp))
      f=[1.0_dp-exp(-sum((x-c)**2)),1.0_dp-exp(-sum((x+c)**2))]
   end subroutine
   subroutine gianna(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[1.0_dp/(sqrt(10.0_dp-x(1))+sqrt(x(1)-5.0_dp)),0.04_dp*(x(1)-8.0_dp)**2+0.3_dp]
   end subroutine
   subroutine hanne1(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); f=x
   end subroutine
   subroutine hanne1_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:); g(1)=sum(x)-5.0_dp
   end subroutine
   subroutine hanne2(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); f=x*x
   end subroutine
   subroutine hanne2_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:); g(1)=sum(x)-5.0_dp
   end subroutine
   subroutine hanne3(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); f=sqrt(x)
   end subroutine
   subroutine hanne3_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:); g(1)=sum(x)-5.0_dp
   end subroutine
   subroutine hanne4(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); f=x
   end subroutine
   subroutine hanne4_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
      g(1)=x(2)-5.0_dp+0.5_dp*x(1)*sin(4.0_dp*x(1))
   end subroutine
   subroutine hanne5(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      real(dp)::a,b
      a=2.0_dp*pi*(x(2)-aint(x(2))); b=x(1)-aint(x(1))
      f=[aint(x(1))+0.5_dp+b*sin(a),aint(x(2))+0.5_dp+b*cos(a)]
   end subroutine
   subroutine hanne5_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:); g(1)=sum(x)-5.0_dp
   end subroutine
   subroutine jimenez(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=-[5.0_dp*x(1)+3.0_dp*x(2),2.0_dp*x(1)+8.0_dp*x(2)]
   end subroutine
   subroutine jimenez_constraints(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
      g=-[x(1)+4.0_dp*x(2)-100.0_dp,3.0_dp*x(1)+2.0_dp*x(2)-150.0_dp, &
          200.0_dp-5.0_dp*x(1)-3.0_dp*x(2),75.0_dp-2.0_dp*x(1)-8.0_dp*x(2)]
   end subroutine
   subroutine vnt(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      real(dp)::xn
      xn=sum(x*x)
      f(1)=xn/2.0_dp+sin(xn)
      f(2)=(3.0_dp*x(1)-2.0_dp*x(2)+4.0_dp)**2/8.0_dp+(x(1)-x(2)+1.0_dp)**2/27.0_dp+15.0_dp
      f(3)=1.0_dp/(xn+1.0_dp)-1.1_dp*exp(-xn)
   end subroutine
   subroutine zdt1(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); real(dp)::g
      g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp); f=[x(1),g*(1.0_dp-sqrt(x(1)/g))]
   end subroutine
   subroutine zdt2(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); real(dp)::g
      g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp); f=[x(1),g*(1.0_dp-(x(1)/g)**2)]
   end subroutine
   subroutine zdt3(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); real(dp)::g
      g=1.0_dp+9.0_dp*sum(x(2:))/real(size(x)-1,dp)
      f=[x(1),g*(1.0_dp-sqrt(x(1)/g)-(x(1)/g)*sin(10.0_dp*pi*x(1)))]
   end subroutine
end module mco_test_functions
