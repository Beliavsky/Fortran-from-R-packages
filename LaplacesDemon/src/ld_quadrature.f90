module ld_quadrature
use ld_kinds, only: dp, pi
implicit none
private
public :: hermite_phys, gauss_hermite_rule, gauss_hermite_cube
contains
pure function hermite_phys(n,x) result(h)
   integer,intent(in)::n; real(dp),intent(in)::x; real(dp)::h,h0,h1,h2; integer::k
   if(n==0) then; h=1.0_dp; return; end if
   if(n==1) then; h=2.0_dp*x; return; end if
   h0=1.0_dp; h1=2.0_dp*x
   do k=2,n; h2=2.0_dp*x*h1-2.0_dp*real(k-1,dp)*h0; h0=h1; h1=h2; end do
   h=h1
end function hermite_phys

subroutine gauss_hermite_rule(n,nodes,weights)
   integer,intent(in)::n
   real(dp),intent(out)::nodes(n),weights(n)
   integer::i,j,m,it
   real(dp)::z,z1,p1,p2,p3,pp,eps,logfac
   eps=1e-14_dp; m=(n+1)/2
   do i=1,m
      if(i==1) then; z=sqrt(2.0_dp*real(n,dp)+1.0_dp)-1.85575_dp*(2.0_dp*real(n,dp)+1.0_dp)**(-1.0_dp/6.0_dp)
      else if(i==2) then; z=z-1.14_dp*real(n,dp)**0.426_dp/z
      else if(i==3) then; z=1.86_dp*z-0.86_dp*nodes(1)
      else if(i==4) then; z=1.91_dp*z-0.91_dp*nodes(2)
      else; z=2.0_dp*z-nodes(i-2); end if
      do it=1,50
         p1=pi**(-0.25_dp); p2=0.0_dp
         do j=1,n
            p3=p2; p2=p1; p1=z*sqrt(2.0_dp/real(j,dp))*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
         end do
         pp=sqrt(2.0_dp*real(n,dp))*p2; z1=z; z=z1-p1/pp; if(abs(z-z1)<eps) exit
      end do
      nodes(i)=z; nodes(n+1-i)=-z; weights(i)=2.0_dp/(pp*pp); weights(n+1-i)=weights(i)
   end do
   ! Above normalized-Hermite weights integrate exp(-x^2); normalize small roundoff to sqrt(pi).
   weights=weights*(sqrt(pi)/sum(weights))
end subroutine gauss_hermite_rule

subroutine gauss_hermite_cube(nodes1,weights1,dims,nodes,weights)
   real(dp),intent(in)::nodes1(:),weights1(:)
   integer,intent(in)::dims
   real(dp),allocatable,intent(out)::nodes(:,:),weights(:)
   integer::n1,total,i,d,idx,tmp
   n1=size(nodes1); total=n1**dims; allocate(nodes(total,dims),weights(total))
   do i=0,total-1
      tmp=i; weights(i+1)=1.0_dp
      do d=1,dims
         idx=mod(tmp,n1)+1; tmp=tmp/n1; nodes(i+1,d)=nodes1(idx); weights(i+1)=weights(i+1)*weights1(idx)
      end do
   end do
end subroutine gauss_hermite_cube
end module ld_quadrature
