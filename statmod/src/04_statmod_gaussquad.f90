! Gaussian quadrature derived from statmod R/gaussquad.R and src/gaussq2.f.
! gausq2 lineage: Netlib gaussq.f / EISPACK IMTQL2, Martin-Wilkinson/Dubrulle.
module statmod_gaussquad
use r_compat, only: dp
implicit none
private
public :: quad_rule_t, gauss_quad, gauss_quad_prob
real(dp),parameter::pi=acos(-1.0_dp)
type :: quad_rule_t
   real(dp),allocatable :: nodes(:),weights(:)
   integer :: ierr=0
end type
contains

subroutine gausq2(n,d,e,z,ierr)
integer,intent(in)::n
real(dp),intent(inout)::d(n),e(n),z(n)
integer,intent(out)::ierr
integer::i,j,k,l,m,ii,mml
real(dp)::b,c,f,g,p,r,s,machep,tmp
machep=epsilon(1.0_dp)
ierr=0
if(n==1) return
e(n)=0.0_dp
do l=1,n
   j=0
   do
      m=l
      do while(m<n)
         if(abs(e(m))<=machep*(abs(d(m))+abs(d(m+1)))) exit
         m=m+1
      end do
      p=d(l)
      if(m==l) exit
      if(j==30) then
      ierr=l
      return
      end if
      j=j+1
      g=(d(l+1)-p)/(2.0_dp*e(l))
      r=sqrt(g*g+1.0_dp)
      g=d(m)-p+e(l)/(g+sign(r,g))
      s=1.0_dp
      c=1.0_dp
      p=0.0_dp
      mml=m-l
      do ii=1,mml
         i=m-ii
         f=s*e(i)
         b=c*e(i)
         if(abs(f)>=abs(g)) then
            c=g/f
            r=sqrt(c*c+1.0_dp)
            e(i+1)=f*r
            s=1.0_dp/r
            c=c*s
         else
            s=f/g
            r=sqrt(s*s+1.0_dp)
            e(i+1)=g*r
            c=1.0_dp/r
            s=s*c
         end if
         g=d(i+1)-p
         r=(d(i)-g)*s+2.0_dp*c*b
         p=s*r
         d(i+1)=g+p
         g=c*r-b
         f=z(i+1)
         z(i+1)=s*z(i)+c*f
         z(i)=c*z(i)-s*f
      end do
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.0_dp
   end do
end do
! Order eigenvalues and first eigenvector components.
do ii=2,n
   i=ii-1
   k=i
   p=d(i)
   do j=ii,n
      if(d(j)<p) then
      k=j
      p=d(j)
      end if
   end do
   if(k/=i) then
      d(k)=d(i)
      d(i)=p
      tmp=z(i)
      z(i)=z(k)
      z(k)=tmp
   end if
end do
end subroutine gausq2

function gauss_quad(n,kind,alpha,beta) result(out)
integer,intent(in)::n
character(len=*),intent(in),optional::kind
real(dp),intent(in),optional::alpha,beta
type(quad_rule_t)::out
character(len=16)::knd
real(dp)::aa,bb,ab,abi,lnmuzero
real(dp),allocatable::a(:),b(:),z(:)
integer::i
if(n<0) then
out%ierr=-1
allocate(out%nodes(0),out%weights(0))
return
end if
allocate(out%nodes(n),out%weights(n))
if(n==0) return
knd='legendre'
if(present(kind)) knd=trim(kind)
aa=0.0_dp
if(present(alpha)) aa=alpha
bb=0.0_dp
if(present(beta)) bb=beta
allocate(a(n),b(n),z(n))
a=0
b=0
z=0
z(1)=1
select case(trim(knd))
case('legendre')
   lnmuzero=log(2.0_dp)
   do i=1,n-1
   b(i)=real(i,dp)/sqrt(4.0_dp*real(i*i,dp)-1.0_dp)
   end do
case('chebyshev1')
   lnmuzero=log(pi)
   if(n>1) b(1)=sqrt(0.5_dp)
   if(n>2) b(2:n-1)=0.5_dp
case('chebyshev2')
   lnmuzero=log(pi/2.0_dp)
   if(n>1) b(1:n-1)=0.5_dp
case('hermite')
   lnmuzero=log(pi)/2.0_dp
   do i=1,n-1
   b(i)=sqrt(real(i,dp)/2.0_dp)
   end do
case('jacobi')
   ab=aa+bb
   lnmuzero=(ab+1.0_dp)*log(2.0_dp)+log_gamma(aa+1.0_dp)+log_gamma(bb+1.0_dp)-log_gamma(ab+2.0_dp)
   a(1)=(bb-aa)/(ab+2.0_dp)
   do i=2,n
      abi=ab+2.0_dp*real(i,dp)
      a(i)=(bb*bb-aa*aa)/(abi-2.0_dp)/abi
   end do
   if(n>1) b(1)=sqrt(4.0_dp*(aa+1)*(bb+1)/(ab+2)**2/(ab+3))
   do i=2,n-1
      abi=ab+2.0_dp*real(i,dp)
      b(i)=sqrt(4.0_dp*i*(i+aa)*(i+bb)*(i+ab)/(abi*abi-1.0_dp)/(abi*abi))
   end do
case('laguerre')
   do i=1,n
   a(i)=2.0_dp*i-1.0_dp+aa
   end do
   do i=1,n-1
   b(i)=sqrt(real(i,dp)*(real(i,dp)+aa))
   end do
   lnmuzero=log_gamma(aa+1.0_dp)
case default
   out%ierr=-2
   out%nodes=0
   out%weights=0
   return
end select
call gausq2(n,a,b,z,out%ierr)
out%nodes=a
out%weights=exp(lnmuzero+2.0_dp*log(abs(z)))
end function gauss_quad

function gauss_quad_prob(n,dist,l,u,mu,sigma,alpha,beta) result(out)
integer,intent(in)::n
character(len=*),intent(in),optional::dist
real(dp),intent(in),optional::l,u,mu,sigma,alpha,beta
type(quad_rule_t)::out
character(len=16)::dst
real(dp)::ll,uu,mm,ss,aa,bb,ab,abi
real(dp),allocatable::a(:),b(:),z(:)
integer::i
if(n<0) then
out%ierr=-1
allocate(out%nodes(0),out%weights(0))
return
end if
allocate(out%nodes(n),out%weights(n))
if(n==0) return
dst='uniform'
if(present(dist)) dst=trim(dist)
ll=0
if(present(l)) ll=l
uu=1
if(present(u)) uu=u
mm=0
if(present(mu)) mm=mu
ss=1
if(present(sigma)) ss=sigma
aa=1
if(present(alpha)) aa=alpha
bb=1
if(present(beta)) bb=beta
if(n==1) then
   select case(trim(dst))
   case('uniform'); out%nodes(1)=(ll+uu)/2
   case('beta1','beta2','beta'); out%nodes(1)=aa/(aa+bb)
   case('normal'); out%nodes(1)=mm
   case('gamma'); out%nodes(1)=aa*bb
   case default
   out%ierr=-2
   out%nodes=0
   out%weights=0
   return
   end select
   out%weights=1
   return
end if
if(trim(dst)=='beta' .and. aa==0.5_dp .and. bb==0.5_dp) dst='beta1'
if(trim(dst)=='beta' .and. aa==1.5_dp .and. bb==1.5_dp) dst='beta2'
allocate(a(n),b(n),z(n))
a=0
b=0
z=0
z(1)=1
select case(trim(dst))
case('uniform')
   do i=1,n-1
   b(i)=real(i,dp)/sqrt(4.0_dp*real(i*i,dp)-1.0_dp)
   end do
case('beta1')
   b(1)=sqrt(0.5_dp)
   if(n>2)b(2:n-1)=0.5_dp
case('beta2')
   b(1:n-1)=0.5_dp
case('normal')
   do i=1,n-1
   b(i)=sqrt(real(i,dp)/2.0_dp)
   end do
case('beta')
   ab=aa+bb
   a(1)=(aa-bb)/ab
   do i=2,n
      abi=ab-2.0_dp+2.0_dp*i
      a(i)=((aa-1)**2-(bb-1)**2)/(abi-2.0_dp)/abi
   end do
   b(1)=sqrt(4.0_dp*aa*bb/ab**2/(ab+1.0_dp))
   do i=2,n-1
      abi=ab-2.0_dp+2.0_dp*i
      b(i)=sqrt(4.0_dp*i*(i+aa-1)*(i+bb-1)*(i+ab-2)/(abi**2-1)/abi**2)
   end do
case('gamma')
   do i=1,n
   a(i)=2.0_dp*i+aa-2.0_dp
   end do
   do i=1,n-1
   b(i)=sqrt(real(i,dp)*(real(i,dp)+aa-1.0_dp))
   end do
case default
   out%ierr=-2
   out%nodes=0
   out%weights=0
   return
end select
call gausq2(n,a,b,z,out%ierr)
out%nodes=a
out%weights=z*z
select case(trim(dst))
case('uniform'); out%nodes=ll+(uu-ll)*(out%nodes+1.0_dp)/2.0_dp
case('beta1','beta2','beta'); out%nodes=(out%nodes+1.0_dp)/2.0_dp
case('normal'); out%nodes=mm+sqrt(2.0_dp)*ss*out%nodes
case('gamma'); out%nodes=bb*out%nodes
end select
end function gauss_quad_prob

end module statmod_gaussquad
