module statmod_utils
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use r_compat, only: dp, median, pbinom, dbinom, pchisq, phyper, rbinom, fisher_test, fisher_test_result_t, runif1
use statmod_linalg, only: least_squares, mean_vec, variance_vec
use statmod_gaussquad, only: quad_rule_t, gauss_quad_prob
implicit none
private
public :: matvec, vecmat, mscale, rho_hampel, psi_hampel, forward_select
public :: hommel_test, permp, sage_test, power_fisher_test
public :: mean_t, compare_two_growth_curves
contains

pure function matvec(m,v) result(out)
real(dp),intent(in)::m(:,:),v(:)
real(dp)::out(size(m,1),size(m,2))
integer::j
out=m
do j=1,size(m,2)
out(:,j)=m(:,j)*v(j)
end do
end function

pure function vecmat(v,m) result(out)
real(dp),intent(in)::v(:),m(:,:)
real(dp)::out(size(m,1),size(m,2))
integer::i
out=m
do i=1,size(m,1)
out(i,:)=v(i)*m(i,:)
end do
end function

pure elemental function rho_hampel(u,a,b,c) result(rho)
real(dp),intent(in)::u
real(dp),intent(in),optional::a,b,c
real(dp)::rho,aa,bb,cc,z
aa=1.5_dp
if(present(a))aa=a
bb=3.5_dp
if(present(b))bb=b
cc=8.0_dp
if(present(c))cc=c
z=abs(u)
if(z<=aa)then
rho=z*z/2
else if(z<=bb)then
rho=aa*(z-aa/2)
else if(z<=cc)then
rho=aa*(bb-aa/2)+aa*(z-bb)*(1-(z-bb)/(cc-bb)/2)
else
rho=aa*(bb-aa+cc)/2
end if
end function

pure elemental function psi_hampel(u,a,b,c) result(psi)
real(dp),intent(in)::u
real(dp),intent(in),optional::a,b,c
real(dp)::psi,aa,bb,cc,z
aa=1.5_dp
if(present(a))aa=a
bb=3.5_dp
if(present(b))bb=b
cc=8.0_dp
if(present(c))cc=c
z=abs(u)
if(z<=aa)then
psi=u
else if(z<=bb)then
psi=sign(aa,u)
else if(z<=cc)then
psi=sign(aa*(cc-z)/(cc-bb),u)
else
psi=0
end if
end function

function mscale(u) result(s)
real(dp),intent(in)::u(:)
real(dp)::s,d1,d2
real(dp),allocatable::z(:),au(:)
integer::iter
if(size(u)==0)then
s=0
return
end if
if(real(count(u==0.0_dp),dp)/size(u)>=0.5_dp)then
s=0
return
end if
allocate(au(size(u)))
au=abs(u)
s=median(au)/0.6744898_dp
do iter=1,51
   z=u/(0.212_dp*s)
   d1=sum(rho_hampel(z))/size(z)-3.75_dp
   d2=sum(z*psi_hampel(z))/real(size(z),dp)
   if(d2==0)exit
   s=s*(1+d1/d2)
   if(abs(d1/d2)<1e-13_dp)exit
end do
end function

function forward_select(y,x,xkept,intercept,nvar) result(orderin)
real(dp),intent(in)::y(:),x(:,:)
real(dp),intent(in),optional::xkept(:,:)
logical,intent(in),optional::intercept
integer,intent(in),optional::nvar
integer,allocatable::orderin(:)
real(dp),allocatable::yy(:),xx(:,:),keep(:,:),beta(:),fit(:),res(:),bestx(:),corr(:),norms(:),tmp(:,:)
integer,allocatable::cand(:)
integer::n,nv,nin,j,k,best,rank,info,kc
logical::inc
n=size(y)
inc=.true.
if(present(intercept))inc=intercept
if(present(xkept))then
   if(inc)then
   allocate(keep(n,size(xkept,2)+1))
   keep(:,1)=1
   keep(:,2:)=xkept
   else
   keep=xkept
   end if
else if(inc)then
allocate(keep(n,1))
keep=1
else
allocate(keep(n,0))
end if
yy=y
xx=x
if(size(keep,2)>0)then
   call least_squares(keep,yy,beta,fit,res,rank,info)
   yy=res
   allocate(tmp(n,size(xx,2)))
   do j=1,size(xx,2)
   call least_squares(keep,xx(:,j),beta,fit,res,rank,info)
   tmp(:,j)=res
   end do
   call move_alloc(tmp,xx)
end if
nv=size(x,2)
if(present(nvar))nv=min(nv,nvar)
nv=min(nv,n-size(keep,2))
nv=max(0,nv)
allocate(orderin(nv),cand(size(x,2)))
cand=[(j,j=1,size(x,2))]
do nin=1,nv
   kc=size(cand)
   if(kc==1)then
   orderin(nin)=cand(1)
   exit
   end if
   yy=yy/sqrt(sum(yy*yy))
   allocate(norms(kc),corr(kc))
   do j=1,kc
   norms(j)=sqrt(sum(xx(:,j)**2))
   if(norms(j)>0)xx(:,j)=xx(:,j)/norms(j)
   corr(j)=dot_product(xx(:,j),yy)
   end do
   best=maxloc(abs(corr),dim=1)
   bestx=xx(:,best)
   orderin(nin)=cand(best)
   yy=yy-corr(best)*bestx
   if(kc>1)then
      allocate(tmp(n,kc-1))
      j=0
      do k=1,kc
         if(k==best)cycle
         j=j+1
         tmp(:,j)=xx(:,k)-dot_product(xx(:,k),bestx)*bestx
      end do
      call move_alloc(tmp,xx)
      cand=[cand(:best-1),cand(best+1:)]
   end if
   deallocate(norms,corr)
end do
end function

function hommel_test(p,alpha) result(reject)
real(dp),intent(in)::p(:)
real(dp),intent(in),optional::alpha
logical::reject(size(p))
real(dp)::a,po(size(p)),tmp
integer::n,j,k,i
n=size(p)
a=0.05_dp
if(present(alpha))a=alpha
po=p
! simple insertion sort
do i=2,n
   tmp=po(i)
   k=i-1
   do while(k>=1)
      if(po(k)<=tmp) exit
      po(k+1)=po(k)
      k=k-1
   end do
   po(k+1)=tmp
end do
j=n
do
   if(all([(po(n-j+k)>real(k,dp)*a/real(j,dp),k=1,j)]))exit
   j=j-1
   if(j==0)exit
end do
if(j==0)then
reject=.true.
else
reject=p>=a/real(j,dp)
end if
end function

function permp(x,nperm,n1,n2,total_nperm,method_exact,twosided) result(pout)
integer,intent(in)::x(:),nperm,n1,n2
integer(kind=8),intent(in),optional::total_nperm
logical,intent(in),optional::method_exact,twosided
real(dp)::pout(size(x))
integer(kind=8)::tot,k
logical::exact,two
integer::i,j
real(dp)::pr,intv
real(dp),allocatable::nodes(:),weights(:)
type(quad_rule_t)::q
two=.true.
if(present(twosided))two=twosided
if(present(total_nperm))then
tot=total_nperm
else
tot=nchoosek_int8(n1+n2,n1)
if(n1==n2.and.two)tot=tot/2
end if
if(present(method_exact))then
exact=method_exact
else
exact=tot<=10000
end if
if(exact)then
   do i=1,size(x)
      pout(i)=0
      do k=1,tot
      pr=real(k,dp)/real(tot,dp)
      pout(i)=pout(i)+pbinom(real(x(i),dp),nperm,pr)
      end do
      pout(i)=pout(i)/real(tot,dp)
   end do
else
   q=gauss_quad_prob(128,'uniform',0.0_dp,0.5_dp/real(tot,dp))
   nodes=q%nodes
   weights=q%weights
   do i=1,size(x)
      intv=0
      do j=1,128
      intv=intv+weights(j)*pbinom(real(x(i),dp),nperm,nodes(j))
      end do
      intv=0.5_dp/real(tot,dp)*intv
      pout(i)=real(x(i)+1,dp)/real(nperm+1,dp)-intv
   end do
end if
end function

pure function nchoosek_int8(n,k) result(v)
integer,intent(in)::n,k
integer(kind=8)::v
integer::i,kk
kk=min(k,n-k)
v=1
if(kk<0)then
v=0
return
end if
do i=1,kk
v=v*int(n-kk+i,8)/int(i,8)
end do
end function

function sage_test(x,y,n1,n2) result(pv)
integer,intent(in)::x(:),y(:)
integer,intent(in),optional::n1,n2
real(dp)::pv(size(x)),prob,pobs,pcur,chi,den
integer::nn1,nn2,i,k,sz,a,b,c,d
nn1=sum(x)
if(present(n1))nn1=n1
nn2=sum(y)
if(present(n2))nn2=n2
pv=1
if(nn1==nn2)then
   do i=1,size(x)
   sz=x(i)+y(i)
   if(sz>0)pv(i)=min(2.0_dp*pbinom(real(min(x(i),y(i)),dp),sz,0.5_dp),1.0_dp)
   end do
   return
end if
prob=real(nn1,dp)/real(nn1+nn2,dp)
do i=1,size(x)
   sz=x(i)+y(i)
   if(sz<=0)cycle
   if(sz>10000)then
      a=x(i)
      b=y(i)
      c=nn1-a
      d=nn2-b
      den=real((a+b)*(c+d)*(a+c)*(b+d),dp)
      if(den>0)then
      chi=real((a*d-b*c),dp)**2*real(a+b+c+d,dp)/den
      pv(i)=1.0_dp-pchisq(chi,1.0_dp)
      end if
   else
      pobs=dbinom(real(x(i),dp),sz,prob)
      pv(i)=0
      do k=0,sz
      pcur=dbinom(real(k,dp),sz,prob)
      if(pcur<=pobs*(1+1000*epsilon(1.0_dp)))pv(i)=pv(i)+pcur
      end do
      pv(i)=min(1.0_dp,pv(i))
   end if
end do
end function

function power_fisher_test(p1,p2,n1,n2,alpha,nsim,alternative) result(power)
real(dp),intent(in)::p1,p2
integer,intent(in)::n1,n2
real(dp),intent(in),optional::alpha
integer,intent(in),optional::nsim
character(len=*),intent(in),optional::alternative
real(dp)::power,a,pv
integer::ns,i,mwhite,mblack
integer,allocatable::y1(:),y2(:)
integer::tab(2,2)
character(len=16)::alt
type(fisher_test_result_t)::ft
a=0.05_dp
if(present(alpha))a=alpha
ns=100
if(present(nsim))ns=nsim
alt='two.sided'
if(present(alternative))alt=trim(alternative)
y1=rbinom(ns,n1,p1)
y2=rbinom(ns,n2,p2)
power=0.0_dp
do i=1,ns
   mwhite=y1(i)+y2(i)
   mblack=n1+n2-mwhite
   select case(trim(alt))
   case('less')
      pv=phyper(real(y1(i),dp),mwhite,mblack,n1)
   case('greater')
      if(y1(i)<=0)then
         pv=1.0_dp
      else
         pv=1.0_dp-phyper(real(y1(i)-1,dp),mwhite,mblack,n1)
      end if
   case default
      tab=reshape([y1(i),y2(i),n1-y1(i),n2-y2(i)],[2,2])
      ft=fisher_test(tab)
      pv=ft%p_value
   end select
   if(pv<a)power=power+1.0_dp
end do
power=power/real(ns,dp)
end function

function mean_t(y1,y2) result(stat)
real(dp),intent(in)::y1(:,:),y2(:,:)
real(dp)::stat,m1,m2,v1,v2,s,t,w,sw
integer::j,n1,n2,i
stat=0
sw=0
do j=1,size(y1,2)
   n1=count(ieee_is_finite(y1(:,j)))
   n2=count(ieee_is_finite(y2(:,j)))
   if(n1<2.or.n2<2)cycle
   m1=sum(y1(:,j),mask=ieee_is_finite(y1(:,j)))/n1
   m2=sum(y2(:,j),mask=ieee_is_finite(y2(:,j)))/n2
   v1=sum((y1(:,j)-m1)**2,mask=ieee_is_finite(y1(:,j)))/(n1-1)
   v2=sum((y2(:,j)-m2)**2,mask=ieee_is_finite(y2(:,j)))/(n2-1)
   s=((n1-1)*v1+(n2-1)*v2)/(n1+n2-2)
   if(s<=0)cycle
   t=(m1-m2)/sqrt(s*(1.0_dp/n1+1.0_dp/n2))
   w=real(n1+n2-2,dp)/real(n1+n2,dp)
   stat=stat+w*t
   sw=sw+w
end do
if(sw>0)stat=stat/sw
end function

subroutine compare_two_growth_curves(group,y,nsim,stat,p_value,n0)
integer,intent(in)::group(:)
real(dp),intent(in)::y(:,:)
integer,intent(in),optional::nsim
real(dp),intent(out)::stat,p_value
real(dp),intent(in),optional::n0
integer::ns,i,j,k,g1,g2,n
integer,allocatable::idx(:),pg(:)
real(dp)::s,base,u,asbig
ns=100
if(present(nsim))ns=nsim
base=0.5_dp
if(present(n0))base=n0
g1=group(1)
g2=g1
do i=1,size(group)
   if(group(i)/=g1)then
   g2=group(i)
   exit
   end if
end do
idx=pack([(i,i=1,size(group))],group==g1)
stat=mean_t(y(idx,:),y(pack([(i,i=1,size(group))],group==g2),:))
asbig=0.0_dp
n=size(group)
allocate(pg(n))
pg=group
do i=1,ns
   pg=group
   do j=n,2,-1
      u=runif1()
      k=1+int(u*real(j,dp))
      if(k>j)k=j
      call swap_int(pg(j),pg(k))
   end do
   idx=pack([(j,j=1,n)],pg==g1)
   s=mean_t(y(idx,:),y(pack([(j,j=1,n)],pg==g2),:))
   if(abs(s)>abs(stat))asbig=asbig+1.0_dp
   if(abs(s)==abs(stat))asbig=asbig+0.5_dp
end do
p_value=(asbig+base)/(real(ns,dp)+base)
end subroutine

subroutine swap_int(a,b)
integer,intent(inout)::a,b
integer::t
t=a
a=b
b=t
end subroutine

end module statmod_utils
