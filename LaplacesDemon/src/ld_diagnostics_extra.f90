module ld_diagnostics_extra
use ld_kinds, only: dp
use ld_distributions, only: normal_quantile
implicit none
private
public :: raftery_result_t, raftery_diagnostic, hangartner_chisq

type :: raftery_result_t
   integer, allocatable :: burn(:), total(:), minimum(:), kthinning(:)
   real(dp), allocatable :: dependence_factor(:)
   logical :: enough_samples=.true.
end type raftery_result_t

contains

subroutine raftery_diagnostic(chain,res,q,r,s,eps,original_thinning)
   real(dp), intent(in) :: chain(:,:)
   type(raftery_result_t), intent(out) :: res
   real(dp), intent(in), optional :: q,r,s,eps
   integer, intent(in), optional :: original_thinning
   real(dp) :: qq,rr,ss,ee,phi,quant,alpha,beta,tempburn,tempprec,iratio
   integer :: n,p,j,nmin,kthin,thin0,newdim,nburn,nkeep,i
   logical, allocatable :: dichot(:),testres(:)
   real(dp) :: bic,g2
   integer :: trans3(2,2,2),trans2(2,2)
   qq=0.025_dp; if(present(q)) qq=q
   rr=0.005_dp; if(present(r)) rr=r
   ss=0.95_dp; if(present(s)) ss=s
   ee=0.001_dp; if(present(eps)) ee=eps
   thin0=1; if(present(original_thinning)) thin0=max(1,original_thinning)
   n=size(chain,1); p=size(chain,2); phi=normal_quantile(0.5_dp*(1.0_dp+ss),0.0_dp,1.0_dp)
   nmin=ceiling((qq*(1.0_dp-qq)*phi*phi)/(rr*rr))
   allocate(res%burn(p),res%total(p),res%minimum(p),res%kthinning(p),res%dependence_factor(p))
   res%burn=0; res%total=0; res%minimum=nmin; res%kthinning=thin0; res%dependence_factor=0.0_dp
   if(nmin>n) then; res%enough_samples=.false.; return; end if
   allocate(dichot(n))
   do j=1,p
      quant=sample_quantile(chain(:,j),qq); dichot=(chain(:,j)<=quant); kthin=0; bic=1.0_dp
      do while(bic>=0.0_dp)
         kthin=kthin+thin0
         call thin_logical(dichot,kthin,testres); newdim=size(testres)
         if(newdim<6) then; bic=-1.0_dp; exit; end if
         trans3=0
         do i=1,newdim-2
            trans3(li(testres(i)),li(testres(i+1)),li(testres(i+2)))=&
                 trans3(li(testres(i)),li(testres(i+1)),li(testres(i+2)))+1
         end do
         g2=markov_g2(trans3); bic=g2-log(real(newdim-2,dp))*2.0_dp
         if(kthin>max(1,n/4)) then; bic=-1.0_dp; exit; end if
      end do
      trans2=0
      do i=1,newdim-1
         trans2(li(testres(i)),li(testres(i+1)))=trans2(li(testres(i)),li(testres(i+1)))+1
      end do
      alpha=real(trans2(1,2),dp)/real(max(1,trans2(1,1)+trans2(1,2)),dp)
      beta =real(trans2(2,1),dp)/real(max(1,trans2(2,1)+trans2(2,2)),dp)
      if(alpha<=0.0_dp .or. beta<=0.0_dp .or. abs(1.0_dp-alpha-beta)<=1.0e-12_dp) then
         nburn=0; nkeep=nmin
      else
         tempburn=log(ee*(alpha+beta)/max(alpha,beta))/log(abs(1.0_dp-alpha-beta))
         tempprec=((2.0_dp-alpha-beta)*alpha*beta*phi*phi)/(((alpha+beta)**3)*rr*rr)
         nburn=max(0,ceiling(tempburn)*kthin); nkeep=max(1,ceiling(tempprec)*kthin)
      end if
      iratio=real(nburn+nkeep,dp)/real(max(1,nmin),dp)
      res%burn(j)=nburn; res%total(j)=nburn+nkeep; res%kthinning(j)=kthin; res%dependence_factor(j)=iratio
      if(allocated(testres)) deallocate(testres)
   end do
end subroutine raftery_diagnostic

subroutine hangartner_chisq(x,groups,statistic,df)
   integer, intent(in) :: x(:),groups
   real(dp), intent(out) :: statistic
   integer, intent(out) :: df
   integer :: n,j,g,k,minx,maxx,ncat,idx
   integer, allocatable :: tab(:,:),row(:),col(:)
   real(dp) :: expected
   n=size(x); if(groups<2 .or. mod(n,groups)/=0) then; statistic=huge(1.0_dp); df=0; return; end if
   minx=minval(x); maxx=maxval(x); ncat=maxx-minx+1
   allocate(tab(ncat,groups),row(ncat),col(groups)); tab=0
   do j=1,n
      g=1+(j-1)/(n/groups); idx=x(j)-minx+1; tab(idx,g)=tab(idx,g)+1
   end do
   row=sum(tab,dim=2); col=sum(tab,dim=1); statistic=0.0_dp
   do k=1,ncat
      if(row(k)==0) cycle
      do g=1,groups
         expected=real(row(k)*col(g),dp)/real(n,dp)
         if(expected>0.0_dp) statistic=statistic+(real(tab(k,g),dp)-expected)**2/expected
      end do
   end do
   df=max(0,(count(row>0)-1)*(groups-1))
end subroutine hangartner_chisq

pure integer function li(x)
   logical, intent(in) :: x
   if(x) then; li=2; else; li=1; end if
end function li

subroutine thin_logical(x,by,y)
   logical, intent(in) :: x(:)
   integer, intent(in) :: by
   logical, allocatable, intent(out) :: y(:)
   integer :: n,i,k
   n=1+(size(x)-1)/max(1,by); allocate(y(n)); k=0
   do i=1,size(x),max(1,by); k=k+1; y(k)=x(i); end do
end subroutine thin_logical

pure function markov_g2(tab) result(g2)
   integer, intent(in) :: tab(2,2,2)
   real(dp) :: g2,fitted,den
   integer :: i1,i2,i3
   g2=0.0_dp
   do i1=1,2; do i2=1,2; do i3=1,2
      if(tab(i1,i2,i3)>0) then
         den=real(sum(tab(:,i2,:)),dp)
         if(den>0.0_dp) then
            fitted=real(sum(tab(i1,i2,:))*sum(tab(:,i2,i3)),dp)/den
            if(fitted>0.0_dp) g2=g2+2.0_dp*real(tab(i1,i2,i3),dp)*log(real(tab(i1,i2,i3),dp)/fitted)
         end if
      end if
   end do; end do; end do
end function markov_g2

function sample_quantile(x,p) result(q)
   real(dp), intent(in) :: x(:),p
   real(dp) :: q
   real(dp), allocatable :: y(:)
   real(dp) :: h,w
   integer :: n,k
   y=x; call insertion_sort(y); n=size(y)
   if(p<=0.0_dp) then; q=y(1); return; end if
   if(p>=1.0_dp) then; q=y(n); return; end if
   h=1.0_dp+(real(n-1,dp))*p; k=floor(h); w=h-real(k,dp)
   if(k>=n) then; q=y(n); else; q=(1.0_dp-w)*y(k)+w*y(k+1); end if
end function sample_quantile

subroutine insertion_sort(x)
   real(dp), intent(inout) :: x(:)
   integer :: i,j
   real(dp) :: key
   do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
         if(x(j)<=key) exit
         x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
   end do
end subroutine insertion_sort

end module ld_diagnostics_extra
