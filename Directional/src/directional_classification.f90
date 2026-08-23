module directional_classification
   use directional_kinds, only : dp
   use directional_inference, only : vmf_mle, vmf_mle_result
   use directional_special, only : log_bessel_i
   implicit none
   private
   public :: dirknn, dirda_vmf
contains
   subroutine dirknn(xnew,x,groups,k,pred)
      real(dp),intent(in)::xnew(:,:),x(:,:);integer,intent(in)::groups(:),k;integer,intent(out)::pred(size(xnew,1))
      real(dp)::dist(size(x,1));integer::idx(size(x,1)),counts(maxval(groups)),i,j,kk,t,tmp
      do i=1,size(xnew,1)
         do j=1,size(x,1);dist(j)=1.0_dp-dot_product(xnew(i,:),x(j,:));idx(j)=j;end do
         do j=1,min(k,size(x,1));t=j;do kk=j+1,size(x,1);if(dist(idx(kk))<dist(idx(t)))t=kk;end do;tmp=idx(j);idx(j)=idx(t);idx(t)=tmp;end do
         counts=0;do j=1,min(k,size(x,1));counts(groups(idx(j)))=counts(groups(idx(j)))+1;end do
         pred(i)=maxloc(counts,dim=1)
      end do
   end subroutine

   subroutine dirda_vmf(xnew,x,groups,pred,scores)
      real(dp),intent(in)::xnew(:,:),x(:,:);integer,intent(in)::groups(:);integer,intent(out)::pred(size(xnew,1));real(dp),intent(out),optional::scores(size(xnew,1),maxval(groups))
      real(dp)::sc(size(xnew,1),maxval(groups)),prior,nu;integer::g,j,nj;type(vmf_mle_result)::fit
      nu=0.5_dp*size(x,2)-1
      do g=1,maxval(groups)
         nj=count(groups==g);fit=vmf_mle(pack_rows(x,groups==g));prior=real(nj,dp)/size(x,1)
         do j=1,size(xnew,1);sc(j,g)=log(prior)+nu*log(max(fit%kappa,tiny(1.0_dp)))-0.5_dp*size(x,2)*log(2*acos(-1.0_dp))-log_bessel_i(nu,max(fit%kappa,tiny(1.0_dp)))+fit%kappa*dot_product(xnew(j,:),fit%mu);end do
      end do
      do j=1,size(xnew,1);pred(j)=maxloc(sc(j,:),dim=1);end do;if(present(scores))scores=sc
   end subroutine
   function pack_rows(x,mask) result(y)
      real(dp),intent(in)::x(:,:);logical,intent(in)::mask(:);real(dp),allocatable::y(:,:);integer::i,j
      allocate(y(count(mask),size(x,2)));j=0;do i=1,size(x,1);if(mask(i))then;j=j+1;y(j,:)=x(i,:);end if;end do
   end function
end module directional_classification
