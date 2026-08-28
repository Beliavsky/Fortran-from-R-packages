module lavaan_models
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model
   implicit none
   private
   public :: ram_from_lisrel
contains
   function ram_from_lisrel(lambda,beta,psi,theta,nu,alpha) result(model)
      real(dp),intent(in)::lambda(:,:),beta(:,:),psi(:,:),theta(:,:)
      real(dp),intent(in),optional::nu(:),alpha(:)
      type(ram_model)::model
      integer::p,q,i
      p=size(lambda,1)
      q=size(lambda,2)
      allocate(model%a(p+q,p+q),model%s(p+q,p+q),model%observed(p))
      model%a=0
      model%s=0
      if(present(nu) .or. present(alpha)) then
         allocate(model%m(p+q))
         model%m=0
      end if
      model%a(1:p,p+1:p+q)=lambda
      model%a(p+1:p+q,p+1:p+q)=beta
      model%s(1:p,1:p)=theta
      model%s(p+1:p+q,p+1:p+q)=psi
      if(present(nu)) model%m(1:p)=nu
      if(present(alpha)) model%m(p+1:p+q)=alpha
      model%observed=[(i,i=1,p)]
   end function ram_from_lisrel
end module lavaan_models
