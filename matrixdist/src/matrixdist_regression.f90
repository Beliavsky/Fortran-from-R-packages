! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_regression
   use r_compat, only: dp
   use matrixdist_iph, only: iph_density, iph_cdf
   implicit none
   private
   public :: sph_density, sph_survival, sph_loglik
contains
   function sph_density(y,x,b,alpha,s,kind,beta,model_type,c) result(f)
      real(dp),intent(in)::y,x(:),b(:),alpha(:),s(:,:),beta(:)
      character(len=*),intent(in)::kind,model_type
      real(dp),intent(in),optional::c(:)
      real(dp)::f,eta,sc
      real(dp),allocatable::ss(:,:),bb(:)
      if(size(x)/=size(b))error stop 'sph_density: x/B mismatch'
      eta=dot_product(x,b)
      sc=exp(eta)
      bb=beta
      select case(trim(model_type))
      case('reg')
         ss=sc*s
         f=iph_density(y,alpha,ss,kind,bb)
      case('aft')
         f=iph_density(y,alpha,s,kind,bb,scale=sc)
      case('reg2')
         if(.not.present(c))error stop 'sph_density: reg2 requires C coefficients'
         if(size(c)/=size(x))error stop 'sph_density: x/C mismatch'
         bb=beta*exp(dot_product(x,c))
         ss=sc*s
         f=iph_density(y,alpha,ss,kind,bb)
      case default
         error stop 'sph_density: model_type must be reg, reg2, or aft'
      end select
   end function sph_density

   function sph_survival(y,x,b,alpha,s,kind,beta,model_type,c) result(q)
      real(dp),intent(in)::y,x(:),b(:),alpha(:),s(:,:),beta(:)
      character(len=*),intent(in)::kind,model_type
      real(dp),intent(in),optional::c(:)
      real(dp)::q,eta,sc
      real(dp),allocatable::ss(:,:),bb(:)
      if(size(x)/=size(b))error stop 'sph_survival: x/B mismatch'
      eta=dot_product(x,b)
      sc=exp(eta)
      bb=beta
      select case(trim(model_type))
      case('reg')
         ss=sc*s
         q=iph_cdf(y,alpha,ss,kind,bb,.false.)
      case('aft')
         q=iph_cdf(y,alpha,s,kind,bb,.false.,scale=sc)
      case('reg2')
         if(.not.present(c))error stop 'sph_survival: reg2 requires C coefficients'
         if(size(c)/=size(x))error stop 'sph_survival: x/C mismatch'
         bb=beta*exp(dot_product(x,c))
         ss=sc*s
         q=iph_cdf(y,alpha,ss,kind,bb,.false.)
      case default
         error stop 'sph_survival: model_type must be reg, reg2, or aft'
      end select
   end function sph_survival

   function sph_loglik(y,x,delta,b,alpha,s,kind,beta,model_type,weight,c) result(ll)
      real(dp),intent(in)::y(:),x(:,:),b(:),alpha(:),s(:,:),beta(:)
      logical,intent(in)::delta(:)
      character(len=*),intent(in)::kind,model_type
      real(dp),intent(in),optional::weight(:),c(:)
      real(dp)::ll,f,w
      integer::i
      if(size(x,1)/=size(y) .or. size(delta)/=size(y))error stop 'sph_loglik: dimension mismatch'
      ll=0.0_dp
      do i=1,size(y)
         w=1.0_dp
         if(present(weight))w=weight(i)
         if(delta(i))then
            if(present(c))then
               f=sph_density(y(i),x(i,:),b,alpha,s,kind,beta,model_type,c)
            else
               f=sph_density(y(i),x(i,:),b,alpha,s,kind,beta,model_type)
            end if
         else
            if(present(c))then
               f=sph_survival(y(i),x(i,:),b,alpha,s,kind,beta,model_type,c)
            else
               f=sph_survival(y(i),x(i,:),b,alpha,s,kind,beta,model_type)
            end if
         end if
         if(f<=0.0_dp)then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function sph_loglik
end module matrixdist_regression
