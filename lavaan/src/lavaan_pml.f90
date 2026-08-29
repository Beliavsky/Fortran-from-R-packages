module lavaan_pml
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma
   use lavaan_fit, only : sem_fit_result
   use lavaan_ordinal, only : ordinal_thresholds, bvn_rectangle
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_linalg, only : inverse_general
   use numderiv, only : hessian, nd_success
   implicit none
   private
   public :: fit_ram_pml_ordinal
contains

   subroutine fit_ram_pml_ordinal(template,map,data,result)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      integer,intent(in)::data(:,:)
      type(sem_fit_result),intent(out)::result
      type(ram_model)::work
      integer::n,p,k,i,j,r,info,status,maxcat,iter
      integer,allocatable::ncat(:),counts(:),tab(:,:)
      real(dp),allocatable::threshold(:,:),x(:),hess(:,:),hi(:,:)
      real(dp)::fval,llsat,pr,gradmax,hg
      logical::conv

      n=size(data,1)
      p=size(data,2)
      x=ram_get_free(template,map)
      k=size(x)
      if(n<2 .or. p<2 .or. minval(data)<1) then
      result%status=-1
      return
      end if
      allocate(ncat(p))
      do j=1,p
      ncat(j)=maxval(data(:,j))
      if(ncat(j)<2) then
      result%status=-2
      return
      end if
      end do
      maxcat=maxval(ncat)
      allocate(threshold(maxcat-1,p))
      threshold=0.0_dp
      do j=1,p
         allocate(counts(ncat(j)))
         counts=0
         do r=1,n
         counts(data(r,j))=counts(data(r,j))+1
         end do
         block
            real(dp),allocatable::th(:)
            th=ordinal_thresholds(counts)
            if(size(th)>0) threshold(1:size(th),j)=th
         end block
         deallocate(counts)
      end do

      call bfgs_minimize(nll,x,fval,conv,iter,maxiter=1400,tol=1.0e-7_dp)
      gradmax=0.0_dp
      block
         real(dp),allocatable :: xp(:),xm(:)
         allocate(xp(k),xm(k))
         do i=1,k
            hg=1.0e-4_dp*max(1.0_dp,abs(x(i)))
            xp=x
            xm=x
            xp(i)=xp(i)+hg
            xm(i)=xm(i)-hg
            gradmax=max(gradmax,abs((nll(xp)-nll(xm))/(2.0_dp*hg)))
         end do
      end block
      result%converged=conv .or. gradmax<1.0e-5_dp
      result%iterations=iter
      result%par=x
      result%loglik=-fval
      work=template
      call ram_set_free(work,map,x)
      call ram_sigma(work,result%sigma,info)
      allocate(result%mu(p))
      result%mu=0.0_dp
      if(info/=0) then
      result%status=info
      return
      end if
      llsat=0.0_dp
      result%df=-real(k,dp)
      do j=1,p-1
         do i=j+1,p
            allocate(tab(ncat(j),ncat(i)))
            tab=0
            do r=1,n
            tab(data(r,j),data(r,i))=tab(data(r,j),data(r,i))+1
            end do
            do r=1,ncat(j)
               do status=1,ncat(i)
                  if(tab(r,status)>0) then
                     pr=real(tab(r,status),dp)/real(n,dp)
                     llsat=llsat+real(tab(r,status),dp)*log(pr)
                  end if
               end do
            end do
            result%df=result%df+real((ncat(j)-1)*(ncat(i)-1),dp)
            deallocate(tab)
         end do
      end do
      result%chisq=max(0.0_dp,2.0_dp*(llsat-result%loglik))
      result%objective=2.0_dp*fval/real(n,dp)
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(n,dp))*real(k,dp)
      call hessian(nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
            result%vcov=hi
            do i=1,k
            if(hi(i,i)>=0.0_dp) result%se(i)=sqrt(hi(i,i))
            end do
         end if
      end if
      result%status=0

   contains
      function nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v,rho,l1,u1,l2,u2,pp
         real(dp),allocatable::s(:,:)
         integer::a,b,rr,istat
         work=template
         call ram_set_free(work,map,z)
         call ram_sigma(work,s,istat)
         if(istat/=0 .or. any([(s(a,a)<=0.0_dp,a=1,size(s,1))])) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         v=0.0_dp
         do a=1,p-1
            do b=a+1,p
               rho=s(a,b)/sqrt(s(a,a)*s(b,b))
               if(abs(rho)>=0.9999_dp) then
               v=huge(1.0_dp)/100.0_dp
               return
               end if
               do rr=1,n
                  if(data(rr,a)==1) then
                  l1=-huge(1.0_dp)
                  else
                  l1=threshold(data(rr,a)-1,a)
                  end if
                  if(data(rr,a)==ncat(a)) then
                  u1=huge(1.0_dp)
                  else
                  u1=threshold(data(rr,a),a)
                  end if
                  if(data(rr,b)==1) then
                  l2=-huge(1.0_dp)
                  else
                  l2=threshold(data(rr,b)-1,b)
                  end if
                  if(data(rr,b)==ncat(b)) then
                  u2=huge(1.0_dp)
                  else
                  u2=threshold(data(rr,b),b)
                  end if
                  pp=bvn_rectangle(l1,u1,l2,u2,rho)
                  v=v-log(max(pp,1.0e-300_dp))
               end do
            end do
         end do
      end function nll
   end subroutine fit_ram_pml_ordinal
end module lavaan_pml
