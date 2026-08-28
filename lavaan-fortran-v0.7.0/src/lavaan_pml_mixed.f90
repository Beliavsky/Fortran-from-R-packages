module lavaan_pml_mixed
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_fit, only : sem_fit_result
   use lavaan_ordinal, only : ordinal_thresholds, bvn_rectangle
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_linalg, only : inverse_general
   use numderiv, only : hessian, nd_success
   implicit none
   private
   public :: fit_ram_pml_mixed
contains
   subroutine fit_ram_pml_mixed(template,map,data,is_ordinal,result)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      logical,intent(in)::is_ordinal(:)
      type(sem_fit_result),intent(out)::result
      type(ram_model)::work
      integer::n,p,k,j,r,maxcat,info,status,i,iter
      integer,allocatable::ncat(:),cnt(:)
      real(dp),allocatable::threshold(:,:),x(:),hess(:,:),hi(:,:)
      real(dp)::fval
      logical::conv
      n=size(data,1)
      p=size(data,2)
      if(size(is_ordinal)/=p .or. n<2 .or. p<2) then
      result%status=-1
      return
      end if
      allocate(ncat(p))
      ncat=0
      maxcat=1
      do j=1,p
         if(is_ordinal(j)) then
            ncat(j)=nint(maxval(data(:,j)))
            if(ncat(j)<2 .or. minval(data(:,j))<1.0_dp) then
            result%status=-2
            return
            end if
            if(maxval(abs(data(:,j)-real(nint(data(:,j)),dp)))>1.0e-10_dp) then
            result%status=-3
            return
            end if
            maxcat=max(maxcat,ncat(j))
         end if
      end do
      allocate(threshold(max(1,maxcat-1),p))
      threshold=0.0_dp
      do j=1,p
         if(is_ordinal(j)) then
            allocate(cnt(ncat(j)))
            cnt=0
            do r=1,n
            cnt(nint(data(r,j)))=cnt(nint(data(r,j)))+1
            end do
            if(any(cnt==0)) then
            result%status=-4
            return
            end if
            block
               real(dp),allocatable::th(:)
               th=ordinal_thresholds(cnt)
               threshold(1:size(th),j)=th
            end block
            deallocate(cnt)
         end if
      end do
      x=ram_get_free(template,map)
      k=size(x)
      call bfgs_minimize(nll,x,fval,conv,iter,maxiter=1600,tol=1.0e-7_dp)
      result%par=x
      result%converged=conv
      result%iterations=iter
      result%loglik=-fval
      result%objective=2.0_dp*fval/real(n,dp)
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(n,dp))*real(k,dp)
      work=template
      call ram_set_free(work,map,x)
      call ram_sigma(work,result%sigma,info)
      call ram_mu(work,result%mu,info)
      if(info/=0) then
      result%status=info
      return
      end if
      result%chisq=huge(1.0_dp)
      result%df=0.0_dp
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
         real(dp)::v
         real(dp),allocatable::sg(:,:),mu(:)
         real(dp)::rho,s1,s2,z1,z2,pr,ll,l,u,condmu,condsd
         integer::a,b,rr,ca,istat
         work=template
         call ram_set_free(work,map,z)
         call ram_sigma(work,sg,istat)
         call ram_mu(work,mu,istat)
         if(istat/=0 .or. any([(sg(a,a)<=0.0_dp,a=1,p)])) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         v=0.0_dp
         do a=1,p-1
         do b=a+1,p
            rho=sg(a,b)/sqrt(sg(a,a)*sg(b,b))
            if(abs(rho)>=0.9999_dp) then
            v=huge(1.0_dp)/100.0_dp
            return
            end if
            s1=sqrt(sg(a,a))
            s2=sqrt(sg(b,b))
            do rr=1,n
               if(.not.is_ordinal(a) .and. .not.is_ordinal(b)) then
                  z1=(data(rr,a)-mu(a))/s1
                  z2=(data(rr,b)-mu(b))/s2
                  ll=-log(2.0_dp*acos(-1.0_dp))-log(s1*s2)-0.5_dp*log(1.0_dp-rho*rho) &
                     -0.5_dp*(z1*z1-2.0_dp*rho*z1*z2+z2*z2)/(1.0_dp-rho*rho)
                  v=v-ll
               else if(is_ordinal(a) .and. is_ordinal(b)) then
                  ca=nint(data(rr,a))
                  if(ca==1) then
                  l=-huge(1.0_dp)
                  else
                  l=threshold(ca-1,a)
                  end if
                  if(ca==ncat(a)) then
                  u=huge(1.0_dp)
                  else
                  u=threshold(ca,a)
                  end if
                  block
                     integer::cb
                     real(dp)::l2,u2
                     cb=nint(data(rr,b))
                     if(cb==1) then
                     l2=-huge(1.0_dp)
                     else
                     l2=threshold(cb-1,b)
                     end if
                     if(cb==ncat(b)) then
                     u2=huge(1.0_dp)
                     else
                     u2=threshold(cb,b)
                     end if
                     pr=bvn_rectangle(l,u,l2,u2,rho)
                  end block
                  v=v-log(max(pr,1.0e-300_dp))
               else if(.not.is_ordinal(a)) then
                  z1=(data(rr,a)-mu(a))/s1
                  ca=nint(data(rr,b))
                  condmu=rho*z1
                  condsd=sqrt(1.0_dp-rho*rho)
                  if(ca==1) then
                  l=-huge(1.0_dp)
                  else
                  l=(threshold(ca-1,b)-condmu)/condsd
                  end if
                  if(ca==ncat(b)) then
                  u=huge(1.0_dp)
                  else
                  u=(threshold(ca,b)-condmu)/condsd
                  end if
                  pr=max(normal_cdf(u)-normal_cdf(l),1.0e-300_dp)
                  ll=-0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(s1)-0.5_dp*z1*z1+log(pr)
                  v=v-ll
               else
                  z2=(data(rr,b)-mu(b))/s2
                  ca=nint(data(rr,a))
                  condmu=rho*z2
                  condsd=sqrt(1.0_dp-rho*rho)
                  if(ca==1) then
                  l=-huge(1.0_dp)
                  else
                  l=(threshold(ca-1,a)-condmu)/condsd
                  end if
                  if(ca==ncat(a)) then
                  u=huge(1.0_dp)
                  else
                  u=(threshold(ca,a)-condmu)/condsd
                  end if
                  pr=max(normal_cdf(u)-normal_cdf(l),1.0e-300_dp)
                  ll=-0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(s2)-0.5_dp*z2*z2+log(pr)
                  v=v-ll
               end if
            end do
         end do
         end do
      end function nll
   end subroutine fit_ram_pml_mixed

   pure function normal_cdf(z) result(p)
      real(dp),intent(in)::z
      real(dp)::p
      if(z>8.0_dp) then
      p=1.0_dp
      else if(z< -8.0_dp) then
      p=0.0_dp
      else
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      end if
   end function normal_cdf
end module lavaan_pml_mixed
