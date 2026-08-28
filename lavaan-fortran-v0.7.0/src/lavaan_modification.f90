module lavaan_modification
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_objectives, only : objective_ml
   use lavaan_linalg, only : inverse_general
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: modification_result
      real(dp), allocatable :: mi(:), epc(:), score(:), efficient_info(:, :)
      real(dp) :: joint_score = 0.0_dp
      integer :: status = 0
   end type modification_result

   public :: modification_indices_cov

contains

   subroutine modification_indices_cov(template,free_map,candidate_map,data_cov,data_mean,nobs,result)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::free_map,candidate_map
      real(dp),intent(in)::data_cov(:,:),data_mean(:)
      integer,intent(in)::nobs
      type(modification_result),intent(out)::result
      type(ram_model)::work
      type(ram_free_map)::joint
      real(dp),allocatable::xf(:),xc(:),x(:),grad(:),hess(:,:),hff(:,:),hfc(:,:),hcf(:,:),hcc(:,:)
      real(dp),allocatable::hffi(:,:),effi(:,:),eiinv(:,:),xp(:),xm(:)
      real(dp)::fval,h
      integer::kf,kc,k,i,info,status,it
      logical::conv

      xf=ram_get_free(template,free_map)
      xc=ram_get_free(template,candidate_map)
      kf=size(xf)
      kc=size(xc)
      k=kf+kc
      call combine_maps(free_map,candidate_map,joint)
      call bfgs_minimize(free_only_objective,xf,fval,conv,it,maxiter=1000,tol=1.0e-8_dp)
      x=[xf,xc]
      allocate(grad(k),xp(k),xm(k))
      grad=0.0_dp
      do i=1,k
         h=1.0e-5_dp*max(1.0_dp,abs(x(i)))
         xp=x
         xm=x
         xp(i)=xp(i)+h
         xm(i)=xm(i)-h
         grad(i)=(joint_objective(xp)-joint_objective(xm))/(2.0_dp*h)
      end do
      call hessian(joint_objective,x,hess,status=status)
      allocate(result%mi(kc),result%epc(kc),result%score(kc),result%efficient_info(kc,kc))
      result%mi=0.0_dp
      result%epc=0.0_dp
      result%score=-grad(kf+1:k)
      if(status/=nd_success) then
      result%status=1
      return
      end if
      if(kf>0) then
         hff=hess(1:kf,1:kf)
         hfc=hess(1:kf,kf+1:k)
         hcf=hess(kf+1:k,1:kf)
         hcc=hess(kf+1:k,kf+1:k)
         call inverse_general(hff,hffi,info)
         if(info/=0) then
         result%status=info
         return
         end if
         effi=hcc-matmul(hcf,matmul(hffi,hfc))
      else
         effi=hess
      end if
      result%efficient_info=effi
      do i=1,kc
         if(effi(i,i)>1.0e-12_dp) then
            result%mi(i)=result%score(i)**2/effi(i,i)
            result%epc(i)=result%score(i)/effi(i,i)
         end if
      end do
      call inverse_general(effi,eiinv,info)
      if(info==0) result%joint_score=dot_product(result%score,matmul(eiinv,result%score))
      result%status=0

   contains
      function free_only_objective(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::sg(:,:),mu(:)
         integer::istat,jj
         work=template
         call ram_set_free(work,free_map,z)
         call ram_sigma(work,sg,istat)
         if(istat/=0 .or. any([(sg(jj,jj)<=0.0_dp,jj=1,size(sg,1))])) then
            v=huge(1.0_dp)/100.0_dp
            return
         end if
         call ram_mu(work,mu,istat)
         v=0.5_dp*real(nobs,dp)*objective_ml(sg,mu,data_cov,data_mean,allocated(template%m),istat)
      end function free_only_objective
      function joint_objective(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::sg(:,:),mu(:)
         integer::istat,jj
         work=template
         call ram_set_free(work,joint,z)
         call ram_sigma(work,sg,istat)
         if(istat/=0 .or. any([(sg(jj,jj)<=0.0_dp,jj=1,size(sg,1))])) then
            v=huge(1.0_dp)/100.0_dp
            return
         end if
         call ram_mu(work,mu,istat)
         v=0.5_dp*real(nobs,dp)*objective_ml(sg,mu,data_cov,data_mean,allocated(template%m),istat)
      end function joint_objective
   end subroutine modification_indices_cov

   subroutine combine_maps(a,b,c)
      type(ram_free_map),intent(in)::a,b
      type(ram_free_map),intent(out)::c
      integer::na,nb
      na=size(a%matrix_id)
      nb=size(b%matrix_id)
      allocate(c%matrix_id(na+nb),c%row(na+nb),c%col(na+nb))
      if(na>0) then
         c%matrix_id(1:na)=a%matrix_id
         c%row(1:na)=a%row
         c%col(1:na)=a%col
      end if
      if(nb>0) then
         c%matrix_id(na+1:na+nb)=b%matrix_id
         c%row(na+1:na+nb)=b%row
         c%col(na+1:na+nb)=b%col
      end if
   end subroutine combine_maps
end module lavaan_modification
