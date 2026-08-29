module lavaan_robust
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_set_free, ram_sigma, ram_mu
   use lavaan_fit, only : sem_fit_result
   use lavaan_linalg, only : inverse_general, inverse_spd, logdet_spd, vech, trace_matrix
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: robust_sem_result
      real(dp), allocatable :: score(:, :), vcov(:, :), se(:)
      real(dp) :: sb_scaling = 1.0_dp
      real(dp) :: chisq_scaled = huge(1.0_dp)
      real(dp) :: df = 0.0_dp
      integer :: status = 0, ncluster = 0
   end type robust_sem_result

   public :: robust_ml_inference

contains

   subroutine robust_ml_inference(template,map,data,fit,result,cluster)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      type(sem_fit_result),intent(in)::fit
      type(robust_sem_result),intent(out)::result
      integer,intent(in),optional::cluster(:)
      real(dp),allocatable::hess(:,:),hinv(:,:),bmat(:,:),sg(:,:),mu(:),sagg(:,:),g(:)
      real(dp),allocatable::gamma_emp(:,:),gamma_norm(:,:),w(:,:),delta(:,:),mid(:,:),midi(:,:),u(:,:)
      real(dp),allocatable::z(:,:),zbar(:),d(:),vv(:),xp(:),xm(:),sp(:,:),sm(:,:),vvp(:),vvm(:)
      integer,allocatable::clabel(:)
      integer::n,p,k,i,j,a,c,info,status,nc,ci,qidx,row,ii,jj
      real(dp)::h,correction
      type(ram_model)::work

      n=size(data,1)
      p=size(data,2)
      k=size(fit%par)
      if(any(shape(data)/=[n,p]) .or. k/=size(map%matrix_id)) then
      result%status=-1
      return
      end if
      if(present(cluster)) then
         if(size(cluster)/=n) then
         result%status=-2
         return
         end if
      end if
      work=template
      call ram_set_free(work,map,fit%par)
      call ram_sigma(work,sg,info)
      call ram_mu(work,mu,info)
      if(info/=0) then
      result%status=info
      return
      end if

      allocate(result%score(n,k),g(k))
      result%score=0.0_dp
      do i=1,n
         call case_score(data(i,:),fit%par,g)
         result%score(i,:)=g
      end do

      call hessian(total_nll,fit%par,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status/=nd_success) then
      result%status=3
      return
      end if
      call inverse_general(hess,hinv,info)
      if(info/=0) then
      result%status=info
      return
      end if
      allocate(bmat(k,k))
      bmat=0.0_dp
      if(present(cluster)) then
         call unique_labels(cluster,clabel)
         nc=size(clabel)
         result%ncluster=nc
         allocate(sagg(nc,k))
         sagg=0.0_dp
         do i=1,n
            ci=label_index(clabel,cluster(i))
            sagg(ci,:)=sagg(ci,:)+result%score(i,:)
         end do
         do i=1,nc
            bmat=bmat+spread(sagg(i,:),2,k)*spread(sagg(i,:),1,k)
         end do
         if(nc>1 .and. n>k) then
            correction=real(nc,dp)/real(nc-1,dp)*real(n-1,dp)/real(n-k,dp)
            bmat=bmat*correction
         end if
      else
         result%ncluster=n
         do i=1,n
            bmat=bmat+spread(result%score(i,:),2,k)*spread(result%score(i,:),1,k)
         end do
      end if
      result%vcov=matmul(hinv,matmul(bmat,hinv))
      do i=1,k
         if(result%vcov(i,i)>=0.0_dp) result%se(i)=sqrt(result%vcov(i,i))
      end do

      ! Satorra-Bentler scaling for covariance-structure ML.
      qidx=p*(p+1)/2
      allocate(z(n,qidx),zbar(qidx),d(p))
      do i=1,n
         d=data(i,:)-sum(data,dim=1)/real(n,dp)
         z(i,:)=vech(spread(d,2,p)*spread(d,1,p))
      end do
      zbar=sum(z,dim=1)/real(n,dp)
      allocate(gamma_emp(qidx,qidx))
      gamma_emp=0.0_dp
      do i=1,n
         vv=z(i,:)-zbar
         gamma_emp=gamma_emp+spread(vv,2,qidx)*spread(vv,1,qidx)
      end do
      gamma_emp=gamma_emp/real(n,dp)
      allocate(gamma_norm(qidx,qidx))
      gamma_norm=0.0_dp
      row=0
      do j=1,p
         do i=j,p
            row=row+1
            c=0
            do jj=1,p
               do ii=jj,p
                  c=c+1
                  gamma_norm(row,c)=sg(i,ii)*sg(j,jj)+sg(i,jj)*sg(j,ii)
               end do
            end do
         end do
      end do
      call inverse_general(gamma_norm,w,info)
      if(info==0 .and. fit%df>0.0_dp) then
         allocate(delta(qidx,k),xp(k),xm(k))
         do a=1,k
            h=1.0e-5_dp*max(1.0_dp,abs(fit%par(a)))
            xp=fit%par
            xm=fit%par
            xp(a)=xp(a)+h
            xm(a)=xm(a)-h
            work=template
            call ram_set_free(work,map,xp)
            call ram_sigma(work,sp,info)
            vvp=vech(sp)
            work=template
            call ram_set_free(work,map,xm)
            call ram_sigma(work,sm,info)
            vvm=vech(sm)
            delta(:,a)=(vvp-vvm)/(2.0_dp*h)
         end do
         mid=matmul(transpose(delta),matmul(w,delta))
         call inverse_general(mid,midi,info)
         if(info==0) then
            u=w-matmul(w,matmul(delta,matmul(midi,matmul(transpose(delta),w))))
            result%sb_scaling=trace_matrix(matmul(u,gamma_emp))/fit%df
            if(result%sb_scaling>1.0e-12_dp) then
               result%chisq_scaled=fit%chisq/result%sb_scaling
            else
               result%sb_scaling=1.0_dp
               result%chisq_scaled=fit%chisq
            end if
         else
            result%chisq_scaled=fit%chisq
         end if
      else
         result%chisq_scaled=fit%chisq
      end if
      result%df=fit%df
      result%status=0

   contains
      function case_loglik(obs,zpar) result(ll)
         real(dp),intent(in)::obs(:),zpar(:)
         real(dp)::ll,ld,quad
         real(dp),allocatable::ss(:,:),mmi(:),sii(:,:),dd(:)
         integer::istat
         work=template
         call ram_set_free(work,map,zpar)
         call ram_sigma(work,ss,istat)
         if(istat/=0) then
         ll=-huge(1.0_dp)/100.0_dp
         return
         end if
         call ram_mu(work,mmi,istat)
         call inverse_spd(ss,sii,istat)
         if(istat/=0) then
         ll=-huge(1.0_dp)/100.0_dp
         return
         end if
         ld=logdet_spd(ss,istat)
         dd=obs-mmi
         quad=dot_product(dd,matmul(sii,dd))
         ll=-0.5_dp*(real(size(obs),dp)*log(2.0_dp*acos(-1.0_dp))+ld+quad)
      end function case_loglik
      subroutine case_score(obs,zpar,sc)
         real(dp),intent(in)::obs(:),zpar(:)
         real(dp),intent(out)::sc(:)
         real(dp),allocatable::zp(:),zm(:)
         real(dp)::hh
         integer::kk
         allocate(zp(size(zpar)),zm(size(zpar)))
         do kk=1,size(zpar)
            hh=1.0e-5_dp*max(1.0_dp,abs(zpar(kk)))
            zp=zpar
            zm=zpar
            zp(kk)=zp(kk)+hh
            zm(kk)=zm(kk)-hh
            sc(kk)=(case_loglik(obs,zp)-case_loglik(obs,zm))/(2.0_dp*hh)
         end do
      end subroutine case_score
      function total_nll(zpar) result(v)
         real(dp),intent(in)::zpar(:)
         real(dp)::v
         integer::rr
         v=0.0_dp
         do rr=1,n
         v=v-case_loglik(data(rr,:),zpar)
         end do
      end function total_nll
   end subroutine robust_ml_inference

   subroutine unique_labels(x,u)
      integer,intent(in)::x(:)
      integer,allocatable,intent(out)::u(:)
      integer,allocatable::tmp(:)
      integer::i,m
      allocate(tmp(size(x)))
      m=0
      do i=1,size(x)
         if(m==0 .or. .not.any(tmp(1:m)==x(i))) then
         m=m+1
         tmp(m)=x(i)
         end if
      end do
      allocate(u(m))
      u=tmp(1:m)
   end subroutine unique_labels

   integer function label_index(u,x) result(idx)
      integer,intent(in)::u(:),x
      integer::i
      idx=0
      do i=1,size(u)
      if(u(i)==x) then
      idx=i
      return
      end if
      end do
   end function label_index
end module lavaan_robust
