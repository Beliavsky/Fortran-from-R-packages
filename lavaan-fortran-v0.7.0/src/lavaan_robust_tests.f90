module lavaan_robust_tests
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_set_free, ram_sigma
   use lavaan_fit, only : sem_fit_result
   use lavaan_linalg, only : vech, inverse_general, trace_matrix
   implicit none
   private

   type, public :: scaled_test_result
      real(dp) :: chisq_standard=huge(1.0_dp), df_standard=0.0_dp
      real(dp) :: sb_scaling=1.0_dp, chisq_sb=huge(1.0_dp), df_sb=0.0_dp
      real(dp) :: mv_scaling=1.0_dp, chisq_mv=huge(1.0_dp), df_mv=0.0_dp
      real(dp) :: ss_scaling=1.0_dp, ss_shift=0.0_dp, chisq_ss=huge(1.0_dp), df_ss=0.0_dp
      real(dp) :: yb_scaling=1.0_dp, chisq_yb=huge(1.0_dp), yb_h1_scaling=huge(1.0_dp)
      real(dp) :: yb_h0_scaling=huge(1.0_dp)
      real(dp) :: trace1=0.0_dp, trace2=0.0_dp
      integer :: status=0
   end type scaled_test_result

   public :: covariance_scaled_tests, scaled_tests_from_ugamma, yuan_bentler_from_traces

contains

   subroutine covariance_scaled_tests(template,map,data,fit,result)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      type(sem_fit_result),intent(in)::fit
      type(scaled_test_result),intent(out)::result
      type(ram_model)::work
      real(dp),allocatable::gamma_emp(:,:),gamma_norm(:,:),w(:,:),delta(:,:),mid(:,:),midi(:,:),u(:,:),ug(:,:)
      real(dp),allocatable::z(:,:),zbar(:),d(:),v0(:),xp(:),xm(:),sp(:,:),sm(:,:),vp(:),vm(:)
      integer::n,p,q,k,i,j,ii,jj,row,col,a,info
      real(dp)::h
      n=size(data,1)
      p=size(data,2)
      q=p*(p+1)/2
      k=size(fit%par)
      if(size(data,2)/=size(fit%sigma,1) .or. fit%df<=0.0_dp) then
      result%status=-1
      return
      end if
      allocate(z(n,q),zbar(q),d(p))
      z=0.0_dp
      do i=1,n
         d=data(i,:)-sum(data,dim=1)/real(n,dp)
         z(i,:)=vech(spread(d,2,p)*spread(d,1,p))
      end do
      zbar=sum(z,dim=1)/real(n,dp)
      allocate(gamma_emp(q,q))
      gamma_emp=0.0_dp
      do i=1,n
         v0=z(i,:)-zbar
         gamma_emp=gamma_emp+spread(v0,2,q)*spread(v0,1,q)
      end do
      gamma_emp=gamma_emp/real(n,dp)
      allocate(gamma_norm(q,q))
      gamma_norm=0.0_dp
      row=0
      do j=1,p
      do i=j,p
      row=row+1
      col=0
         do jj=1,p
         do ii=jj,p
         col=col+1
            gamma_norm(row,col)=fit%sigma(i,ii)*fit%sigma(j,jj)+fit%sigma(i,jj)*fit%sigma(j,ii)
         end do
         end do
      end do
      end do
      call inverse_general(gamma_norm,w,info)
      if(info/=0) then
      result%status=info
      return
      end if
      allocate(delta(q,k),xp(k),xm(k))
      work=template
      do a=1,k
         h=1.0e-5_dp*max(1.0_dp,abs(fit%par(a)))
         xp=fit%par
         xm=fit%par
         xp(a)=xp(a)+h
         xm(a)=xm(a)-h
         work=template
         call ram_set_free(work,map,xp)
         call ram_sigma(work,sp,info)
         if(info/=0) then
            result%status=info
            return
         end if
         work=template
         call ram_set_free(work,map,xm)
         call ram_sigma(work,sm,info)
         if(info/=0) then
            result%status=info
            return
         end if
         vp=vech(sp)
         vm=vech(sm)
         delta(:,a)=(vp-vm)/(2.0_dp*h)
      end do
      mid=matmul(transpose(delta),matmul(w,delta))
      call inverse_general(mid,midi,info)
      if(info/=0) then
      result%status=info
      return
      end if
      u=w-matmul(w,matmul(delta,matmul(midi,matmul(transpose(delta),w))))
      ug=matmul(u,gamma_emp)
      call scaled_tests_from_ugamma(fit%chisq,fit%df,ug,result)
   end subroutine covariance_scaled_tests

   subroutine scaled_tests_from_ugamma(chisq,df,ug,result)
      real(dp),intent(in)::chisq,df,ug(:,:)
      type(scaled_test_result),intent(out)::result
      real(dp)::tr1,tr2
      if(size(ug,1)/=size(ug,2) .or. df<=0.0_dp) then
      result%status=-1
      return
      end if
      tr1=trace_matrix(ug)
      tr2=trace_matrix(matmul(ug,ug))
      result%chisq_standard=chisq
      result%df_standard=df
      result%trace1=tr1
      result%trace2=tr2
      if(tr1<=1.0e-14_dp .or. tr2<=1.0e-14_dp) then
      result%status=-2
      return
      end if
      result%sb_scaling=tr1/df
      result%chisq_sb=chisq/result%sb_scaling
      result%df_sb=df
      result%mv_scaling=tr2/tr1
      result%chisq_mv=chisq/result%mv_scaling
      result%df_mv=tr1*tr1/tr2
      result%df_ss=df
      result%ss_scaling=sqrt(tr2/df)
      result%ss_shift=df-sqrt(df/tr2)*tr1
      result%chisq_ss=chisq/result%ss_scaling+result%ss_shift
      result%yb_scaling=result%sb_scaling
      result%chisq_yb=result%chisq_sb
      result%status=0
   end subroutine scaled_tests_from_ugamma

   subroutine yuan_bentler_from_traces(chisq,df,trace_ugamma,h1_trace,h1_dim,h0_trace,h0_dim,result)
      real(dp),intent(in)::chisq,df,trace_ugamma,h1_trace,h1_dim,h0_trace,h0_dim
      type(scaled_test_result),intent(inout)::result
      if(df<=0.0_dp .or. trace_ugamma<=0.0_dp) then
      result%status=-1
      return
      end if
      result%chisq_standard=chisq
      result%df_standard=df
      result%trace1=trace_ugamma
      result%yb_scaling=trace_ugamma/df
      result%chisq_yb=chisq/result%yb_scaling
      if(h1_dim>0.0_dp) result%yb_h1_scaling=h1_trace/h1_dim
      if(h0_dim>0.0_dp) result%yb_h0_scaling=h0_trace/h0_dim
      result%status=0
   end subroutine yuan_bentler_from_traces

end module lavaan_robust_tests
