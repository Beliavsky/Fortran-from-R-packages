module lmtest_nonnested
   use lmtest_kinds, only : dp
   use lmtest_types, only : lm_result, test_result, pair_test_result, nonnested_test_result
   use lmtest_regression, only : lm_fit
   use lmtest_inference, only : nested_linear_test
   use lmtest_distributions, only : normal_sf, student_t_sf
   use lmtest_linalg, only : independent_union
   implicit none
   private
   public :: cox_test, j_test, pe_test, encompassing_test

contains

   function cox_test(y, x, z) result(out)
      real(dp), intent(in) :: y(:), x(:,:), z(:,:)
      type(nonnested_test_result) :: out
      type(lm_result) :: xfit, zfit, zxfit, xzfit, xzxfit, zxzfit
      real(dp) :: s2x, s2z, s2zx, s2xz, s2xzx, s2zxz, nreal

      if (size(x,1) /= size(y) .or. size(z,1) /= size(y)) return
      xfit = lm_fit(x,y)
      zfit = lm_fit(z,y)
      if (xfit%info /= 0 .or. zfit%info /= 0) return
      s2x = xfit%rss/real(size(y),dp)
      s2z = zfit%rss/real(size(y),dp)
      zxfit = lm_fit(z,xfit%fitted)
      xzfit = lm_fit(x,zfit%fitted)
      s2zx = s2x + zxfit%rss/real(size(y),dp)
      s2xz = s2z + xzfit%rss/real(size(y),dp)
      xzxfit = lm_fit(x,zxfit%residuals)
      zxzfit = lm_fit(z,xzfit%residuals)
      s2xzx = xzxfit%rss/real(size(y),dp)
      s2zxz = zxzfit%rss/real(size(y),dp)
      nreal = real(size(y),dp)

      out%estimate(1) = 0.5_dp*nreal*log(s2z/s2zx)
      out%std_error(1) = sqrt(nreal*(s2x/(s2zx*s2zx))*s2xzx)
      out%estimate(2) = 0.5_dp*nreal*log(s2x/s2xz)
      out%std_error(2) = sqrt(nreal*(s2z/(s2xz*s2xz))*s2zxz)
      where (out%std_error > 0.0_dp)
         out%statistic = out%estimate/out%std_error
      elsewhere
         out%statistic = 0.0_dp
      end where
      out%p_value(1) = 2.0_dp*normal_sf(abs(out%statistic(1)))
      out%p_value(2) = 2.0_dp*normal_sf(abs(out%statistic(2)))
   end function cox_test

   function j_test(y, x, z) result(out)
      real(dp), intent(in) :: y(:), x(:,:), z(:,:)
      type(nonnested_test_result) :: out
      type(lm_result) :: fitx, fitz, auxx, auxz
      real(dp), allocatable :: ax(:,:), az(:,:)
      integer :: n, px, pz

      n=size(y)
      px=size(x,2)
      pz=size(z,2)
      if(size(x,1)/=n .or. size(z,1)/=n)return
      fitx=lm_fit(x,y)
      fitz=lm_fit(z,y)
      if(fitx%info/=0 .or. fitz%info/=0)return
      allocate(ax(n,px+1),az(n,pz+1))
      ax(:,1:px)=x
      ax(:,px+1)=fitz%fitted
      az(:,1:pz)=z
      az(:,pz+1)=fitx%fitted
      auxx=lm_fit(ax,y)
      auxz=lm_fit(az,y)
      if(auxx%info/=0 .or. auxz%info/=0)return
      out%estimate=[auxx%beta(px+1),auxz%beta(pz+1)]
      out%std_error=[sqrt(max(0.0_dp,auxx%vcov(px+1,px+1))), &
                     sqrt(max(0.0_dp,auxz%vcov(pz+1,pz+1)))]
      if(out%std_error(1)>0.0_dp)out%statistic(1)=out%estimate(1)/out%std_error(1)
      if(out%std_error(2)>0.0_dp)out%statistic(2)=out%estimate(2)/out%std_error(2)
      out%df=[real(auxx%df_resid,dp),real(auxz%df_resid,dp)]
      out%p_value(1)=2.0_dp*student_t_sf(abs(out%statistic(1)),out%df(1))
      out%p_value(2)=2.0_dp*student_t_sf(abs(out%statistic(2)),out%df(2))
   end function j_test

   function pe_test(y1, x, y2, z, islog1, islog2) result(out)
      real(dp),intent(in)::y1(:),x(:,:),y2(:),z(:,:)
      logical,intent(in)::islog1,islog2
      type(nonnested_test_result)::out
      type(lm_result)::f1,f2,a1,a2
      real(dp),allocatable::auxv1(:),auxv2(:),ax(:,:),az(:,:)
      integer::n,px,pz
      n=size(y1)
      px=size(x,2)
      pz=size(z,2)
      if(size(y2)/=n .or. size(x,1)/=n .or. size(z,1)/=n)return
      f1=lm_fit(x,y1)
      f2=lm_fit(z,y2)
      if(f1%info/=0 .or. f2%info/=0)return
      allocate(auxv1(n),auxv2(n))
      if(islog1)then
         auxv1=f2%fitted-exp(f1%fitted)
      else
         if(any(f1%fitted<=0.0_dp))return
         auxv1=log(f1%fitted)-f2%fitted
      end if
      if(islog2)then
         auxv2=f1%fitted-exp(f2%fitted)
      else
         if(any(f2%fitted<=0.0_dp))return
         auxv2=log(f2%fitted)-f1%fitted
      end if
      allocate(ax(n,px+1),az(n,pz+1))
      ax(:,1:px)=x
      ax(:,px+1)=auxv1
      az(:,1:pz)=z
      az(:,pz+1)=auxv2
      a1=lm_fit(ax,y1)
      a2=lm_fit(az,y2)
      if(a1%info/=0 .or. a2%info/=0)return
      out%estimate=[a1%beta(px+1),a2%beta(pz+1)]
      out%std_error=[sqrt(max(0.0_dp,a1%vcov(px+1,px+1))),sqrt(max(0.0_dp,a2%vcov(pz+1,pz+1)))]
      if(out%std_error(1)>0.0_dp)out%statistic(1)=out%estimate(1)/out%std_error(1)
      if(out%std_error(2)>0.0_dp)out%statistic(2)=out%estimate(2)/out%std_error(2)
      out%df=[real(a1%df_resid,dp),real(a2%df_resid,dp)]
      out%p_value(1)=2.0_dp*student_t_sf(abs(out%statistic(1)),out%df(1))
      out%p_value(2)=2.0_dp*student_t_sf(abs(out%statistic(2)),out%df(2))
   end function pe_test

   function encompassing_test(y, x, z, use_f) result(out)
      real(dp),intent(in)::y(:),x(:,:),z(:,:)
      logical,intent(in),optional::use_f
      type(pair_test_result)::out
      real(dp),allocatable::u(:,:)
      logical::as_f
      as_f=.true.
      if(present(use_f))as_f=use_f
      call independent_union(x,z,u)
      if(size(u,1)==0)return
      out%first=nested_linear_test(y,u,x,as_f)
      out%second=nested_linear_test(y,u,z,as_f)
   end function encompassing_test

end module lmtest_nonnested
