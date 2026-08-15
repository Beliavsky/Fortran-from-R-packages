! Additional numerical diagnostics for matrix-first GAMLSS fits.
! Includes worm-plot polynomial summaries and leverage-based influence measures.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_diagnostics_v04
   use gamlss_kinds, only : dp
   use gamlss_special, only : normal_quantile
   use gamlss_linalg, only : solve_linear
   implicit none
   private
   public :: worm_result_t, worm_plot_diagnostics, influence_from_hat, jarque_bera_statistic

   type,public :: worm_result_t
      real(dp),allocatable :: bin_center(:)
      real(dp),allocatable :: coefficients(:,:)
      integer,allocatable :: n_bin(:)
      integer :: status=0
   end type worm_result_t
contains

   subroutine worm_plot_diagnostics(residuals,covariate,n_intervals,result,status)
      real(dp),intent(in)::residuals(:),covariate(:)
      integer,intent(in)::n_intervals
      type(worm_result_t),intent(out)::result
      integer,intent(out),optional::status
      integer,allocatable::idx(:)
      real(dp),allocatable::r(:),z(:),d(:),a(:,:),rhs(:),coef(:)
      integer::n,b,lo,hi,m,i,j,k,istat
      n=size(residuals);istat=0
      if(size(covariate)/=n.or.n<4.or.n_intervals<1.or.n_intervals>n/4)then
         result%status=1;if(present(status))status=1;return
      end if
      call order_index(covariate,idx)
      allocate(result%bin_center(n_intervals),result%coefficients(4,n_intervals),result%n_bin(n_intervals))
      result%coefficients=0.0_dp
      do b=1,n_intervals
         lo=1+(b-1)*n/n_intervals
         hi=b*n/n_intervals
         if(b==n_intervals)hi=n
         m=hi-lo+1;result%n_bin(b)=m
         result%bin_center(b)=sum(covariate(idx(lo:hi)))/real(m,dp)
         allocate(r(m),z(m),d(m));r=residuals(idx(lo:hi));call sort_real(r)
         do i=1,m
            z(i)=normal_quantile((real(i,dp)-0.375_dp)/(real(m,dp)+0.25_dp))
         end do
         d=r-z;allocate(a(4,4),rhs(4));a=0.0_dp;rhs=0.0_dp
         do i=1,m
            do j=1,4
               rhs(j)=rhs(j)+d(i)*z(i)**(j-1)
               do k=1,4;a(j,k)=a(j,k)+z(i)**(j+k-2);end do
            end do
         end do
         call solve_linear(a,rhs,coef,istat)
         if(istat==0.and.size(coef)==4)then
            result%coefficients(:,b)=coef
         else
            result%status=2
         end if
         deallocate(r,z,d,a,rhs)
         if(allocated(coef))deallocate(coef)
      end do
      if(result%status==0)result%status=istat
      if(present(status))status=result%status
   end subroutine worm_plot_diagnostics

   subroutine influence_from_hat(residuals,hat,n_parameters,studentized,cooks_distance,status)
      real(dp),intent(in)::residuals(:),hat(:)
      integer,intent(in)::n_parameters
      real(dp),allocatable,intent(out)::studentized(:),cooks_distance(:)
      integer,intent(out),optional::status
      real(dp)::den
      integer::i,n
      n=size(residuals)
      if(size(hat)/=n.or.n_parameters<1)then
         allocate(studentized(0),cooks_distance(0));if(present(status))status=1;return
      end if
      allocate(studentized(n),cooks_distance(n))
      do i=1,n
         den=max(1.0e-12_dp,1.0_dp-min(0.999999999_dp,max(0.0_dp,hat(i))))
         studentized(i)=residuals(i)/sqrt(den)
         cooks_distance(i)=residuals(i)**2*max(0.0_dp,hat(i))/(real(n_parameters,dp)*den*den)
      end do
      if(present(status))status=0
   end subroutine influence_from_hat

   real(dp) function jarque_bera_statistic(residuals) result(jb)
      real(dp),intent(in)::residuals(:)
      real(dp)::m,s2,s3,s4,n
      n=real(size(residuals),dp)
      if(size(residuals)<3)then;jb=0.0_dp;return;end if
      m=sum(residuals)/n;s2=sum((residuals-m)**2)/n
      if(s2<=tiny(1.0_dp))then;jb=0.0_dp;return;end if
      s3=sum((residuals-m)**3)/n/s2**1.5_dp
      s4=sum((residuals-m)**4)/n/(s2*s2)-3.0_dp
      jb=n*(s3*s3/6.0_dp+s4*s4/24.0_dp)
   end function jarque_bera_statistic

   subroutine order_index(x,idx)
      real(dp),intent(in)::x(:)
      integer,allocatable,intent(out)::idx(:)
      integer::i,j,key
      allocate(idx(size(x)));idx=[(i,i=1,size(x))]
      do i=2,size(x)
         key=idx(i);j=i-1
         do while(j>=1)
            if(x(idx(j))<=x(key))exit
            idx(j+1)=idx(j);j=j-1
         end do
         idx(j+1)=key
      end do
   end subroutine order_index

   subroutine sort_real(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::i,j
      do i=2,size(x)
         key=x(i);j=i-1
         do while(j>=1)
            if(x(j)<=key)exit
            x(j+1)=x(j);j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real
end module gamlss_diagnostics_v04
