module rfast_batch
   use rfast_special, only : dp, f_cdf, chisq_cdf
   use rfast_arrays, only : colmeans, colvars, mean_r, variance_r, rank_average
   use rfast_mle, only : mle_result, normal_mle, gamma_mle, laplace_mle, pareto_mle, rayleigh_mle, &
                         poisson_mle, geometric_mle, exponential_mle, exponential2_mle, lognormal_mle, &
                         invgauss_mle, lindley_mle, maxboltz_mle
   use rfast_extra_mle, only : weibull_mle
   use rfast_tests, only : test_result, ttest1, ttest2, ftest_variance, vartest_chisq
   implicit none
   private
   public :: colnormal_mle, colgamma_mle, collaplace_mle, colpareto_mle, colrayleigh_mle
   public :: colpoisson_mle, colgeom_mle
   public :: colexponential_mle, colexponential2_mle, collognormal_mle, colinvgauss_mle
   public :: collindley_mle, colmaxboltz_mle, colweibull_mle
   public :: column_ttests, column_ftests, column_vartests, column_aucs
   public :: one_way_anova, one_way_anovas, correlation_pairs, column_correlations

contains

   function colnormal_mle(x) result(par)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: par(2,size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=normal_mle(x(:,j)); par(:,j)=r%param
      end do
   end function colnormal_mle

   function colgamma_mle(x) result(par)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: par(2,size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=gamma_mle(x(:,j)); par(:,j)=r%param
      end do
   end function colgamma_mle

   function collaplace_mle(x) result(par)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: par(2,size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=laplace_mle(x(:,j)); par(:,j)=r%param
      end do
   end function collaplace_mle

   function colpareto_mle(x) result(par)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: par(2,size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=pareto_mle(x(:,j)); par(:,j)=r%param
      end do
   end function colpareto_mle

   function colrayleigh_mle(x) result(par)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: par(size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=rayleigh_mle(x(:,j)); par(j)=r%param(1)
      end do
   end function colrayleigh_mle

   function colexponential_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=exponential_mle(x(:,j));par(j)=r%param(1);end do
   end function colexponential_mle

   function colexponential2_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(2,size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=exponential2_mle(x(:,j));par(:,j)=r%param;end do
   end function colexponential2_mle

   function collognormal_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(2,size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=lognormal_mle(x(:,j));par(:,j)=r%param;end do
   end function collognormal_mle

   function colinvgauss_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(2,size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=invgauss_mle(x(:,j));par(:,j)=r%param;end do
   end function colinvgauss_mle

   function collindley_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=lindley_mle(x(:,j));par(j)=r%param(1);end do
   end function collindley_mle

   function colmaxboltz_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=maxboltz_mle(x(:,j));par(j)=r%param(1);end do
   end function colmaxboltz_mle

   function colweibull_mle(x) result(par)
      real(dp), intent(in) :: x(:,:);real(dp)::par(2,size(x,2));type(mle_result)::r;integer::j
      do j=1,size(x,2);r=weibull_mle(x(:,j));par(:,j)=r%param;end do
   end function colweibull_mle

   function colpoisson_mle(x) result(par)
      integer, intent(in) :: x(:,:)
      real(dp) :: par(size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=poisson_mle(x(:,j)); par(j)=r%param(1)
      end do
   end function colpoisson_mle

   function colgeom_mle(x) result(par)
      integer, intent(in) :: x(:,:)
      real(dp) :: par(size(x,2))
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=geometric_mle(x(:,j)); par(j)=r%param(1)
      end do
   end function colgeom_mle

   function column_ttests(x,mu) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: mu
      real(dp) :: out(3,size(x,2))
      type(test_result) :: r
      real(dp) :: m0
      integer :: j
      m0=0.0_dp; if(present(mu))m0=mu
      do j=1,size(x,2)
         r=ttest1(x(:,j),m0);out(:,j)=[r%statistic,r%df1,r%pvalue]
      end do
   end function column_ttests

   function column_ftests(x,y) result(out)
      real(dp), intent(in) :: x(:,:),y(:,:)
      real(dp) :: out(4,min(size(x,2),size(y,2)))
      type(test_result) :: r
      integer :: j
      do j=1,size(out,2)
         r=ftest_variance(x(:,j),y(:,j));out(:,j)=[r%statistic,r%df1,r%df2,r%pvalue]
      end do
   end function column_ftests

   function column_vartests(x,var0) result(out)
      real(dp), intent(in) :: x(:,:),var0
      real(dp) :: out(3,size(x,2))
      type(test_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=vartest_chisq(x(:,j),var0);out(:,j)=[r%statistic,r%df1,r%pvalue]
      end do
   end function column_vartests

   function column_aucs(x,label) result(auc)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: label(:)
      real(dp) :: auc(size(x,2)),r(size(x,1))
      integer :: j,n1,n0
      n1=count(label/=0);n0=size(label)-n1
      do j=1,size(x,2)
         r=rank_average(x(:,j))
         auc(j)=(sum(r,mask=label/=0)-real(n1*(n1+1),dp)/2.0_dp)/real(n1*n0,dp)
      end do
   end function column_aucs

   function one_way_anova(y,group) result(res)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: group(:)
      type(test_result) :: res
      integer, allocatable :: ug(:)
      integer :: k,g,ng,n
      real(dp) :: gm,ssb,ssw,m
      ug=unique_groups(group);n=size(y);gm=mean_r(y);ssb=0.0_dp;ssw=0.0_dp
      do k=1,size(ug)
         g=ug(k);ng=count(group==g);m=sum(y,mask=group==g)/real(ng,dp)
         ssb=ssb+real(ng,dp)*(m-gm)**2;ssw=ssw+sum((y-m)**2,mask=group==g)
      end do
      res%df1=real(size(ug)-1,dp);res%df2=real(n-size(ug),dp)
      res%statistic=(ssb/res%df1)/(ssw/res%df2);res%pvalue=1.0_dp-f_cdf(res%statistic,res%df1,res%df2)
   end function one_way_anova

   function one_way_anovas(x,group) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:)
      real(dp) :: out(4,size(x,2))
      type(test_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=one_way_anova(x(:,j),group);out(:,j)=[r%statistic,r%df1,r%df2,r%pvalue]
      end do
   end function one_way_anovas

   function correlation_pairs(x) result(r)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: r(:)
      integer :: p,i,j,k
      real(dp) :: xi(size(x,1)),xj(size(x,1)),mi,mj,si,sj
      p=size(x,2);allocate(r(p*(p-1)/2));k=0
      do j=1,p-1
         xj=x(:,j);mj=mean_r(xj);sj=sqrt(sum((xj-mj)**2))
         do i=j+1,p
            xi=x(:,i);mi=mean_r(xi);si=sqrt(sum((xi-mi)**2));k=k+1
            if(si>0.and.sj>0)then;r(k)=sum((xi-mi)*(xj-mj))/(si*sj);else;r(k)=0;end if
         end do
      end do
   end function correlation_pairs

   function column_correlations(x,y) result(r)
      real(dp), intent(in) :: x(:,:),y(:)
      real(dp) :: r(size(x,2));integer::j;real(dp)::my,sy,mx,sx
      my=mean_r(y);sy=sqrt(sum((y-my)**2))
      do j=1,size(x,2)
         mx=mean_r(x(:,j));sx=sqrt(sum((x(:,j)-mx)**2))
         if(sx>0.and.sy>0)then;r(j)=sum((x(:,j)-mx)*(y-my))/(sx*sy);else;r(j)=0;end if
      end do
   end function column_correlations

   function unique_groups(g) result(u)
      integer,intent(in)::g(:);integer,allocatable::u(:),tmp(:);integer::i,j,n,key
      tmp=g
      do i=2,size(tmp)
         key=tmp(i);j=i-1
         do while(j>=1)
            if(tmp(j)<=key)exit
            tmp(j+1)=tmp(j);j=j-1
         end do
         tmp(j+1)=key
      end do
      n=0
      do i=1,size(tmp)
         if(i==1)then;n=n+1
         else if(tmp(i)/=tmp(i-1))then;n=n+1
         end if
      end do
      allocate(u(n));n=0
      do i=1,size(tmp)
         if(i==1)then;n=n+1;u(n)=tmp(i)
         else if(tmp(i)/=tmp(i-1))then;n=n+1;u(n)=tmp(i)
         end if
      end do
   end function unique_groups

end module rfast_batch
