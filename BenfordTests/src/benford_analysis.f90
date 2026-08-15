module benford_analysis
   use benford_kinds, only: dp
   use benford_math, only: normal_cdf, normal_quantile
   use benford_core, only: pbenf, benford_frequencies
   implicit none
   private
   public :: significant_digit_analysis_t, signifd_analysis

   type :: significant_digit_analysis_t
      integer :: n=0
      integer :: digits=1
      integer, allocatable :: digit(:)
      real(dp), allocatable :: frequency(:)
      real(dp), allocatable :: p_value(:)
      real(dp), allocatable :: quantile_probability(:)
      real(dp), allocatable :: confidence(:,:)
   end type
contains
   subroutine signifd_analysis(x,digits,result,alpha,relative)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits
      type(significant_digit_analysis_t),intent(out)::result
      real(dp),intent(in),optional::alpha(:);logical,intent(in),optional::relative
      integer,allocatable::counts(:);real(dp),allocatable::rel(:),p(:),a(:)
      real(dp)::ev,varv,zv,cdfv;integer::k,n,i,j,nq,lo;logical::relout
      k=1;if(present(digits))k=digits;relout=.true.;if(present(relative))relout=relative
      call benford_frequencies(x,k,counts,rel,n);p=pbenf(k);lo=10**(k-1)
      result%n=n;result%digits=k;allocate(result%digit(size(p)),result%frequency(size(p)),result%p_value(size(p)))
      do i=1,size(p);result%digit(i)=lo+i-1;end do
      if(relout)then;result%frequency=rel;else;result%frequency=real(counts,dp);end if
      do i=1,size(p)
         ev=real(n,dp)*p(i);varv=real(n,dp)*p(i)*(1.0_dp-p(i))
         if(varv>0.0_dp)then
            cdfv=normal_cdf((real(counts(i),dp)-ev)/sqrt(varv));result%p_value(i)=2.0_dp*min(cdfv,1.0_dp-cdfv)
         else;result%p_value(i)=1.0_dp;end if
      end do
      if(present(alpha))then;a=alpha;else;allocate(a(1));a=0.05_dp;end if
      nq=2*size(a)+1;allocate(result%quantile_probability(nq),result%confidence(nq,size(p)))
      do j=1,size(a);result%quantile_probability(j)=a(j)/2.0_dp;end do
      result%quantile_probability(size(a)+1)=0.5_dp
      do j=1,size(a);result%quantile_probability(size(a)+1+j)=1.0_dp-a(size(a)-j+1)/2.0_dp;end do
      do j=1,nq
         zv=normal_quantile(result%quantile_probability(j))
         do i=1,size(p)
            ev=real(n,dp)*p(i);varv=real(n,dp)*p(i)*(1.0_dp-p(i));result%confidence(j,i)=ev+zv*sqrt(varv)
            if(relout) result%confidence(j,i)=result%confidence(j,i)/real(n,dp)
         end do
      end do
   end subroutine
end module benford_analysis
