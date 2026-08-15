module benford_tests
   use benford_kinds, only: dp
   use benford_math, only: normal_cdf, chi_square_sf, symmetric_eigen
   use benford_core, only: pbenf, qbenf, signifd_seq, benford_frequencies, benford_sample_frequencies
   implicit none
   private
   public :: benford_test_result_t, jointdigit_result_t
   public :: chisq_benftest, ks_benftest, mdist_benftest, edist_benftest
   public :: usq_benftest, meandigit_benftest, jpsq_benftest
   public :: jointdigit_benftest, jointdigit_benftest_indices, simulate_h0
   public :: stat_chisq, stat_ks, stat_mdist, stat_edist, stat_usq, stat_meandigit, stat_jpsq

   type :: benford_test_result_t
      real(dp) :: statistic=0.0_dp
      real(dp) :: p_value=1.0_dp
      integer :: n=0
      integer :: digits=1
   end type

   type :: jointdigit_result_t
      real(dp) :: statistic=0.0_dp
      real(dp) :: p_value=1.0_dp
      integer :: df=0
      integer :: n=0
      integer :: digits=1
      integer, allocatable :: eigenvalues_tested(:)
      real(dp), allocatable :: all_eigenvalues(:)
   end type

contains

   pure real(dp) function stat_chisq(f,p,n) result(s)
      real(dp), intent(in) :: f(:),p(:)
      integer, intent(in) :: n
      s=real(n,dp)*sum((f-p)**2/p)
   end function

   pure real(dp) function stat_ks(f,q,n) result(s)
      real(dp), intent(in) :: f(:),q(:)
      integer, intent(in) :: n
      real(dp) :: c(size(f))
      integer :: i
      c(1)=f(1)-q(1)
      do i=2,size(f)
         c(i)=c(i-1)+f(i)-f(i-1)-q(i)+q(i-1) ! overwritten below for clarity
      end do
      c=0.0_dp
      c(1)=f(1)-q(1)
      do i=2,size(f)
         c(i)=c(i-1)+f(i)-q(i)+q(i-1)
      end do
      s=sqrt(real(n,dp))*maxval(abs(c))
   end function

   pure real(dp) function stat_mdist(f,p,n) result(s)
      real(dp), intent(in) :: f(:),p(:)
      integer, intent(in) :: n
      s=sqrt(real(n,dp))*maxval(abs(f-p))
   end function

   pure real(dp) function stat_edist(f,p,n) result(s)
      real(dp), intent(in) :: f(:),p(:)
      integer, intent(in) :: n
      s=sqrt(real(n,dp)*sum((f-p)**2))
   end function

   pure real(dp) function stat_usq(f,p,n) result(s)
      real(dp), intent(in) :: f(:),p(:)
      integer, intent(in) :: n
      real(dp) :: c(size(f))
      integer :: i,m
      m=size(f); c(1)=f(1)-p(1)
      do i=2,m
         c(i)=c(i-1)+f(i)-p(i)
      end do
      s=real(n,dp)/real(m,dp)*(sum(c*c)-sum(c)**2/real(m,dp))
   end function

   pure real(dp) function stat_meandigit(f,p,digits) result(s)
      real(dp), intent(in) :: f(:),p(:)
      integer, intent(in) :: digits
      integer :: i,lo
      real(dp) :: mu,mu0,maxd
      lo=10**(digits-1); mu=0.0_dp; mu0=0.0_dp
      do i=1,size(f)
         mu=mu+f(i)*real(lo+i-1,dp)
         mu0=mu0+p(i)*real(lo+i-1,dp)
      end do
      maxd=real(lo+size(f)-1,dp)
      s=abs(mu-mu0)/(maxd-mu0)
   end function

   pure real(dp) function stat_jpsq(f,p) result(s)
      real(dp), intent(in) :: f(:),p(:)
      real(dp) :: fm,pm,num,den,r
      fm=sum(f)/real(size(f),dp); pm=sum(p)/real(size(p),dp)
      num=sum((f-fm)*(p-pm))
      den=sqrt(sum((f-fm)**2)*sum((p-pm)**2))
      if (den <= 0.0_dp) then
         s=0.0_dp
      else
         r=num/den
         s=sign(1.0_dp,r)*r*r
      end if
   end function

   subroutine observed(x,digits,f,p,q,n)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: digits
      real(dp), allocatable, intent(out) :: f(:),p(:),q(:)
      integer, intent(out) :: n
      integer, allocatable :: counts(:)
      call benford_frequencies(x,digits,counts,f,n)
      p=pbenf(digits); q=qbenf(digits)
   end subroutine

   subroutine chisq_benftest(x,digits,res,simulate,pvalsims)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: digits,pvalsims
      type(benford_test_result_t), intent(out) :: res
      logical, intent(in), optional :: simulate
      real(dp), allocatable :: f(:),p(:),q(:),h0(:)
      integer :: k,n,ns
      logical :: sim
      k=1;if(present(digits))k=digits; ns=10000;if(present(pvalsims))ns=pvalsims
      sim=.false.;if(present(simulate))sim=simulate
      call observed(x,k,f,p,q,n); res%statistic=stat_chisq(f,p,n)
      if(sim)then
         call simulate_h0('chisq',n,k,ns,h0); res%p_value=count(h0>res%statistic)/real(ns,dp)
      else
         res%p_value=chi_square_sf(res%statistic,size(p)-1)
      end if
      res%n=n;res%digits=k
   end subroutine

   subroutine ks_benftest(x,digits,res,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      call observed(x,k,f,p,q,n);res%statistic=stat_ks(f,q,n)
      call simulate_h0('ks',n,k,ns,h0);res%p_value=count(h0>res%statistic)/real(ns,dp);res%n=n;res%digits=k
   end subroutine

   subroutine mdist_benftest(x,digits,res,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      call observed(x,k,f,p,q,n);res%statistic=stat_mdist(f,p,n)
      call simulate_h0('mdist',n,k,ns,h0);res%p_value=count(h0>res%statistic)/real(ns,dp);res%n=n;res%digits=k
   end subroutine

   subroutine edist_benftest(x,digits,res,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      call observed(x,k,f,p,q,n);res%statistic=stat_edist(f,p,n)
      call simulate_h0('edist',n,k,ns,h0);res%p_value=count(h0>res%statistic)/real(ns,dp);res%n=n;res%digits=k
   end subroutine

   subroutine usq_benftest(x,digits,res,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      call observed(x,k,f,p,q,n);res%statistic=stat_usq(f,p,n)
      call simulate_h0('usq',n,k,ns,h0);res%p_value=count(h0>res%statistic)/real(ns,dp);res%n=n;res%digits=k
   end subroutine

   subroutine meandigit_benftest(x,digits,res,simulate,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res;logical,intent(in),optional::simulate
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns,i,lo
      real(dp)::mu0,var0,sd0,maxd;logical::sim
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      sim=.false.;if(present(simulate))sim=simulate
      call observed(x,k,f,p,q,n);res%statistic=stat_meandigit(f,p,k)
      if(sim)then
         call simulate_h0('meandigit',n,k,ns,h0);res%p_value=count(h0>res%statistic)/real(ns,dp)
      else
         lo=10**(k-1);mu0=0.0_dp
         do i=1,size(p);mu0=mu0+real(lo+i-1,dp)*p(i);end do
         var0=0.0_dp
         do i=1,size(p);var0=var0+(real(lo+i-1,dp)-mu0)**2*p(i);end do
         maxd=real(lo+size(p)-1,dp);sd0=sqrt(var0/real(n,dp))/(maxd-mu0)
         res%p_value=2.0_dp*(1.0_dp-normal_cdf(res%statistic/sd0))
      end if
      res%n=n;res%digits=k
   end subroutine

   subroutine jpsq_benftest(x,digits,res,pvalsims)
      real(dp),intent(in)::x(:);integer,intent(in),optional::digits,pvalsims
      type(benford_test_result_t),intent(out)::res
      real(dp),allocatable::f(:),p(:),q(:),h0(:);integer::k,n,ns
      k=1;if(present(digits))k=digits;ns=10000;if(present(pvalsims))ns=pvalsims
      call observed(x,k,f,p,q,n);res%statistic=stat_jpsq(f,p)
      call simulate_h0('jpsq',n,k,ns,h0);res%p_value=count(h0<=res%statistic)/real(ns,dp);res%n=n;res%digits=k
   end subroutine

   subroutine simulate_h0(teststatistic,n,digits,pvalsims,h0)
      character(len=*),intent(in)::teststatistic
      integer,intent(in)::n,digits,pvalsims
      real(dp),allocatable,intent(out)::h0(:)
      real(dp),allocatable::f(:),p(:),q(:)
      integer::i
      p=pbenf(digits);q=qbenf(digits);allocate(h0(pvalsims))
      do i=1,pvalsims
         call benford_sample_frequencies(n,digits,f)
         select case(trim(teststatistic))
         case('chisq');h0(i)=stat_chisq(f,p,n)
         case('ks');h0(i)=stat_ks(f,q,n)
         case('mdist');h0(i)=stat_mdist(f,p,n)
         case('edist');h0(i)=stat_edist(f,p,n)
         case('usq');h0(i)=stat_usq(f,p,n)
         case('meandigit');h0(i)=stat_meandigit(f,p,digits)
         case('jpsq');h0(i)=stat_jpsq(f,p)
         case default;h0(i)=0.0_dp
         end select
      end do
   end subroutine

   subroutine jointdigit_benftest(x,digits,res,mode,tol)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: digits
      type(jointdigit_result_t), intent(out) :: res
      character(len=*), intent(in), optional :: mode
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: f(:),p(:),q(:),cov(:,:),vals(:),vecs(:,:)
      real(dp), allocatable :: nzvals(:),nzvecs(:,:),pc(:),pc0(:)
      integer, allocatable :: idx(:)
      integer :: k,n,m,i,j,nnz,nkeep
      real(dp) :: tolerance,meanv
      character(len=16) :: md
      k=1
      if(present(digits)) k=digits
      tolerance=1.0d-15
      if(present(tol)) tolerance=tol
      md='all'
      if(present(mode)) md=adjustl(mode)
      call observed(x,k,f,p,q,n)
      m=size(p)
      allocate(cov(m,m))
      cov=-spread(p,2,m)*spread(p,1,m)
      do i=1,m
         cov(i,i)=p(i)*(1.0_dp-p(i))
      end do
      call symmetric_eigen(cov,vals,vecs)
      res%all_eigenvalues=vals
      nnz=count(abs(vals)>tolerance)
      allocate(nzvals(nnz),nzvecs(m,nnz))
      j=0
      do i=1,m
         if(abs(vals(i))>tolerance)then
            j=j+1
            nzvals(j)=vals(i)
            nzvecs(:,j)=vecs(:,i)
         end if
      end do
      if(trim(md)=='kaiser')then
         meanv=sum(nzvals)/real(nnz,dp)
         nkeep=count(nzvals>=meanv)
         allocate(idx(nkeep))
         j=0
         do i=1,nnz
            if(nzvals(i)>=meanv)then
               j=j+1
               idx(j)=i
            end if
         end do
      else
         nkeep=nnz
         allocate(idx(nkeep))
         idx=[(i,i=1,nkeep)]
      end if
      allocate(pc(nkeep),pc0(nkeep))
      pc=matmul(f,nzvecs(:,idx))
      pc0=matmul(p,nzvecs(:,idx))
      res%statistic=real(n,dp)*sum((pc-pc0)**2/nzvals(idx))
      res%df=nkeep
      res%p_value=chi_square_sf(res%statistic,nkeep)
      res%n=n
      res%digits=k
      res%eigenvalues_tested=idx
   end subroutine jointdigit_benftest

   subroutine jointdigit_benftest_indices(x,digits,indices,res,tol)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: digits,indices(:)
      type(jointdigit_result_t), intent(out) :: res
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: f(:),p(:),q(:),cov(:,:),vals(:),vecs(:,:)
      real(dp), allocatable :: nzvals(:),nzvecs(:,:),pc(:),pc0(:)
      integer, allocatable :: idx(:)
      integer :: n,m,i,j,nnz,nkeep
      real(dp) :: tolerance
      tolerance=1.0d-15
      if(present(tol)) tolerance=tol
      call observed(x,digits,f,p,q,n)
      m=size(p)
      allocate(cov(m,m))
      cov=-spread(p,2,m)*spread(p,1,m)
      do i=1,m
         cov(i,i)=p(i)*(1.0_dp-p(i))
      end do
      call symmetric_eigen(cov,vals,vecs)
      res%all_eigenvalues=vals
      nnz=count(abs(vals)>tolerance)
      allocate(nzvals(nnz),nzvecs(m,nnz))
      j=0
      do i=1,m
         if(abs(vals(i))>tolerance)then
            j=j+1
            nzvals(j)=vals(i)
            nzvecs(:,j)=vecs(:,i)
         end if
      end do
      nkeep=count(indices>=1 .and. indices<=nnz)
      allocate(idx(nkeep))
      j=0
      do i=1,size(indices)
         if(indices(i)>=1 .and. indices(i)<=nnz)then
            j=j+1
            idx(j)=indices(i)
         end if
      end do
      allocate(pc(nkeep),pc0(nkeep))
      pc=matmul(f,nzvecs(:,idx))
      pc0=matmul(p,nzvecs(:,idx))
      res%statistic=real(n,dp)*sum((pc-pc0)**2/nzvals(idx))
      res%df=nkeep
      res%p_value=chi_square_sf(res%statistic,nkeep)
      res%n=n
      res%digits=digits
      res%eigenvalues_tested=idx
   end subroutine jointdigit_benftest_indices
end module benford_tests
