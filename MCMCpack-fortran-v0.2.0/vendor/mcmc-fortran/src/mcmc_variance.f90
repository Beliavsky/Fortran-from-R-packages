module mcmc_variance
   use mcmc_kinds, only : dp
   implicit none
   private

   type, public :: initseq_result
      real(dp) :: gamma0 = 0.0_dp
      real(dp), allocatable :: gamma_pos(:)
      real(dp), allocatable :: gamma_dec(:)
      real(dp), allocatable :: gamma_con(:)
      real(dp) :: var_pos = 0.0_dp
      real(dp) :: var_dec = 0.0_dp
      real(dp) :: var_con = 0.0_dp
      integer :: status = 0
   end type initseq_result

   public :: initseq, olbm

contains

   function initseq(x) result(res)
      real(dp), intent(in) :: x(:)
      type(initseq_result) :: res
      real(dp), allocatable :: xc(:),buff(:),puff(:)
      integer, allocatable :: nuff(:)
      real(dp) :: mu,gam1,gam2,muff
      integer :: n,i,j,lag1,lag2,nstep,jstep,k,ng

      n=size(x)
      if (n < 2) then
         res%status=1
         return
      end if
      mu=sum(x)/real(n,dp)
      xc=x-mu
      allocate(buff(n/2))
      res%gamma0=0.0_dp
      ng=0
      do i=0,n/2-1
         lag1=2*i
         gam1=0.0_dp
         do j=1,n-lag1
            gam1=gam1+xc(j)*xc(j+lag1)
         end do
         gam1=gam1/real(n,dp)
         if (i==0) res%gamma0=gam1
         lag2=lag1+1
         gam2=0.0_dp
         do j=1,n-lag2
            gam2=gam2+xc(j)*xc(j+lag2)
         end do
         gam2=gam2/real(n,dp)
         ng=ng+1
         buff(ng)=gam1+gam2
         if (buff(ng)<0.0_dp) then
            buff(ng)=0.0_dp
            exit
         end if
      end do

      allocate(res%gamma_pos(ng),res%gamma_dec(ng),res%gamma_con(ng))
      res%gamma_pos=buff(1:ng)
      do j=2,ng
         if (buff(j)>buff(j-1)) buff(j)=buff(j-1)
      end do
      res%gamma_dec=buff(1:ng)

      do j=ng,2,-1
         buff(j)=buff(j)-buff(j-1)
      end do
      allocate(puff(ng),nuff(ng))
      nstep=0
      do j=2,ng
         nstep=nstep+1
         puff(nstep)=buff(j)
         nuff(nstep)=1
         do while(nstep>1)
            if (puff(nstep)/real(nuff(nstep),dp) >= &
                puff(nstep-1)/real(nuff(nstep-1),dp)) exit
            puff(nstep-1)=puff(nstep-1)+puff(nstep)
            nuff(nstep-1)=nuff(nstep-1)+nuff(nstep)
            nstep=nstep-1
         end do
      end do
      j=2
      do jstep=1,nstep
         muff=puff(jstep)/real(nuff(jstep),dp)
         do k=1,nuff(jstep)
            buff(j)=buff(j-1)+muff
            j=j+1
         end do
      end do
      res%gamma_con=buff(1:ng)
      res%var_pos=2.0_dp*sum(res%gamma_pos)-res%gamma0
      res%var_dec=2.0_dp*sum(res%gamma_dec)-res%gamma0
      res%var_con=2.0_dp*sum(res%gamma_con)-res%gamma0
   end function initseq

   subroutine olbm(x,batch_length,var,mean,demean,status)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: batch_length
      real(dp), allocatable, intent(out) :: var(:,:)
      real(dp), allocatable, intent(out), optional :: mean(:)
      logical, intent(in), optional :: demean
      integer, intent(out), optional :: status
      integer :: n,p,i,j,k,l,nb
      real(dp), allocatable :: mu(:),work(:)
      logical :: dm

      if (present(status)) status=0
      n=size(x,1); p=size(x,2)
      if (batch_length<=0 .or. batch_length>n .or. p<=0) then
         allocate(var(0,0))
         if (present(mean)) allocate(mean(0))
         if (present(status)) status=1
         return
      end if
      dm=.true.; if (present(demean)) dm=demean
      allocate(mu(p),work(p),var(p,p))
      if (dm) then
         do i=1,p
            mu(i)=sum(x(:,i))/real(n,dp)
         end do
      else
         mu=0.0_dp
      end if
      if (present(mean)) then
         allocate(mean(p)); mean=mu
      end if
      mu=mu*real(batch_length,dp)
      var=0.0_dp
      work=0.0_dp
      do i=1,p
         work(i)=sum(x(1:batch_length,i))
      end do
      do i=1,p
         do j=1,i
            var(i,j)=(work(i)-mu(i))*(work(j)-mu(j))
         end do
      end do
      do k=1,n-batch_length
         l=k+batch_length
         work=work-x(k,:)+x(l,:)
         do i=1,p
            do j=1,i
               var(i,j)=var(i,j)+(work(i)-mu(i))*(work(j)-mu(j))
            end do
         end do
      end do
      nb=n-batch_length+1
      do i=1,p
         do j=1,i
            var(i,j)=var(i,j)/(real(nb*n*batch_length,dp))
            var(j,i)=var(i,j)
         end do
      end do
   end subroutine olbm

end module mcmc_variance
