! Computational RNG tests from randtoolbox R/testRNG.R and src/testrng.c.
module randtoolbox_tests
   use, intrinsic :: iso_fortran_env, only : real64, int64
   use randtoolbox_math, only : chi_square_survival
   implicit none
   private
   type, public :: rng_test_result
      real(real64) :: statistic=0.0_real64
      integer :: df=0
      real(real64) :: p_value=1.0_real64
      real(real64), allocatable :: observed(:), expected(:), residuals(:)
   end type rng_test_result
   public :: gap_test, frequency_test, serial_test, poker_test, order_test
   public :: collision_count, collision_test_counts, stirling_second, stirling_divided_by_k, permutations
contains
   function make_result(obs,exp,df) result(r)
      real(real64), intent(in) :: obs(:), exp(:)
      integer, intent(in) :: df
      type(rng_test_result) :: r
      integer :: i
      if(size(obs)/=size(exp)) error stop 'randtoolbox: test array size mismatch'
      allocate(r%observed(size(obs)),r%expected(size(obs)),r%residuals(size(obs)))
      r%observed=obs; r%expected=exp; r%df=df; r%statistic=0.0_real64
      do i=1,size(obs)
         if(exp(i)>0.0_real64) then
            r%residuals(i)=(obs(i)-exp(i))/sqrt(exp(i))
            r%statistic=r%statistic+r%residuals(i)**2
         else
            r%residuals(i)=0.0_real64
         end if
      end do
      r%p_value=chi_square_survival(r%statistic,df)
   end function make_result

   function frequency_test(u,bins,first_class) result(r)
      real(real64), intent(in) :: u(:)
      integer, intent(in), optional :: bins, first_class
      type(rng_test_result) :: r
      integer :: b,lo,i,k
      real(real64), allocatable :: obs(:),exp(:)
      b=16; if(present(bins)) b=bins; lo=0; if(present(first_class)) lo=first_class
      if(b<2) error stop 'randtoolbox: frequency bins must be >= 2'
      allocate(obs(b),exp(b)); obs=0.0_real64; exp=real(size(u),real64)/real(b,real64)
      do i=1,size(u)
         k=floor(u(i)*real(b,real64))+lo
         if(k>=lo .and. k<lo+b) obs(k-lo+1)=obs(k-lo+1)+1.0_real64
      end do
      r=make_result(obs,exp,b-1)
   end function frequency_test

   function serial_test(u,d) result(r)
      real(real64), intent(in) :: u(:)
      integer, intent(in), optional :: d
      type(rng_test_result) :: r
      integer :: dd,np,i,a,b,k
      real(real64), allocatable :: obs(:),exp(:)
      dd=8; if(present(d)) dd=d
      if(dd<2 .or. modulo(size(u),2)/=0) error stop 'randtoolbox: serial test requires even sample and d>=2'
      np=size(u)/2; allocate(obs(dd*dd),exp(dd*dd)); obs=0.0_real64
      exp=real(np,real64)/real(dd*dd,real64)
      ! R matrix(u,n/2,2): first half is first column, second half is second.
      do i=1,np
         a=floor(u(i)*dd); b=floor(u(i+np)*dd); k=a*dd+b+1
         if(k>=1 .and. k<=dd*dd) obs(k)=obs(k)+1.0_real64
      end do
      r=make_result(obs,exp,dd*dd-1)
   end function serial_test

   function poker_test(u,nbcard) result(r)
      real(real64), intent(in) :: u(:)
      integer, intent(in), optional :: nbcard
      type(rng_test_result) :: r
      integer :: d,nh,i,j,v,nd
      integer, allocatable :: counts(:)
      real(real64), allocatable :: obs(:),exp(:),st(:)
      real(real64) :: factd,falling
      d=5; if(present(nbcard)) d=nbcard
      if(d<2 .or. modulo(size(u),d)/=0) error stop 'randtoolbox: poker sample must be multiple of nbcard'
      nh=size(u)/d; allocate(obs(d),exp(d),counts(d),st(0:d)); obs=0.0_real64
      do i=1,nh
         counts=0
         do j=1,d
            v=floor(u(i+(j-1)*nh)*d)
            if(v>=0 .and. v<d) counts(v+1)=counts(v+1)+1
         end do
         nd=count(counts>0); if(nd>=1) obs(nd)=obs(nd)+1.0_real64
      end do
      st=stirling_second(d); factd=gamma(real(d+1,real64))
      do i=1,d
         falling=factd/gamma(real(d-i+1,real64))
         exp(i)=real(nh,real64)*falling*st(i)/real(d,real64)**d
      end do
      r=make_result(obs,exp,d-1)
   end function poker_test

   function order_test(u,d) result(r)
      real(real64), intent(in) :: u(:)
      integer, intent(in), optional :: d
      type(rng_test_result) :: r
      integer :: dd,nh,nperm,h,j,k,rank
      integer, allocatable :: p(:,:),ord(:)
      real(real64), allocatable :: obs(:),exp(:)
      dd=3; if(present(d)) dd=d
      if(dd<2 .or. dd>8 .or. modulo(size(u),dd)/=0) error stop 'randtoolbox: order test requires d=2..8 and compatible sample'
      p=permutations(dd); nperm=size(p,1); nh=size(u)/dd
      allocate(obs(nperm),exp(nperm),ord(dd)); obs=0.0_real64; exp=real(nh,real64)/real(nperm,real64)
      do h=1,nh
         do j=1,dd; ord(j)=j; end do
         ! Stable insertion sort of indices by the R column-major hand values.
         do j=2,dd
            rank=ord(j); k=j-1
            do while(k>=1)
               if(u(h+(ord(k)-1)*nh)<=u(h+(rank-1)*nh)) exit
               ord(k+1)=ord(k); k=k-1
            end do
            ord(k+1)=rank
         end do
         do k=1,nperm
            if(all(p(k,:)==ord)) then; obs(k)=obs(k)+1.0_real64; exit; end if
         end do
      end do
      r=make_result(obs,exp,nperm-1)
   end function order_test

   function gap_test(u,lower,upper) result(r)
      real(real64), intent(in) :: u(:)
      real(real64), intent(in), optional :: lower,upper
      type(rng_test_result) :: r
      real(real64) :: lo,hi,p,targ
      integer :: i,last,g,maxlen,ng
      integer, allocatable :: gaps(:),tmp(:)
      real(real64), allocatable :: obs(:),exp(:)
      lo=0.0_real64; hi=0.5_real64; if(present(lower)) lo=lower; if(present(upper)) hi=upper
      if(hi<=lo .or. lo<0.0_real64 .or. hi>1.0_real64) error stop 'randtoolbox: invalid gap interval'
      allocate(gaps(size(u)+1)); ng=0; last=0
      do i=1,size(u)
         if(u(i)<lo .or. u(i)>hi) then
            g=i-last-1; if(g>0) then; ng=ng+1; gaps(ng)=g; end if; last=i
         end if
      end do
      g=size(u)+1-last-1; if(g>0) then; ng=ng+1; gaps(ng)=g; end if
      p=hi-lo
      if(ng>0) then; maxlen=maxval(gaps(1:ng)); else; maxlen=1; end if
      targ=floor((log(0.1_real64)-2.0_real64*log(1.0_real64-p)-log(max(1.0_real64,real(size(u),real64))))/log(p))
      maxlen=max(maxlen,int(targ),1)
      allocate(obs(maxlen),exp(maxlen)); obs=0.0_real64
      do i=1,ng; if(gaps(i)<=maxlen) obs(gaps(i))=obs(gaps(i))+1.0_real64; end do
      do i=1,maxlen; exp(i)=(1.0_real64-p)**2*p**i*real(size(u),real64); end do
      r=make_result(obs,exp,maxlen-1)
      if(allocated(tmp)) deallocate(tmp)
   end function gap_test

   integer function collision_count(num,nb_urns) result(nc)
      integer, intent(in) :: num(:),nb_urns
      integer, allocatable :: urn(:)
      integer :: i
      if(nb_urns<1) error stop 'randtoolbox: number of urns must be positive'
      allocate(urn(nb_urns)); urn=0; nc=0
      do i=1,size(num)
         if(num(i)<0 .or. num(i)>=nb_urns) error stop 'randtoolbox: urn index out of range'
         if(urn(num(i)+1)/=0) nc=nc+1
         urn(num(i)+1)=urn(num(i)+1)+1
      end do
   end function collision_count

   function collision_test_counts(observed,len_sample,nb_cells) result(r)
      integer, intent(in) :: observed(:),len_sample
      integer(int64), intent(in) :: nb_cells
      type(rng_test_result) :: r
      integer :: cmin,cmax,c,i,k
      real(real64) :: lambda,prob,prodv
      real(real64), allocatable :: obs(:),expc(:),st(:)
      cmin=minval(observed); cmax=maxval(observed); allocate(obs(cmax-cmin+1),expc(cmax-cmin+1)); obs=0.0_real64; expc=0.0_real64
      do i=1,size(observed); obs(observed(i)-cmin+1)=obs(observed(i)-cmin+1)+1.0_real64; end do
      if(real(len_sample,real64)/real(nb_cells,real64)>1.0_real64/32.0_real64 .and. len_sample<=256) then
         st=stirling_second(len_sample)
         do c=cmin,cmax
            k=len_sample-c; prodv=1.0_real64
            do i=0,k-1; prodv=prodv*(1.0_real64-real(i,real64)/real(nb_cells,real64)); end do
            prob=prodv*st(k)/real(nb_cells,real64)**c
            expc(c-cmin+1)=prob*real(size(observed),real64)
         end do
      else
         lambda=real(len_sample,real64)**2/(2.0_real64*real(nb_cells,real64))
         do c=cmin,cmax
            expc(c-cmin+1)=exp(-lambda)*lambda**c/gamma(real(c+1,real64))*real(size(observed),real64)
         end do
      end if
      r=make_result(obs,expc,size(obs)-1)
   end function collision_test_counts

   function stirling_second(n) result(s)
      integer, intent(in) :: n
      real(real64), allocatable :: s(:)
      real(real64), allocatable :: old(:)
      integer :: i,k
      if(n<0) error stop 'randtoolbox: n must be nonnegative'
      allocate(s(0:n)); s=0.0_real64; s(0)=1.0_real64
      do i=1,n
         old=s; s=0.0_real64
         do k=1,i; s(k)=real(k,real64)*old(k)+old(k-1); end do
      end do
   end function stirling_second

   function stirling_divided_by_k(n,cmax,cste) result(s)
      integer, intent(in) :: n,cmax
      real(real64), intent(in) :: cste
      real(real64), allocatable :: s(:),old(:)
      integer :: i,k
      if(n<0 .or. cmax>=n .or. abs(cste)<=tiny(1.0_real64)) error stop 'randtoolbox: invalid Stirling arguments'
      allocate(s(0:n)); s=0.0_real64; s(0)=1.0_real64
      do i=1,n
         old=s; s=0.0_real64
         do k=1,i
            s(k)=real(k,real64)*old(k)+old(k-1)
            if(i>cmax+1) s(k)=real(k,real64)*old(k)+old(k-1)/cste
         end do
      end do
   end function stirling_divided_by_k

   function permutations(n) result(p)
      integer, intent(in) :: n
      integer, allocatable :: p(:,:)
      integer, allocatable :: old(:,:),q(:,:)
      integer :: m,r,pos,nold,rr
      if(n<1 .or. n>10) error stop 'randtoolbox: permutation size must be 1..10'
      allocate(p(1,1)); p=1
      do m=2,n
         old=p; nold=size(old,1); allocate(q(nold*m,m)); rr=0
         ! Matches upstream permut(): insert new symbol first, then each interior position, then last.
         do pos=1,m
            do r=1,nold
               rr=rr+1
               if(pos>1) q(rr,1:pos-1)=old(r,1:pos-1)
               q(rr,pos)=m
               if(pos<=m-1) q(rr,pos+1:m)=old(r,pos:m-1)
            end do
         end do
         call move_alloc(q,p)
      end do
   end function permutations
end module randtoolbox_tests
