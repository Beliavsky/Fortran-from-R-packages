module randtoolbox_quasi
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox_base, only : frac_part, two32, system_seed
   use randtoolbox_primes, only : get_primes
   use randtoolbox_math, only : normal_quantile
   use randtoolbox_sfmt, only : sfmt_rng, sfmt_supported
   use randtoolbox_sobol_data, only : sobol_initial, sobol_degree, sobol_a
   implicit none
   private
   integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64)
   public :: torus, halton, sobol, halton_radical_inverse, sobol_directions
contains
   function torus(n,dim,start,prime,normal,mixed,seed,mexp) result(x)
      integer, intent(in) :: n,dim
      integer(int64), intent(in), optional :: start,seed
      integer, intent(in), optional :: prime(:)
      logical, intent(in), optional :: normal,mixed
      integer, intent(in), optional :: mexp
      real(real64), allocatable :: x(:,:)
      integer, allocatable :: p(:)
      integer(int64) :: s0, raw
      integer :: i,j,me,d
      logical :: norm,mix
      type(sfmt_rng) :: g
      if(n<0 .or. dim<1 .or. dim>100000) error stop 'randtoolbox: invalid torus dimensions'
      d=dim
      if(present(prime)) then
         if(size(prime)<1) error stop 'randtoolbox: empty prime vector'
         d=size(prime); allocate(p(d)); p=prime
      else
         p=get_primes(d)
      end if
      if(any(p<=0)) error stop 'randtoolbox: torus bases must be positive'
      s0=1_int64; if(present(start))s0=start
      if(s0<0_int64) error stop 'randtoolbox: torus start must be nonnegative'
      norm=.false.; if(present(normal))norm=normal
      mix=.false.; if(present(mixed))mix=mixed
      me=19937; if(present(mexp))me=mexp
      allocate(x(n,d))
      if(mix) then
         if(.not.sfmt_supported(me)) error stop 'randtoolbox: invalid SFMT exponent'
         if(present(seed)) then; call g%init(me,seed,1); else; call g%init(me,system_seed(),1); end if
         do j=1,d
            do i=1,n
               raw=g%next_uint32()
               x(i,j)=frac_part(real(raw,real64)*sqrt(real(p(j),real64)))
            end do
         end do
      else
         do j=1,d
            do i=1,n
               x(i,j)=frac_part(real(s0+int(i-1,int64),real64)*sqrt(real(p(j),real64)))
            end do
         end do
      end if
      if(norm) x=normal_quantile(x)
   end function torus

   pure real(real64) function halton_radical_inverse(index,base) result(q)
      integer(int64), intent(in) :: index
      integer, intent(in) :: base
      integer(int64) :: iter,digit
      real(real64) :: f
      if(index<0_int64 .or. base<=1) error stop 'randtoolbox: invalid Halton arguments'
      iter=index; q=0.0_real64; f=1.0_real64/real(base,real64)
      do
         digit=modulo(iter,int(base,int64)); q=q+real(digit,real64)*f
         iter=(iter-digit)/int(base,int64); f=f/real(base,real64)
         if(iter==0_int64) exit
      end do
   end function halton_radical_inverse

   function halton(n,dim,start,normal,mixed,seed,mexp) result(x)
      integer, intent(in) :: n,dim
      integer(int64), intent(in), optional :: start,seed
      logical, intent(in), optional :: normal,mixed
      integer, intent(in), optional :: mexp
      real(real64), allocatable :: x(:,:)
      integer, allocatable :: p(:)
      integer(int64) :: s0,raw
      integer :: i,j,me
      logical :: norm,mix
      type(sfmt_rng)::g
      if(n<0 .or. dim<1 .or. dim>100000) error stop 'randtoolbox: invalid Halton dimensions'
      p=get_primes(dim); s0=1_int64; if(present(start))s0=start
      if(s0<0_int64) error stop 'randtoolbox: Halton start must be nonnegative'
      norm=.false.; if(present(normal))norm=normal
      mix=.false.; if(present(mixed))mix=mixed
      me=19937; if(present(mexp))me=mexp
      allocate(x(n,dim))
      if(mix) then
         if(present(seed)) then; call g%init(me,seed,1); else; call g%init(me,system_seed(),1); end if
         do j=1,dim
            do i=1,n
               raw=g%next_uint32(); x(i,j)=halton_radical_inverse(raw,p(j))
            end do
         end do
      else
         do j=1,dim
            do i=1,n
               x(i,j)=halton_radical_inverse(s0+int(i-1,int64),p(j))
            end do
         end do
      end if
      if(norm)x=normal_quantile(x)
   end function halton

   function sobol_directions(dim,maxbit) result(v)
      integer, intent(in) :: dim,maxbit
      integer(int64), allocatable :: v(:,:)
      integer :: i,j,k,s,a
      integer(int64)::t
      if(dim<1 .or. dim>1111) error stop 'randtoolbox: Sobol dimension must be 1..1111'
      if(maxbit<1 .or. maxbit>32) error stop 'randtoolbox: Sobol bit count must be 1..32'
      allocate(v(maxbit,dim)); v=0_int64
      do i=1,maxbit
         v(i,1)=shiftl(1_int64,32-i)
      end do
      do j=2,dim
         s=sobol_degree(j); a=sobol_a(j)
         do i=1,min(s,maxbit)
            v(i,j)=iand(shiftl(sobol_initial(i,j),32-i),mask32)
         end do
         do i=s+1,maxbit
            t=ieor(v(i-s,j),shiftr(v(i-s,j),s))
            do k=1,s-1
               if(btest(a,s-1-k)) t=ieor(t,v(i-k,j))
            end do
            v(i,j)=iand(t,mask32)
         end do
      end do
   end function sobol_directions

   function sobol(n,dim,start,normal) result(x)
      integer, intent(in) :: n,dim
      integer(int64), intent(in), optional :: start
      logical, intent(in), optional :: normal
      real(real64), allocatable :: x(:,:)
      integer(int64), allocatable :: v(:,:)
      integer(int64) :: s0,idx,g,state
      integer :: i,j,k,maxbit
      logical :: norm
      if(n<0 .or. dim<1 .or. dim>1111) error stop 'randtoolbox: invalid Sobol dimensions'
      s0=1_int64; if(present(start))s0=start
      if(s0<0_int64) error stop 'randtoolbox: Sobol start must be nonnegative'
      if(n==0) then; allocate(x(0,dim)); return; end if
      if(s0+int(n-1,int64)>int(z'FFFFFFFF',int64)) error stop 'randtoolbox: Sobol index exceeds 32-bit implementation'
      idx=s0+int(n-1,int64)
      maxbit=1
      do while(shiftr(idx,maxbit)>0_int64 .and. maxbit<32)
         maxbit=maxbit+1
      end do
      v=sobol_directions(dim,maxbit); allocate(x(n,dim))
      do i=1,n
         idx=s0+int(i-1,int64); g=ieor(idx,shiftr(idx,1))
         do j=1,dim
            state=0_int64
            do k=1,maxbit
               if(btest(g,k-1))state=ieor(state,v(k,j))
            end do
            x(i,j)=real(iand(state,mask32),real64)/two32
         end do
      end do
      norm=.false.; if(present(normal))norm=normal
      if(norm)x=normal_quantile(x)
   end function sobol
end module randtoolbox_quasi
