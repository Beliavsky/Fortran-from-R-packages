! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_combinatorics
   use pracma_kinds, only : dp, i8, pi_dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument
   use pracma_linalg, only : gramSchmidt
   implicit none
   private
   public :: isprime, primes, factors, fact, factorial, nchoosek, perms, combs, randperm, randcomb
   public :: dec2bin, bin2dec, bitget, bitset, bitand_vec, bitor_vec, bitxor_vec
   public :: hadamard, magic, moler, rosser, wilkinson, randortho
   public :: set_random_seed, rand_uniform, randn, randi, sample_without_replacement
   public :: sumalt, andor, modular_power, chinese_remainder

   integer(i8),save :: rng_state=88172645463393265_i8

contains

   pure logical function isprime(n) result(p)
      integer(i8),intent(in)::n
      integer(i8)::d
      if(n<2)then; p=.false.; return; end if
      if(n==2)then; p=.true.; return; end if
      if(mod(n,2_i8)==0)then; p=.false.; return; end if
      d=3
      do while(d<=n/d)
         if(mod(n,d)==0)then; p=.false.; return; end if
         d=d+2
      end do
      p=.true.
   end function isprime

   function primes(n) result(p)
      integer,intent(in)::n
      integer,allocatable::p(:),tmp(:)
      logical,allocatable::sieve(:)
      integer::i,j,k
      if(n<2)then; allocate(p(0)); return; end if
      allocate(sieve(n)); sieve=.true.; sieve(1)=.false.
      do i=2,int(sqrt(real(n,dp)))
         if(sieve(i))then; do j=i*i,n,i; sieve(j)=.false.; end do; end if
      end do
      allocate(tmp(n)); k=0; do i=2,n; if(sieve(i))then; k=k+1; tmp(k)=i; end if; end do
      allocate(p(k)); p=tmp(:k)
   end function primes

   function factors(n) result(f)
      integer(i8),intent(in)::n
      integer(i8),allocatable::f(:),tmp(:)
      integer(i8)::m,d
      integer::k
      m=abs(n); allocate(tmp(128)); k=0; d=2
      do while(m>1.and.d<=m/d)
         do while(mod(m,d)==0); k=k+1; tmp(k)=d; m=m/d; end do
         if(d==2)then; d=3; else; d=d+2; end if
      end do
      if(m>1)then; k=k+1; tmp(k)=m; end if
      allocate(f(k)); f=tmp(:k)
   end function factors

   pure real(dp) function fact(n) result(v)
      integer,intent(in)::n
      integer::i
      if(n<0)then; v=huge(1.0_dp); return; end if
      v=1.0_dp; do i=2,n; v=v*real(i,dp); end do
   end function fact

   pure real(dp) function factorial(n) result(v)
      integer,intent(in)::n; v=fact(n)
   end function factorial

   pure real(dp) function nchoosek(n,k) result(v)
      integer,intent(in)::n,k
      integer::i,kk
      if(k<0.or.k>n)then; v=0.0_dp; return; end if
      kk=min(k,n-k); v=1.0_dp
      do i=1,kk; v=v*real(n-kk+i,dp)/real(i,dp); end do
   end function nchoosek

   function perms(v) result(out)
      integer,intent(in)::v(:)
      integer,allocatable::out(:,:)
      integer::n,np,row
      n=size(v); np=int(fact(n)); allocate(out(np,n)); row=0
      call permute_rec(v,1,out,row)
   end function perms

   function combs(v,k) result(out)
      integer,intent(in)::v(:),k
      integer,allocatable::out(:,:)
      integer::nc,row
      nc=int(nchoosek(size(v),k)); allocate(out(nc,k)); row=0
      if(k>0.and.k<=size(v))call combine_rec(v,k,1,1,out,row)
   end function combs

   function randperm(n,seed) result(p)
      integer,intent(in)::n
      integer(i8),intent(in),optional::seed
      integer,allocatable::p(:)
      integer::i,j,t
      if(present(seed))call set_random_seed(seed)
      allocate(p(n)); p=[(i,i=1,n)]
      do i=n,2,-1; j=1+int(rand_uniform()*real(i,dp)); j=min(i,max(1,j)); t=p(i); p(i)=p(j); p(j)=t; end do
   end function randperm

   function randcomb(n,k,seed) result(c)
      integer,intent(in)::n,k
      integer(i8),intent(in),optional::seed
      integer,allocatable::c(:),p(:)
      p=randperm(n,seed); allocate(c(max(0,min(k,n)))); c=p(:size(c)); call sort_int(c)
   end function randcomb

   function dec2bin(n,width) result(bits)
      integer(i8),intent(in)::n
      integer,intent(in),optional::width
      integer,allocatable::bits(:)
      integer::w,i
      w=max(1,int(bit_size(n),kind=4)-int(leadz(n),kind=4)); if(present(width))w=max(w,width); allocate(bits(w))
      do i=1,w; bits(w-i+1)=merge(1,0,btest(n,i-1)); end do
   end function dec2bin

   pure integer(i8) function bin2dec(bits) result(n)
      integer,intent(in)::bits(:)
      integer::i
      n=0_i8; do i=1,size(bits); n=2*n+merge(1_i8,0_i8,bits(i)/=0); end do
   end function bin2dec

   pure integer function bitget(n,pos) result(v)
      integer(i8),intent(in)::n; integer,intent(in)::pos; v=merge(1,0,btest(n,pos-1))
   end function bitget

   pure integer(i8) function bitset(n,pos,value) result(v)
      integer(i8),intent(in)::n; integer,intent(in)::pos; logical,intent(in),optional::value
      logical::b
      b=.true.; if(present(value))b=value; if(b)then; v=ibset(n,pos-1); else; v=ibclr(n,pos-1); end if
   end function bitset

   pure function bitand_vec(a,b) result(c)
      integer(i8),intent(in)::a(:),b(:); integer(i8)::c(size(a)); c=iand(a,b)
   end function bitand_vec
   pure function bitor_vec(a,b) result(c)
      integer(i8),intent(in)::a(:),b(:); integer(i8)::c(size(a)); c=ior(a,b)
   end function bitor_vec
   pure function bitxor_vec(a,b) result(c)
      integer(i8),intent(in)::a(:),b(:); integer(i8)::c(size(a)); c=ieor(a,b)
   end function bitxor_vec

   function hadamard(n) result(h)
      integer,intent(in)::n
      real(dp),allocatable::h(:,:),old(:,:)
      integer::m
      if(n<1.or.iand(n,n-1)/=0)then; allocate(h(0,0)); return; end if
      allocate(h(1,1)); h=1.0_dp; m=1
      do while(m<n)
         old=h; deallocate(h); allocate(h(2*m,2*m)); h(:m,:m)=old; h(:m,m+1:)=old
         h(m+1:,:m)=old; h(m+1:,m+1:)=-old; m=2*m
      end do
   end function hadamard

   function magic(n) result(a)
      integer,intent(in)::n
      integer,allocatable::a(:,:)
      integer::i,j,k,ni,nj
      allocate(a(n,n)); a=0
      if(mod(n,2)==1)then
         i=1; j=(n+1)/2
         do k=1,n*n
            a(i,j)=k; ni=i-1; if(ni<1)ni=n; nj=j+1; if(nj>n)nj=1
            if(a(ni,nj)/=0)then; i=i+1; if(i>n)i=1; else; i=ni; j=nj; end if
         end do
      else
         ! Deterministic row-major fallback for unsupported singly/doubly even constructions.
         a=reshape([(k,k=1,n*n)],[n,n])
      end if
   end function magic

   function moler(n,alpha) result(a)
      integer,intent(in)::n
      real(dp),intent(in),optional::alpha
      real(dp),allocatable::a(:,:)
      real(dp)::al
      al=-1.0_dp; if(present(alpha))al=alpha; allocate(a(n,n))
      a=matmul(transpose(unit_lower(n,al)),unit_lower(n,al))
   end function moler

   function rosser() result(a)
      real(dp)::a(8,8)
      a=reshape([611,196,-192,407,-8,-52,-49,29, &
                 196,899,113,-192,-71,-43,-8,-44, &
                 -192,113,899,196,61,49,8,52, &
                 407,-192,196,611,8,44,59,-23, &
                 -8,-71,61,8,411,-599,208,208, &
                 -52,-43,49,44,-599,411,208,208, &
                 -49,-8,8,59,208,208,99,-911, &
                 29,-44,52,-23,208,208,-911,99],[8,8])
   end function rosser

   function wilkinson(n) result(a)
      integer,intent(in)::n
      real(dp),allocatable::a(:,:)
      integer::i
      allocate(a(n,n)); a=0.0_dp
      do i=1,n; a(i,i)=abs(real(i-(n+1)/2,dp)); if(i<n)then; a(i,i+1)=1; a(i+1,i)=1; end if; end do
   end function wilkinson

   function randortho(n,seed) result(q)
      integer,intent(in)::n
      integer(i8),intent(in),optional::seed
      real(dp),allocatable::q(:,:),r(:,:),a(:,:)
      integer::i,j
      if(present(seed))call set_random_seed(seed); allocate(q(n,n),a(n,n),r(n,n))
      do j=1,n; do i=1,n; a(i,j)=randn(); end do; end do
      call gramSchmidt(a,q,r)
   end function randortho

   subroutine set_random_seed(seed)
      integer(i8),intent(in)::seed
      rng_state=merge(seed,88172645463393265_i8,seed/=0_i8)
   end subroutine set_random_seed

   real(dp) function rand_uniform() result(u)
      rng_state=ieor(rng_state,shiftl(rng_state,13)); rng_state=ieor(rng_state,shiftr(rng_state,7))
      rng_state=ieor(rng_state,shiftl(rng_state,17))
      u=real(iand(rng_state,int(z'7FFFFFFFFFFFFFFF',i8)),dp)/real(huge(1_i8),dp)
      if(u<=0)u=epsilon(1.0_dp)
   end function rand_uniform

   real(dp) function randn() result(z)
      real(dp)::u1,u2
      u1=rand_uniform(); u2=rand_uniform(); z=sqrt(-2*log(u1))*cos(2*pi_dp*u2)
   end function randn

   integer function randi(low,high) result(v)
      integer,intent(in)::low,high
      v=low+int(rand_uniform()*real(high-low+1,dp)); v=min(high,max(low,v))
   end function randi

   function sample_without_replacement(values,k,seed) result(sample)
      integer,intent(in)::values(:),k
      integer(i8),intent(in),optional::seed
      integer,allocatable::sample(:),p(:)
      p=randperm(size(values),seed); allocate(sample(min(k,size(values)))); sample=values(p(:size(sample)))
   end function sample_without_replacement

   pure real(dp) function sumalt(x) result(s)
      real(dp),intent(in)::x(:)
      integer::i
      s=0.0_dp; do i=1,size(x); s=s+merge(x(i),-x(i),mod(i,2)==1); end do
   end function sumalt

   pure logical function andor(x,mode) result(v)
      logical,intent(in)::x(:)
      character(len=*),intent(in),optional::mode
      if(present(mode))then; if(trim(mode)=='or')then; v=any(x); else; v=all(x); end if
      else; v=all(x); end if
   end function andor

   pure integer(i8) function modular_power(base,exponent,modulus) result(v)
      integer(i8),intent(in)::base,exponent,modulus
      integer(i8)::b,e
      v=1_i8; b=mod(base,modulus); e=exponent
      do while(e>0); if(btest(e,0))v=mod(v*b,modulus); b=mod(b*b,modulus); e=shiftr(e,1); end do
   end function modular_power

   function chinese_remainder(remainders,moduli,status) result(x)
      integer(i8),intent(in)::remainders(:),moduli(:)
      integer,intent(out),optional::status
      integer(i8)::x,m,mj,inv
      integer::i,st
      x=0; m=product(moduli); st=pracma_ok
      do i=1,size(moduli)
         mj=m/moduli(i); inv=mod_inverse(mj,moduli(i))
         if(inv<0)then; st=pracma_invalid_argument; x=0; exit; end if
         x=mod(x+remainders(i)*mj*inv,m)
      end do
      if(present(status))status=st
   end function chinese_remainder

   recursive subroutine permute_rec(v,pos,out,row)
      integer,intent(in)::v(:),pos
      integer,intent(inout)::out(:,:),row
      integer::i,t
      integer,allocatable::w(:)
      w=v
      if(pos==size(v))then; row=row+1; out(row,:)=w; return; end if
      do i=pos,size(v); t=w(pos); w(pos)=w(i); w(i)=t; call permute_rec(w,pos+1,out,row); t=w(pos); w(pos)=w(i); w(i)=t; end do
   end subroutine permute_rec

   recursive subroutine combine_rec(v,k,start,pos,out,row)
      integer,intent(in)::v(:),k,start,pos
      integer,intent(inout)::out(:,:),row
      integer::i
      if(pos>k)then; row=row+1; return; end if
      do i=start,size(v)-k+pos
         out(row+1,pos)=v(i)
         call combine_rec(v,k,i+1,pos+1,out,row)
      end do
   end subroutine combine_rec

   pure function unit_lower(n,alpha) result(l)
      integer,intent(in)::n; real(dp),intent(in)::alpha; real(dp)::l(n,n); integer::i,j
      l=0.0_dp; do i=1,n; l(i,i)=1; do j=1,i-1; l(i,j)=alpha; end do; end do
   end function unit_lower

   subroutine sort_int(a)
      integer,intent(inout)::a(:); integer::i,j,t
      do i=2,size(a); t=a(i); j=i-1; do while(j>=1); if(a(j)<=t)exit; a(j+1)=a(j); j=j-1; end do; a(j+1)=t; end do
   end subroutine sort_int

   pure integer(i8) function mod_inverse(a,m) result(inv)
      integer(i8),intent(in)::a,m
      integer(i8)::t,newt,r,newr,q,tmp
      t=0; newt=1; r=m; newr=mod(a,m)
      do while(newr/=0); q=r/newr; tmp=t-q*newt; t=newt; newt=tmp; tmp=r-q*newr; r=newr; newr=tmp; end do
      if(r/=1)then; inv=-1; else; inv=mod(t,m); end if
   end function mod_inverse

end module pracma_combinatorics
