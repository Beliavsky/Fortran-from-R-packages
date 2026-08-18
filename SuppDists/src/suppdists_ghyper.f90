module suppdists_ghyper
   use suppdists_kinds, only : dp
   use suppdists_stats, only : dist_stats
   implicit none
   private
   integer, parameter, public :: hyper_classic=1, hyper_iai=2, hyper_iaii=3, hyper_ib=4
   integer, parameter, public :: hyper_iia=5, hyper_iib=6, hyper_iiia=7, hyper_iiib=8, hyper_iv=9, hyper_none=0
   public :: hyper_type, hyper_type_name, dghyper, pghyper, qghyper, rghyper, sghyper, hyper_support
contains
   pure function hyper_type_name(a,k,nall) result(name)
      real(dp), intent(in) :: a,k,nall
      character(len=8) :: name
      select case(hyper_type(a,k,nall))
      case(hyper_classic); name='classic'
      case(hyper_iai); name='IAi'
      case(hyper_iaii); name='IAii'
      case(hyper_ib); name='IB'
      case(hyper_iia); name='IIA'
      case(hyper_iib); name='IIB'
      case(hyper_iiia); name='IIIA'
      case(hyper_iiib); name='IIIB'
      case(hyper_iv); name='IV'
      case default; name='no type'
      end select
   end function hyper_type_name

   pure logical function is_integer(x)
      real(dp),intent(in)::x
      is_integer=abs(x-anint(x))<1.0e-12_dp
   end function is_integer

   pure integer function hyper_type(a,k,nall) result(t)
      real(dp),intent(in)::a,k,nall
      if(a>0 .and. nall>0 .and. k>0 .and. is_integer(a) .and. is_integer(nall) .and. is_integer(k))then
         t=hyper_classic
      else if(a>0 .and. nall>0 .and. k>0 .and. is_integer(k) .and. k-1<a .and. a<nall-(k-1))then
         t=hyper_iai
      else if(a>0 .and. nall>0 .and. k>0 .and. is_integer(a) .and. a-1<k .and. k<nall-(a-1))then
         t=hyper_iaii
      else if(a>0 .and. nall>0 .and. k>0 .and. .not.is_integer(a) .and. .not.is_integer(k) .and. &
              a+k-1<nall .and. floor(a)==floor(k))then
         t=hyper_ib
      else if(a<0 .and. nall<k+a-1 .and. k>0 .and. is_integer(k))then
         t=hyper_iia
      else if(a<0 .and. nall>-1 .and. nall<k+a-1 .and. k>0 .and. .not.is_integer(k) .and. &
              floor(k)==floor(k+a-1-nall))then
         t=hyper_iib
      else if(a>0 .and. nall<k-1 .and. k<0 .and. is_integer(a))then
         t=hyper_iiia
      else if(a>0 .and. nall>-1 .and. nall<a+k-1 .and. k<0 .and. .not.is_integer(a) .and. &
              floor(a)==floor(a+k-1-nall))then
         t=hyper_iiib
      else if(a<0 .and. nall>-1 .and. k<0)then
         t=hyper_iv
      else
         t=hyper_none
      end if
   end function hyper_type

   pure subroutine hyper_support(a,k,nall,lo,hi,infinite,ok)
      real(dp),intent(in)::a,k,nall
      integer,intent(out)::lo,hi
      logical,intent(out)::infinite,ok
      integer::t
      t=hyper_type(a,k,nall);ok=t/=hyper_none;infinite=.false.;lo=0;hi=-1
      select case(t)
      case(hyper_classic);lo=max(0,nint(a+k-nall));hi=min(nint(a),nint(k))
      case(hyper_iai,hyper_iia);hi=nint(k)
      case(hyper_iaii,hyper_iiia);hi=nint(a)
      case(hyper_ib,hyper_iib,hyper_iiib,hyper_iv);infinite=.true.;hi=huge(1)
      end select
   end subroutine hyper_support

   pure real(dp) function dclassic(x,a,k,nall) result(p)
      integer,intent(in)::x
      integer,intent(in)::a,k,nall
      integer::lo,hi
      real(dp)::lp
      lo=max(0,a+k-nall);hi=min(a,k)
      if(x<lo .or. x>hi)then;p=0.0_dp;return;end if
      lp=log_gamma(real(k+1,dp))+log_gamma(real(nall-k+1,dp))+log_gamma(real(a+1,dp))+ &
         log_gamma(real(nall-a+1,dp))-log_gamma(real(x+1,dp))-log_gamma(real(k-x+1,dp))- &
         log_gamma(real(a-x+1,dp))-log_gamma(real(nall-a-k+x+1,dp))-log_gamma(real(nall+1,dp))
      p=exp(lp)
   end function dclassic

   pure subroutine gen_norm(a,k,nall,norm,maxx)
      real(dp),intent(in)::a,k,nall
      real(dp),intent(out)::norm
      integer,intent(out)::maxx
      real(dp)::term,sumv,r,b
      integer::i,hi
      logical::inf,ok
      call hyper_support(a,k,nall,i,hi,inf,ok)
      if(.not.ok)then;norm=0.0_dp;maxx=-1;return;end if
      b=nall-a;term=1.0_dp;sumv=1.0_dp;maxx=0
      if(.not.inf)then
         maxx=hi
         do i=0,hi-1
            r=((real(i,dp)-a)*(real(i,dp)-k))/ &
              (real(i+1,dp)*(b-k+real(i+1,dp)))
            term=term*r;sumv=sumv+term
         end do
      else
         do i=0,200000
            r=((real(i,dp)-a)*(real(i,dp)-k))/ &
              (real(i+1,dp)*(b-k+real(i+1,dp)))
            term=term*r
            if(term<0.0_dp .or. .not.(term<huge(1.0_dp)))exit
            sumv=sumv+term;maxx=i+1
            if(abs(term)<1.0e-15_dp*abs(sumv) .and. i>50)exit
         end do
      end if
      norm=sumv
   end subroutine gen_norm

   pure real(dp) function dghyper(x,a,k,nall) result(p)
      integer,intent(in)::x
      real(dp),intent(in)::a,k,nall
      integer::t,lo,hi,i,maxx
      logical::inf,ok
      real(dp)::norm,term,b,r
      t=hyper_type(a,k,nall)
      if(t==hyper_classic)then
         p=dclassic(x,nint(a),nint(k),nint(nall));return
      end if
      call hyper_support(a,k,nall,lo,hi,inf,ok)
      if(.not.ok .or. x<0 .or. (.not.inf .and. x>hi))then;p=0.0_dp;return;end if
      call gen_norm(a,k,nall,norm,maxx)
      if(norm<=0.0_dp .or. x>maxx .and. inf)then;p=0.0_dp;return;end if
      term=1.0_dp;b=nall-a
      do i=0,x-1
         r=((real(i,dp)-a)*(real(i,dp)-k))/(real(i+1,dp)*(b-k+real(i+1,dp)))
         term=term*r
      end do
      p=max(0.0_dp,term/norm)
   end function dghyper

   pure real(dp) function pghyper(x,a,k,nall) result(p)
      integer,intent(in)::x
      real(dp),intent(in)::a,k,nall
      integer::lo,hi,i
      logical::inf,ok
      call hyper_support(a,k,nall,lo,hi,inf,ok)
      if(.not.ok)then;p=0.0_dp;return;end if
      if(x<lo)then;p=0.0_dp;return;end if
      if(.not.inf .and. x>=hi)then;p=1.0_dp;return;end if
      p=0.0_dp
      do i=lo,x
         p=p+dghyper(i,a,k,nall)
      end do
      p=min(1.0_dp,p)
   end function pghyper

   pure integer function qghyper(p,a,k,nall) result(x)
      real(dp),intent(in)::p,a,k,nall
      integer::lo,hi
      logical::inf,ok
      call hyper_support(a,k,nall,lo,hi,inf,ok)
      if(.not.ok)then;x=0;return;end if
      x=lo
      do
         if(pghyper(x,a,k,nall)>=p)return
         x=x+1
         if(.not.inf .and. x>=hi)return
         if(x>200000)return
      end do
   end function qghyper

   integer function rghyper(a,k,nall) result(x)
      real(dp),intent(in)::a,k,nall
      real(dp)::u
      call random_number(u);x=qghyper(u,a,k,nall)
   end function rghyper

   pure function sghyper(a,k,nall) result(s)
      real(dp),intent(in)::a,k,nall
      type(dist_stats)::s
      integer::lo,hi,i,maxx
      logical::inf,ok
      real(dp)::p,d
      call hyper_support(a,k,nall,lo,hi,inf,ok)
      if(.not.ok)return
      if(inf)then;call gen_norm(a,k,nall,d,maxx);hi=maxx;end if
      s%mean=0.0_dp
      do i=lo,hi;s%mean=s%mean+real(i,dp)*dghyper(i,a,k,nall);end do
      do i=lo,hi
         p=dghyper(i,a,k,nall);d=real(i,dp)-s%mean
         s%variance=s%variance+p*d*d;s%third_central=s%third_central+p*d**3
         s%fourth_central=s%fourth_central+p*d**4
         if(p>dghyper(nint(s%mode),a,k,nall))s%mode=real(i,dp)
      end do
      s%median=real(qghyper(0.5_dp,a,k,nall),dp)
   end function sghyper
end module suppdists_ghyper
