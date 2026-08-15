module benford_core
   use benford_kinds, only: dp
   implicit none
   private
   public :: significant_digit, significant_digits, signifd_seq
   public :: pbenf, qbenf, rbenf, benford_frequencies, benford_sample_frequencies

contains

   pure integer function significant_digit(x, digits) result(d)
      real(dp), intent(in) :: x
      integer, intent(in), optional :: digits
      integer :: k, e
      real(dp) :: ax, scale
      k=1; if (present(digits)) k=digits
      if (k < 1 .or. x == 0.0_dp) then
         d=0
         return
      end if
      ax=abs(x)
      e=floor(log10(ax))
      scale=10.0_dp**real(k-1-e,dp)
      d=int(ax*scale)
      d=max(10**(k-1),min(10**k-1,d))
   end function significant_digit

   pure function significant_digits(x, digits) result(d)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: digits
      integer :: d(size(x)), i
      do i=1,size(x)
         d(i)=significant_digit(x(i),digits)
      end do
   end function significant_digits

   pure function signifd_seq(digits) result(seq)
      integer, intent(in), optional :: digits
      integer :: k, lo, hi, i
      integer, allocatable :: seq(:)
      k=1; if (present(digits)) k=digits
      lo=10**(k-1); hi=10**k-1
      allocate(seq(hi-lo+1))
      do i=1,size(seq)
         seq(i)=lo+i-1
      end do
   end function signifd_seq

   pure function pbenf(digits) result(p)
      integer, intent(in), optional :: digits
      integer :: k, lo, hi, i
      real(dp), allocatable :: p(:)
      k=1; if (present(digits)) k=digits
      lo=10**(k-1); hi=10**k-1
      allocate(p(hi-lo+1))
      do i=1,size(p)
         p(i)=log10(1.0_dp+1.0_dp/real(lo+i-1,dp))
      end do
   end function pbenf

   pure function qbenf(digits) result(q)
      integer, intent(in), optional :: digits
      real(dp), allocatable :: q(:), p(:)
      integer :: i
      p=pbenf(digits)
      allocate(q(size(p)))
      q(1)=p(1)
      do i=2,size(p)
         q(i)=q(i-1)+p(i)
      end do
      q(size(q))=1.0_dp
   end function qbenf

   subroutine rbenf(x)
      real(dp), intent(out) :: x(:)
      real(dp) :: u(size(x))
      call random_number(u)
      x=10.0_dp**u
   end subroutine rbenf

   subroutine benford_frequencies(x,digits,counts,relative,n_valid)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: digits
      integer, allocatable, intent(out) :: counts(:)
      real(dp), allocatable, intent(out) :: relative(:)
      integer, intent(out), optional :: n_valid
      integer :: lo, m, i, d, n
      lo=10**(digits-1); m=9*10**(digits-1)
      allocate(counts(m),relative(m)); counts=0; n=0
      do i=1,size(x)
         if (x(i) == 0.0_dp) cycle
         d=significant_digit(x(i),digits)
         if (d >= lo .and. d < lo+m) then
            counts(d-lo+1)=counts(d-lo+1)+1
            n=n+1
         end if
      end do
      if (n > 0) then
         relative=real(counts,dp)/real(n,dp)
      else
         relative=0.0_dp
      end if
      if (present(n_valid)) n_valid=n
   end subroutine benford_frequencies

   subroutine benford_sample_frequencies(n,digits,relative)
      integer, intent(in) :: n,digits
      real(dp), allocatable, intent(out) :: relative(:)
      real(dp), allocatable :: q(:)
      integer :: m, i, j
      real(dp) :: u
      m=9*10**(digits-1)
      q=qbenf(digits)
      allocate(relative(m)); relative=0.0_dp
      do i=1,n
         call random_number(u)
         do j=1,m
            if (u <= q(j)) then
               relative(j)=relative(j)+1.0_dp
               exit
            end if
         end do
      end do
      relative=relative/real(n,dp)
   end subroutine benford_sample_frequencies
end module benford_core
