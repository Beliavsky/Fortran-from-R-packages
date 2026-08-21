! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_ops
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   implicit none
   private
   public :: csr_elem_add, csr_elem_subtract, csr_elem_multiply
   public :: csr_scale_values, csr_divide_scalar, csr_power_scalar
   public :: csr_apply_sqrt, csr_apply_abs, csr_apply_log1p, csr_apply_expm1
   public :: csr_apply_sin, csr_apply_sinh, csr_apply_tan, csr_apply_tanh
   public :: csr_apply_atanh, csr_apply_sign, csr_apply_floor, csr_apply_ceiling
   public :: csr_apply_tanpi
   public :: csr_apply_trunc, csr_apply_round, csr_apply_signif

contains

   subroutine csr_elem_add(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      call merge_rows(a,b,c,1.0_dp,1.0_dp,.false.,info)
   end subroutine csr_elem_add

   subroutine csr_elem_subtract(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      call merge_rows(a,b,c,1.0_dp,-1.0_dp,.false.,info)
   end subroutine csr_elem_subtract

   subroutine csr_elem_multiply(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      call merge_rows(a,b,c,1.0_dp,1.0_dp,.true.,info)
   end subroutine csr_elem_multiply

   subroutine merge_rows(a,b,c,sa,sb,multiply,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      real(dp), intent(in) :: sa,sb
      logical, intent(in) :: multiply
      integer, intent(out), optional :: info
      integer, allocatable :: rows(:),cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i,ka,kb,nz,cap,ca,cb,istat
      real(dp) :: v
      if (a%nrow/=b%nrow .or. a%ncol/=b%ncol) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      cap=max(1,a%nnz()+b%nnz())
      allocate(rows(cap),cols(cap),vals(cap)); nz=0
      do i=1,a%nrow
         ka=a%row_ptr(i); kb=b%row_ptr(i)
         do while (ka<a%row_ptr(i+1) .or. kb<b%row_ptr(i+1))
            if (ka>=a%row_ptr(i+1)) then
               if (multiply) exit
               cb=b%col_ind(kb); v=sb*b%values(kb); kb=kb+1
            else if (kb>=b%row_ptr(i+1)) then
               if (multiply) exit
               cb=a%col_ind(ka); v=sa*a%values(ka); ka=ka+1
            else
               ca=a%col_ind(ka); cb=b%col_ind(kb)
               if (ca==cb) then
                  if (multiply) then
                     v=a%values(ka)*b%values(kb)
                  else
                     v=sa*a%values(ka)+sb*b%values(kb)
                  end if
                  cb=ca; ka=ka+1; kb=kb+1
               else if (ca<cb) then
                  if (multiply) then
                     ka=ka+1; cycle
                  end if
                  cb=ca; v=sa*a%values(ka); ka=ka+1
               else
                  if (multiply) then
                     kb=kb+1; cycle
                  end if
                  v=sb*b%values(kb); kb=kb+1
               end if
            end if
            if (abs(v)>0.0_dp) then
               nz=nz+1; rows(nz)=i; cols(nz)=cb; vals(nz)=v
            end if
         end do
      end do
      call csr_from_triplet(a%nrow,a%ncol,rows(:nz),cols(:nz),vals(:nz),c,istat)
      if (present(info)) info=istat
   end subroutine merge_rows

   subroutine csr_scale_values(a,s)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: s
      if (allocated(a%values)) a%values=s*a%values
   end subroutine csr_scale_values

   subroutine csr_divide_scalar(a,s)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: s
      if (allocated(a%values)) a%values=a%values/s
   end subroutine csr_divide_scalar

   subroutine csr_power_scalar(a,p)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: p
      if (allocated(a%values)) a%values=a%values**p
   end subroutine csr_power_scalar

   subroutine csr_apply_sqrt(a); type(csr_matrix),intent(inout)::a; a%values=sqrt(a%values); end subroutine
   subroutine csr_apply_abs(a); type(csr_matrix),intent(inout)::a; a%values=abs(a%values); end subroutine
   subroutine csr_apply_log1p(a); type(csr_matrix),intent(inout)::a; a%values=log(1.0_dp+a%values); end subroutine
   subroutine csr_apply_expm1(a); type(csr_matrix),intent(inout)::a; a%values=exp(a%values)-1.0_dp; end subroutine
   subroutine csr_apply_sin(a); type(csr_matrix),intent(inout)::a; a%values=sin(a%values); end subroutine
   subroutine csr_apply_sinh(a); type(csr_matrix),intent(inout)::a; a%values=sinh(a%values); end subroutine
   subroutine csr_apply_tan(a); type(csr_matrix),intent(inout)::a; a%values=tan(a%values); end subroutine
   subroutine csr_apply_tanh(a); type(csr_matrix),intent(inout)::a; a%values=tanh(a%values); end subroutine
   subroutine csr_apply_tanpi(a); type(csr_matrix),intent(inout)::a; a%values=tan(acos(-1.0_dp)*a%values); end subroutine
   subroutine csr_apply_atanh(a); type(csr_matrix),intent(inout)::a; a%values=atanh(a%values); end subroutine
   subroutine csr_apply_sign(a)
      type(csr_matrix), intent(inout) :: a
      integer :: k
      do k=1,size(a%values)
         if (a%values(k)>0.0_dp) then
            a%values(k)=1.0_dp
         else if (a%values(k)<0.0_dp) then
            a%values(k)=-1.0_dp
         else
            a%values(k)=0.0_dp
         end if
      end do
   end subroutine csr_apply_sign
   subroutine csr_apply_floor(a); type(csr_matrix),intent(inout)::a; a%values=floor(a%values); end subroutine
   subroutine csr_apply_ceiling(a); type(csr_matrix),intent(inout)::a; a%values=ceiling(a%values); end subroutine
   subroutine csr_apply_trunc(a); type(csr_matrix),intent(inout)::a; a%values=aint(a%values); end subroutine

   subroutine csr_apply_round(a,digits)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in), optional :: digits
      integer :: d
      real(dp) :: s
      d=0; if (present(digits)) d=digits
      s=10.0_dp**real(d,dp)
      a%values=anint(a%values*s)/s
   end subroutine csr_apply_round

   subroutine csr_apply_signif(a,digits)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in), optional :: digits
      integer :: d,k,e
      real(dp) :: s,v
      d=6; if (present(digits)) d=digits
      do k=1,size(a%values)
         v=a%values(k)
         if (abs(v)<=tiny(1.0_dp)) cycle
         e=floor(log10(abs(v)))
         s=10.0_dp**real(d-1-e,dp)
         a%values(k)=anint(v*s)/s
      end do
   end subroutine csr_apply_signif

end module matrixextra_ops
