module lavaan_tests
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general
   implicit none
   private
   public :: wald_test, linear_constraint_estimate
contains
   function wald_test(par,vcov,rmat,rhs,info) result(stat)
      real(dp),intent(in)::par(:),vcov(:,:),rmat(:,:),rhs(:)
      integer,intent(out)::info
      real(dp)::stat
      real(dp),allocatable::d(:),rv(:,:),ri(:,:)
      d=matmul(rmat,par)-rhs
      rv=matmul(rmat,matmul(vcov,transpose(rmat)))
      call inverse_general(rv,ri,info)
      if(info/=0) then
      stat=huge(1.0_dp)
      else
      stat=dot_product(d,matmul(ri,d))
      end if
   end function wald_test
   function linear_constraint_estimate(par,rmat) result(value)
      real(dp),intent(in)::par(:),rmat(:,:)
      real(dp),allocatable::value(:)
      value=matmul(rmat,par)
   end function linear_constraint_estimate
end module lavaan_tests
