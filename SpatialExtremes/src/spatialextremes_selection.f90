module spatialextremes_selection
   use spatialextremes_base, only: dp
   implicit none
   private
   public :: tic_value, aic_value
contains
   pure real(dp) function aic_value(deviance,npar,k) result(v)
      real(dp),intent(in)::deviance
      integer,intent(in)::npar
      real(dp),intent(in),optional::k
      real(dp)::kk
      kk=2.0_dp
      if(present(k))kk=k
      v=deviance+kk*real(npar,dp)
   end function aic_value

   pure real(dp) function tic_value(deviance,ihessian,var_score,k) result(v)
      ! Numerical kernel of TIC.default: deviance + k tr(J H^{-1}),
      ! where upstream stores H^{-1} as ihessian and J as var.score.
      real(dp),intent(in)::deviance,ihessian(:,:),var_score(:,:)
      real(dp),intent(in),optional::k
      real(dp)::kk,pen(size(ihessian,1),size(ihessian,2))
      integer::i
      kk=2.0_dp
      if(present(k))kk=k
      pen=matmul(var_score,ihessian)
      v=deviance
      do i=1,min(size(pen,1),size(pen,2))
      v=v+kk*pen(i,i)
      end do
   end function tic_value
end module spatialextremes_selection
