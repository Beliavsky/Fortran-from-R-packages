module spatialextremes_kriging
   use spatialextremes_base, only: dp,solve_spd
   use spatialextremes_covariance, only: covariance_values
   implicit none
   private
   public :: simple_kriging_weights,simple_kriging_predict
contains
   subroutine simple_kriging_weights(coord,coord_new,covmod,sill,range,smooth,nugget,weights,info,smooth2)
      real(dp),intent(in)::coord(:,:),coord_new(:,:),sill,range,smooth,nugget
      integer,intent(in)::covmod
      real(dp),intent(out)::weights(size(coord,1),size(coord_new,1))
      integer,intent(out)::info
      real(dp),intent(in),optional::smooth2
      real(dp)::c(size(coord,1),size(coord,1)),b(size(coord,1),size(coord_new,1)),h,s2
      integer::i,j
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      do i=1,size(coord,1)
         do j=1,size(coord,1)
            h=sqrt(sum((coord(i,:)-coord(j,:))**2))
            c(i,j)=covariance_values(h,covmod,nugget,sill,range,smooth,s2,size(coord,2))
         end do
         do j=1,size(coord_new,1)
            h=sqrt(sum((coord(i,:)-coord_new(j,:))**2))
            b(i,j)=covariance_values(h,covmod,0.0_dp,sill,range,smooth,s2,size(coord,2))
         end do
      end do
      call solve_spd(c,b,weights,info)
   end subroutine

   function simple_kriging_predict(obs,weights) result(pred)
      real(dp),intent(in)::obs(:),weights(:,:)
      real(dp)::pred(size(weights,2))
      pred=matmul(transpose(weights),obs)
   end function
end module spatialextremes_kriging
