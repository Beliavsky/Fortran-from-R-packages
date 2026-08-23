program test_directional_nb
   use rfast
   implicit none
   real(dp) :: theta(5), xg(6,2), xp(6,2), xm(4,3), xnew(2,2)
   integer :: y(6), cls(2), ym(4)
   type(circular_fit) :: cf
   type(gaussian_nb_model) :: gn
   type(poisson_nb_model) :: pn
   type(multinomial_nb_model) :: mn

   theta=[6.20_dp,0.02_dp,0.05_dp,6.25_dp,0.01_dp]
   cf=vm_mle(theta)
   call assert_true(cf%kappa>1.0_dp,'von Mises concentration')
   call assert_true(min(cf%mu,2.0_dp*pi-cf%mu)<0.10_dp,'von Mises direction')

   xg=reshape([0.0_dp,0.2_dp,-0.1_dp,4.9_dp,5.1_dp,5.0_dp, &
               0.1_dp,-0.1_dp,0.0_dp,5.2_dp,4.8_dp,5.0_dp],[6,2])
   y=[1,1,1,2,2,2]
   gn=fit_gaussian_nb(xg,y)
   xnew=reshape([0.05_dp,5.05_dp,0.0_dp,5.0_dp],[2,2])
   cls=predict_gaussian_nb(gn,xnew)
   call assert_true(all(cls==[1,2]),'gaussian naive Bayes')

   xp=reshape([0.0_dp,1.0_dp,0.0_dp,8.0_dp,9.0_dp,7.0_dp, &
               1.0_dp,0.0_dp,1.0_dp,7.0_dp,8.0_dp,9.0_dp],[6,2])
   pn=fit_poisson_nb(xp,y)
   cls=predict_poisson_nb(pn,reshape([0.0_dp,8.0_dp,1.0_dp,8.0_dp],[2,2]))
   call assert_true(all(cls==[1,2]),'Poisson naive Bayes')

   xm=reshape([8.0_dp,7.0_dp,1.0_dp,1.0_dp, &
               2.0_dp,3.0_dp,7.0_dp,8.0_dp, &
               0.0_dp,0.0_dp,2.0_dp,1.0_dp],[4,3])
   ym=[1,1,2,2]
   mn=fit_multinomial_nb(xm,ym)
   ! Use a three-column prediction to match the trained model.
   cls=predict_multinomial_nb(mn,reshape([9.0_dp,1.0_dp,1.0_dp,9.0_dp,0.0_dp,1.0_dp],[2,3]))
   call assert_true(all(cls==[1,2]),'multinomial naive Bayes')

   print *, 'test_directional_nb: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_directional_nb
