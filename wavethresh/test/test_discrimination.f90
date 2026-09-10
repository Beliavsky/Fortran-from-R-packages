program test_discrimination
   use wavethresh_types, only : dp, lda_model_t, wpst_discrimination_t, wpst_classification_t
   use wavethresh_utilities, only : bm_discr, makewpst_do, wpst_class
   implicit none

   real(dp) :: matrix(8, 2)
   real(dp) :: series(32)
   integer :: groups_small(8)
   integer :: groups(32)
   type(lda_model_t) :: lda
   type(wpst_discrimination_t) :: model
   type(wpst_classification_t) :: classification
   integer :: i

   matrix(1, :) = [0.0_dp, 0.0_dp]
   matrix(2, :) = [1.0_dp, 0.0_dp]
   matrix(3, :) = [0.0_dp, 1.0_dp]
   matrix(4, :) = [1.0_dp, 1.0_dp]
   matrix(5, :) = [4.0_dp, 4.0_dp]
   matrix(6, :) = [5.0_dp, 4.0_dp]
   matrix(7, :) = [4.0_dp, 5.0_dp]
   matrix(8, :) = [5.0_dp, 5.0_dp]
   groups_small = [1, 1, 1, 1, 2, 2, 2, 2]
   lda = bm_discr(matrix, groups_small)
   call require(lda%ok, "BMdiscr should fit a nonsingular two-class example")
   call require(all(lda%classes == [1, 2]), "BMdiscr class labels")
   call require(maxval(abs(lda%prior - 0.5_dp)) < 1.0e-12_dp, "BMdiscr empirical priors")
   call require(maxval(abs(lda%means(1, :) - [0.5_dp, 0.5_dp])) < 1.0e-12_dp, "BMdiscr class-one mean")
   call require(maxval(abs(lda%means(2, :) - [4.5_dp, 4.5_dp])) < 1.0e-12_dp, "BMdiscr class-two mean")
   call require(abs(lda%inverse_covariance(1, 1) - 3.0_dp) < 1.0e-11_dp, "BMdiscr inverse covariance 11")
   call require(abs(lda%inverse_covariance(2, 2) - 3.0_dp) < 1.0e-11_dp, "BMdiscr inverse covariance 22")
   call require(abs(lda%inverse_covariance(1, 2)) < 1.0e-11_dp, "BMdiscr inverse covariance 12")

   do i = 1, size(series)
      if (i <= size(series) / 2) then
         groups(i) = 1
         series(i) = 0.7_dp * (1.7_dp + sin(0.41_dp * real(i, dp)) + 0.13_dp * cos(1.2_dp * real(i, dp)))
      else
         groups(i) = 2
         series(i) = 3.0_dp * (1.7_dp + sin(0.41_dp * real(i, dp)) + 0.13_dp * cos(1.2_dp * real(i, dp)))
      end if
   end do
   model = makewpst_do(series, groups, filter_number=1.0_dp, family="DaubExPhase", mincor=0.4_dp)
   call require(model%ok, "makewpstDO should fit the deterministic packet-discrimination example")
   call require(size(model%level) >= 2, "makewpstDO should retain at least two basis features")
   call require(size(model%matrix, 1) == size(series), "makewpstDO training rows")
   call require(size(model%matrix, 2) == size(model%level), "makewpstDO feature metadata")

   classification = wpst_class(series, model)
   call require(classification%ok, "wpstCLASS should classify the deterministic training series")
   call require(size(classification%predicted_group) == size(groups), "wpstCLASS prediction count")
   call require(count(classification%predicted_group == groups) >= 28, "wpstCLASS should separate the amplitude regimes")
   call require(size(classification%discriminant_score, 2) == 2, "wpstCLASS should return one score per class")

   print *, "test_discrimination: PASS"

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Condition that must hold for the deterministic discrimination test.
      character(len=*), intent(in) :: message !! Failure explanation printed before terminating the test.

      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine require

end program test_discrimination
