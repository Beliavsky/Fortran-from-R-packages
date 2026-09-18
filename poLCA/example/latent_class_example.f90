program latent_class_example
   use polca, only : dp, polca_model, polca_fit
   implicit none
   integer, parameter :: n = 800
   integer :: y(n, 3), n_choices(3), i, j
   real(dp) :: start(2, 2, 3), u
   type(polca_model) :: fit

   n_choices = 2
   start(1, :, 1) = [0.85_dp, 0.15_dp]
   start(1, :, 2) = [0.80_dp, 0.20_dp]
   start(1, :, 3) = [0.90_dp, 0.10_dp]
   start(2, :, 1) = [0.20_dp, 0.80_dp]
   start(2, :, 2) = [0.25_dp, 0.75_dp]
   start(2, :, 3) = [0.15_dp, 0.85_dp]

   call seed_example(20260913)
   do i = 1, n
      call random_number(u)
      if (u < 0.45_dp) then
         do j = 1, 3
            call random_number(u)
            y(i, j) = merge(1, 2, u < start(1, 1, j))
         end do
      else
         do j = 1, 3
            call random_number(u)
            y(i, j) = merge(1, 2, u < start(2, 1, j))
         end do
      end if
   end do

   call polca_fit(y, n_choices, 2, fit, probs_start=start, maxiter=1000, tol=1.0e-10_dp)
   write(*, '(a,f12.4)') 'log likelihood: ', fit%loglik
   write(*, '(a,2f10.5)') 'class shares:   ', fit%class_share
   write(*, '(a,3f10.5)') 'class 1 P(Y=2):', fit%probs(1, 2, :)
   write(*, '(a,3f10.5)') 'class 2 P(Y=2):', fit%probs(2, 2, :)

contains

   subroutine seed_example(seed)
      integer, intent(in) :: seed !! Scalar seed used to make the example reproducible on a given compiler runtime.
      integer, allocatable :: put(:)
      integer :: i, nseed

      call random_seed(size=nseed)
      allocate(put(nseed))
      do i = 1, nseed
         put(i) = modulo(seed + 8191 * i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_example

end program latent_class_example
