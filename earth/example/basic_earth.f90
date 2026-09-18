program basic_earth
   use earth_api, only : dp, earth_model, earth_fit, earth_predict
   implicit none

   integer, parameter :: n = 61
   real(dp) :: x(n, 2)
   real(dp) :: y(n, 1)
   type(earth_model) :: model
   real(dp), allocatable :: prediction(:, :)
   integer :: i
   integer :: ierr
   character(len=:), allocatable :: message

   do i = 1, n
      x(i, 1) = real(i - 1, dp) / real(n - 1, dp)
      x(i, 2) = cos(0.41_dp * real(i, dp))
      y(i, 1) = 0.5_dp + 2.0_dp * max(0.0_dp, x(i, 1) - 0.30_dp)
   end do

   call earth_fit(x, y, model, degree=1, minspan=1, endspan=1, prune=.false., &
                  ierr=ierr, message=message)
   if (ierr /= 0) then
      write(*, '(a)') 'earth_fit failed: ' // trim(message)
      error stop 1
   end if
   call earth_predict(model, x, prediction, ierr)
   if (ierr /= 0) error stop 'earth_predict failed'

   write(*, '(a,i0)') 'forward terms:  ', model%n_forward_terms
   write(*, '(a,i0)') 'selected terms: ', model%n_selected
   write(*, '(a,f10.6)') 'R-squared:      ', model%rsq
   write(*, '(a,es12.4)') 'maximum error: ', maxval(abs(prediction - y))
end program basic_earth
