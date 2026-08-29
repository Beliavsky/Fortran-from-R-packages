program demo
use r_compat, only: dp
use compquadform_mod, only: davies_result_t, farebrother_result_t, &
   imhof_result_t, davies, farebrother, imhof, liu
implicit none
real(kind=dp) :: lambda(3), delta(3)
integer :: h(3)
type(davies_result_t) :: dres
type(farebrother_result_t) :: fres
type(imhof_result_t) :: ires
lambda = [0.5_dp, 1.2_dp, 2.0_dp]
h = [1, 2, 3]
delta = [0.0_dp, 0.7_dp, 0.3_dp]
dres = davies(6.0_dp, lambda, h, delta)
fres = farebrother(6.0_dp, lambda, h, delta)
ires = imhof(6.0_dp, lambda, h, delta)
write(*,'(a,g0)') 'Davies Qq:      ', dres%qq
write(*,'(a,g0)') 'Farebrother Qq: ', fres%qq
write(*,'(a,g0)') 'Imhof Qq:       ', ires%qq
write(*,'(a,g0)') 'Liu Qq:         ', liu(6.0_dp, lambda, h, delta)
end program demo
