! SPDX-License-Identifier: GPL-3.0-or-later
program common_families
  use copula
  implicit none
  type(copula_model) :: models(5)
  real(dp) :: u(2)
  integer :: i

  u = [0.4_dp,0.6_dp]
  models(1) = clayton_copula(1.5_dp)
  models(2) = gumbel_copula(1.8_dp)
  models(3) = frank_copula(4.0_dp)
  models(4) = fgm_copula(0.7_dp)
  models(5) = plackett_copula(3.0_dp)

  print '(a)', ' family        C(u)        density         tau'
  do i = 1, size(models)
    print '(i4,3f14.7)', models(i)%family, pCopula(u,models(i)), dCopula(u,models(i)), tau(models(i))
  end do
end program common_families
