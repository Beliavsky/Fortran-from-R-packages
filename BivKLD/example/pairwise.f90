program pairwise
   use bivkld, only : BIVKLD_SUCCESS, biv_kld_matrix, biv_sample, dp
   implicit none

   type(biv_sample) :: samples(3)
   real(dp), allocatable :: divergence(:, :)
   real(dp) :: base(6, 2), H(2, 2, 3)
   integer :: i, info

   base(:, 1) = [-1.2_dp, -0.7_dp, -0.1_dp, 0.4_dp, 0.9_dp, 1.3_dp]
   base(:, 2) = [-0.4_dp, 0.5_dp, -0.8_dp, 0.9_dp, 0.1_dp, 1.1_dp]
   do i = 1, 3
      allocate(samples(i)%values(6, 2))
   end do
   samples(1)%values = base
   samples(2)%values = base
   samples(2)%values(:, 1) = samples(2)%values(:, 1) + 0.25_dp
   samples(3)%values = base
   samples(3)%values(:, 2) = samples(3)%values(:, 2) - 0.35_dp
   do i = 1, 3
      H(:, :, i) = reshape([0.35_dp, 0.0_dp, 0.0_dp, 0.45_dp], [2, 2])
   end do

   call biv_kld_matrix(samples, divergence, H=H, grid_size=[25, 25], info=info)
   if (info /= BIVKLD_SUCCESS) error stop 'biv_kld_matrix failed'
   do i = 1, size(divergence, 1)
      print '(*(f11.6,1x))', divergence(i, :)
   end do
end program pairwise
