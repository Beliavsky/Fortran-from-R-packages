program test_triplets_psd
   use, intrinsic :: iso_c_binding, only : c_size_t
   use clarabel
   implicit none
   integer :: rows(5), cols(5)
   real(dp) :: values(5), dense_ref(3,2), sym(3,3)
   real(dp), allocatable :: dense(:,:), vector(:), recovered(:,:)
   type(csc_matrix) :: a, canceled

   rows = [1, 3, 1, 2, 1]
   cols = [1, 1, 1, 2, 2]
   values = [1.0_dp, 2.0_dp, 4.0_dp, 3.0_dp, -1.0_dp]
   a = csc_from_triplets(3, 2, rows, cols, values)
   if (a%nnz() /= 4) error stop "duplicate triplets were not aggregated"
   if (any(a%colptr /= [0_c_size_t, 2_c_size_t, 4_c_size_t])) error stop "bad triplet colptr"
   dense_ref = reshape([5.0_dp, 0.0_dp, 2.0_dp, -1.0_dp, 3.0_dp, 0.0_dp], shape(dense_ref))
   dense = a%to_dense()
   canceled = csc_from_triplets(1, 1, [1, 1], [1, 1], [1.0_dp, -1.0_dp])
   if (canceled%nnz() /= 0) error stop "cancelled duplicate remained"
   if (maxval(abs(dense - dense_ref)) > 1.0e-14_dp) error stop "triplet conversion"

   sym = reshape([2.0_dp, 1.0_dp, -3.0_dp, 1.0_dp, 4.0_dp, 0.5_dp, &
                  -3.0_dp, 0.5_dp, 5.0_dp], shape(sym))
   vector = psd_svec_upper(sym)
   recovered = psd_smat_upper(vector)
   if (size(vector) /= 6) error stop "PSD vector length"
   if (maxval(abs(recovered - sym)) > 1.0e-14_dp) error stop "PSD svec roundtrip"
   if (psd_matrix_order(6) /= 3 .or. psd_matrix_order(5) /= -1) error stop "PSD order"
   print *, "test_triplets_psd: PASS"
end program test_triplets_psd
