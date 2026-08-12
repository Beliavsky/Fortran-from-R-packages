program test_options
   use rdsdp
   implicit none
   type(dsdp_control) :: c
   integer :: u
   c=dsdp_control()
   open(newunit=u,file='test.options.tmp',status='replace')
   write(u,'(a)') '-gaptol 1e-9'
   write(u,'(a)') '-maxit 37'
   write(u,'(a)') '-print 2'
   write(u,'(a)') '-penalty 2e7'
   write(u,'(a)') '-usecg .true.'
   write(u,'(a)') '-cgtol 1e-8'
   write(u,'(a)') '-cgmaxit 55'
   write(u,'(a)') '-sparsedensity 0.15'
   write(u,'(a)') '-cgmatrixfree .false.'
   write(u,'(a)') '-usesparsefactor .true.'
   write(u,'(a)') '-sparsefactorthreshold 77'
   write(u,'(a)') '-sparsefactordensity 0.12'
   close(u)
   call read_dsdp_options('test.options.tmp',c)
   if (abs(c%gaptol-1.0e-9_dp)>tiny(1.0_dp)) error stop 'gaptol'
   if (c%maxiter/=37) error stop 'maxit'
   if (c%print/=2) error stop 'print'
   if (abs(c%penalty-2.0e7_dp)>1.0e-6_dp) error stop 'penalty'
   if (.not.c%use_cg) error stop 'usecg'
   if (abs(c%cg_tol-1.0e-8_dp)>tiny(1.0_dp)) error stop 'cgtol'
   if (c%cg_maxiter/=55) error stop 'cgmaxit'
   if (abs(c%sparse_density_threshold-0.15_dp)>tiny(1.0_dp)) error stop 'sparsedensity'
   if (c%cg_matrix_free) error stop 'cgmatrixfree'
   if (.not.c%use_sparse_schur_factor) error stop 'usesparsefactor'
   if (c%sparse_schur_threshold/=77) error stop 'sparsefactorthreshold'
   if (abs(c%sparse_schur_density_limit-0.12_dp)>tiny(1.0_dp)) error stop 'sparsefactordensity'
   open(newunit=u,file='test.options.tmp',status='old'); close(u,status='delete')
   print *, 'test_options: PASS'
end program test_options
