program test_sdpa_io
   use rdsdp
   implicit none
   type(dsdp_problem) :: p,q
   type(dsdp_control) :: ctrl
   type(dsdp_solution) :: sol
   real(dp), allocatable :: xf(:),ap(:,:),aq(:,:)
   real(dp) :: err
   integer :: k,i,u

   call read_sdpa('data/truss1.dat-s',p)
   call write_sdpa('truss1-roundtrip.tmp',p)
   call read_sdpa('truss1-roundtrip.tmp',q)
   err=maxval(abs(p%b-q%b))
   do k=1,size(p%block)
      if (p%block(k)%category==dsdp_sdp_block) then
         call get_data_dense(p%block(k),0,ap); call get_data_dense(q%block(k),0,aq)
         err=max(err,maxval(abs(ap-aq)))
         do i=1,p%m
            call get_data_dense(p%block(k),i,ap); call get_data_dense(q%block(k),i,aq)
            err=max(err,maxval(abs(ap-aq)))
         end do
      else
         err=max(err,maxval(abs(p%block(k)%cdiag-q%block(k)%cdiag)))
         err=max(err,maxval(abs(p%block(k)%adiag-q%block(k)%adiag)))
      end if
   end do
   if (err>1.0e-13_dp) error stop 'SDPA roundtrip mismatch'
   ctrl=dsdp_control(); ctrl%gaptol=1.0e-6_dp; ctrl%pinfeastol=1.0e-6_dp; ctrl%rtol=1.0e-7_dp
   call dsdp_solve(p,sol,ctrl)
   call flatten_primal(sol,xf)
   if (.not.allocated(sol%y) .or. size(sol%y)/=p%m) error stop 'missing y'
   if (size(xf)/=sum([(p%block(k)%n*p%block(k)%n,k=1,size(p%block))])) then
      ! truss1 contains SDP blocks only
      error stop 'flattened X has wrong size'
   end if
   print '(a,es14.6,a,es10.2,a,es10.2)', 'truss1 dobj=',sol%dobj,' gap=',sol%relgap,' pinf=',sol%pinfeas
   open(newunit=u,file='truss1-roundtrip.tmp',status='old'); close(u,status='delete')
   print *, 'test_sdpa_io: PASS'
end program test_sdpa_io
