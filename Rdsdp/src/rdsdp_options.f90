! DSDP-style option-file reader used by the Rdsdp interface.
! DSDP copyright/license notice retained in LICENSE and licenses/DSDP-LICENSE.
module rdsdp_options
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_control
   implicit none
   private
   public :: read_dsdp_options
contains
   subroutine read_dsdp_options(filename,control)
      character(len=*), intent(in) :: filename
      type(dsdp_control), intent(inout) :: control
      integer :: u,ios,ival
      logical :: lval
      real(dp) :: rval
      character(len=512) :: line
      character(len=64) :: key,val
      open(newunit=u,file=filename,status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'read_dsdp_options: cannot open file'
      do
         read(u,'(A)',iostat=ios) line
         if (ios<0) exit
         if (ios>0) error stop 'read_dsdp_options: I/O failure'
         line=adjustl(line)
         if (len_trim(line)==0) cycle
         if (line(1:1)=='%' .or. line(1:1)=='#' .or. line(1:1)=='*') cycle
         key=''; val=''
         read(line,*,iostat=ios) key,val
         if (ios/=0) cycle
         select case(trim(adjustl(key)))
         case('-gaptol','gaptol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%gaptol=rval
         case('-maxit','maxit')
            read(val,*,iostat=ios) ival; if (ios==0 .and. ival>0) control%maxiter=ival
         case('-print','print')
            read(val,*,iostat=ios) ival; if (ios==0) control%print=ival
         case('-penalty','penalty')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%penalty=rval
         case('-pnormtol','pnormtol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%newton_tol=rval
         case('-infptol','infptol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%pinfeastol=rval
         case('-infdtol','infdtol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%rtol=rval
         case('-usesparse','usesparse')
            read(val,*,iostat=ios) lval; if (ios==0) control%use_sparse_data=lval
         case('-sparsedensity','sparsedensity')
            read(val,*,iostat=ios) rval
            if (ios==0 .and. rval>=0.0_dp .and. rval<=1.0_dp) control%sparse_density_threshold=rval
         case('-usecg','usecg')
            read(val,*,iostat=ios) lval; if (ios==0) control%use_cg=lval
         case('-cgtol','cgtol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>0.0_dp) control%cg_tol=rval
         case('-cgmaxit','cgmaxit')
            read(val,*,iostat=ios) ival; if (ios==0 .and. ival>0) control%cg_maxiter=ival
         case('-cgthreshold','cgthreshold')
            read(val,*,iostat=ios) ival; if (ios==0 .and. ival>=1) control%cg_threshold=ival
         case('-cgmatrixfree','cgmatrixfree')
            read(val,*,iostat=ios) lval; if (ios==0) control%cg_matrix_free=lval
         case('-usesparsefactor','usesparsefactor')
            read(val,*,iostat=ios) lval; if (ios==0) control%use_sparse_schur_factor=lval
         case('-sparsefactorthreshold','sparsefactorthreshold')
            read(val,*,iostat=ios) ival; if (ios==0 .and. ival>=1) control%sparse_schur_threshold=ival
         case('-sparsefactordensity','sparsefactordensity')
            read(val,*,iostat=ios) rval
            if (ios==0 .and. rval>=0.0_dp .and. rval<=1.0_dp) control%sparse_schur_density_limit=rval
         case('-sparsefactordroptol','sparsefactordroptol')
            read(val,*,iostat=ios) rval; if (ios==0 .and. rval>=0.0_dp) control%sparse_schur_drop_tol=rval
         case default
            ! Other historical DSDP tuning options that do not affect the
            ! translated numerical path are accepted by omission.
         end select
      end do
      close(u)
   end subroutine read_dsdp_options
end module rdsdp_options
