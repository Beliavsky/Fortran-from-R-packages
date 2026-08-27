!
!     %---------------------------------------------%
!     | Initialize statistic and timing information |
!     | for nonsymmetric Arnoldi code.              |
!     %---------------------------------------------%
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2
!
subroutine dstatn
!
!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
!
include   'stat.h'
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
nopx   = 0
nbx    = 0
nrorth = 0
nitref = 0
nrstrt = 0
!
tnaupd = 0.0D+0
tnaup2 = 0.0D+0
tnaitr = 0.0D+0
tneigh = 0.0D+0
tngets = 0.0D+0
tnapps = 0.0D+0
tnconv = 0.0D+0
titref = 0.0D+0
tgetv0 = 0.0D+0
trvec  = 0.0D+0
!
!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
!
tmvopx = 0.0D+0
tmvbx  = 0.0D+0
!
return
!
!
!     %---------------%
!     | End of dstatn |
!     %---------------%
!
end
