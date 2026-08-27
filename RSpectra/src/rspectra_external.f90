! SPDX-License-Identifier: MPL-2.0
! Pinned pure-Fortran numerical backends for the RSpectra translation.
module rspectra_external
  use rfortran_arpack, only: dsaupd,dseupd,dnaupd,dneupd
  use la_lapack_d, only: dgetrf => la_dgetrf,dgetrs => la_dgetrs
  use la_lapack_d, only: dsyev => la_dsyev,dgeev => la_dgeev,dgesvd => la_dgesvd
  use la_lapack_z, only: zgetrf => la_zgetrf,zgetrs => la_zgetrs
  implicit none
  private
  public :: dsaupd,dseupd,dnaupd,dneupd
  public :: dgetrf,dgetrs,zgetrf,zgetrs,dsyev,dgeev,dgesvd
end module rspectra_external
