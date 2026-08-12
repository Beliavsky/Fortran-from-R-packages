module flsss_parallel
  implicit none
  private
  public :: openmp_enabled
contains
  logical function openmp_enabled() result(enabled)
    enabled = .false.
!$  enabled = .true.
  end function openmp_enabled
end module flsss_parallel
