! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_types
   use r_compat, only: dp
   implicit none
   private
   public :: ph_type, dph_type, bivph_type, bivdph_type, mph_type, mdph_type, mphstar_type

   type :: ph_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s(:,:)
   end type ph_type

   type :: dph_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s(:,:)
   end type dph_type

   type :: bivph_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s11(:,:), s12(:,:), s22(:,:)
   end type bivph_type

   type :: bivdph_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s11(:,:), s12(:,:), s22(:,:)
   end type bivdph_type

   type :: mph_type
      real(dp), allocatable :: alpha(:)
      ! Common transient dimension p, d marginals: s(:,:,j).
      real(dp), allocatable :: s(:,:,:)
   end type mph_type

   type :: mdph_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s(:,:,:)
   end type mdph_type

   type :: mphstar_type
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: s(:,:)
      real(dp), allocatable :: reward(:,:)
   end type mphstar_type

end module matrixdist_types
