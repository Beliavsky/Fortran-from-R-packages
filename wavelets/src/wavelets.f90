! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets
   use wavelets_kinds, only : dp
   use wavelets_types, only : wt_filter_type, coefficient_level_type, wavelet_transform_type, mra_type
   use wavelets_filters, only : wt_filter, wt_filter_named, wt_filter_coefficients, wt_filter_qmf, &
                               wt_filter_equivalent, wt_filter_shift, waveletshift_dwt, scalingshift_dwt
   use wavelets_transform, only : dwt, dwt_named, dwt_with_filter, dwt_forward, dwt_backward, idwt, &
                                 modwt, modwt_named, modwt_with_filter, modwt_forward, modwt_backward, &
                                 imodwt, extend_series, align
   use wavelets_analysis, only : mra, mra_named, mra_with_filter
   implicit none
   private

   public :: dp
   public :: wt_filter_type
   public :: coefficient_level_type
   public :: wavelet_transform_type
   public :: mra_type
   public :: wt_filter
   public :: wt_filter_named
   public :: wt_filter_coefficients
   public :: wt_filter_qmf
   public :: wt_filter_equivalent
   public :: wt_filter_shift
   public :: waveletshift_dwt
   public :: scalingshift_dwt
   public :: dwt
   public :: dwt_named
   public :: dwt_with_filter
   public :: dwt_forward
   public :: dwt_backward
   public :: idwt
   public :: modwt
   public :: modwt_named
   public :: modwt_with_filter
   public :: modwt_forward
   public :: modwt_backward
   public :: imodwt
   public :: extend_series
   public :: align
   public :: mra
   public :: mra_named
   public :: mra_with_filter

end module wavelets
