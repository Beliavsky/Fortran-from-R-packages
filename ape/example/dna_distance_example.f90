! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program dna_distance_example
   use ape, only : dp, dna_a, dna_c, dna_g, dna_t, dna_distance_with_variance
   implicit none

   integer :: sequence_a(8)
   integer :: sequence_b(8)
   real(dp) :: distance
   real(dp) :: variance
   integer :: comparable_sites
   integer :: info

   sequence_a = [dna_a, dna_c, dna_g, dna_t, dna_a, dna_c, dna_g, dna_t]
   sequence_b = [dna_g, dna_c, dna_g, dna_t, dna_a, dna_c, dna_g, dna_t]

   call dna_distance_with_variance(sequence_a, sequence_b, 'K80', distance, variance, comparable_sites, info)
   if (info /= 0) error stop 'K80 distance/variance failed'

   print '(a,i0)', 'comparable sites: ', comparable_sites
   print '(a,f12.8)', 'K80 distance:     ', distance
   print '(a,f12.8)', 'K80 variance:     ', variance
end program dna_distance_example
