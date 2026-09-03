! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! DNA computational routines are translated from ape R/DNA.R and
! src/dist_dna.c (Copyright 2005-2020 Emmanuel Paradis and contributors).
! R's DNAbin byte values are exposed through language-neutral integer state
! constants; internal helpers preserve the upstream ambiguity/bit semantics.
module ape_dna
   use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: dna_unknown = 0
   integer, parameter, public :: dna_a = 1
   integer, parameter, public :: dna_c = 2
   integer, parameter, public :: dna_g = 3
   integer, parameter, public :: dna_t = 4
   integer, parameter, public :: dna_gap = 5
   integer, parameter, public :: dna_r = 6
   integer, parameter, public :: dna_m = 7
   integer, parameter, public :: dna_w = 8
   integer, parameter, public :: dna_s = 9
   integer, parameter, public :: dna_k = 10
   integer, parameter, public :: dna_y = 11
   integer, parameter, public :: dna_v = 12
   integer, parameter, public :: dna_h = 13
   integer, parameter, public :: dna_d = 14
   integer, parameter, public :: dna_b = 15
   integer, parameter, public :: dna_n = 16

   public :: dna_distance
   public :: dna_distance_matrix
   public :: dna_distance_with_variance
   public :: dna_distance_matrix_with_variance
   public :: dna_bh87_matrix
   public :: dna_base_frequencies
   public :: dna_base_proportions
   public :: dna_leading_trailing_gaps_to_n
   public :: dna_global_deletion_mask
   public :: dna_segregating_sites
   public :: dna_contingency_table
   public :: dna_pattern_positions
   public :: translate_dna

contains

   pure subroutine dna_distance(sequence_a, sequence_b, model, distance, comparable_sites, info, &
      gamma_shape, base_frequency)
      !! Computes a pairwise nucleotide distance for ape `dist.dna` models using pairwise deletion.
      integer, intent(in) :: sequence_a(:) !! First sequence encoded A=1, C=2, G=3, T=4, gap=5; other values are unknown.
      integer, intent(in) :: sequence_b(:) !! Second sequence using the same integer encoding and site order.
      character(len=*), intent(in) :: model !! Pairwise model name; one of the non-matrix-only `dist.dna` models.
      real(dp), intent(out) :: distance !! Estimated distance; IEEE NaN when no comparable sites or a correction is undefined.
      integer, intent(out) :: comparable_sites !! Number of sites where both sequences contain unambiguous A/C/G/T states.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid lengths/model/correction domain.
      real(dp), intent(in), optional :: gamma_shape !! Positive gamma shape for JC69, K80, F81, or TN93 rate heterogeneity.
      real(dp), intent(in), optional :: base_frequency(4) !! A,C,G,T frequencies; pair-estimated when absent.
      real(dp) :: bf(4)

      if (present(base_frequency)) then
         bf = base_frequency
      else
         call pair_base_frequencies(sequence_a, sequence_b, bf)
      end if
      call dna_distance_core(sequence_a, sequence_b, model, bf, distance, comparable_sites, info, gamma_shape)
   end subroutine dna_distance

   pure subroutine dna_distance_matrix(sequences, model, distance, info, gamma_shape, base_frequency, pairwise_deletion)
      !! Computes an all-pairs DNA distance matrix with ape-style optional global or pairwise deletion.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site DNA codes: A=1, C=2, G=3, T=4, gap=5.
      character(len=*), intent(in) :: model !! One of the 17 ape `dist.dna` model names, case-insensitive.
      real(dp), allocatable, intent(out) :: distance(:, :) !! Symmetric taxon distance matrix with a zero diagonal.
      integer, intent(out) :: info !! Status code: zero if every pair succeeds, otherwise the first nonzero error code.
      real(dp), intent(in), optional :: gamma_shape !! Optional positive gamma shape for supported correction models.
      real(dp), intent(in), optional :: base_frequency(4) !! Optional global A,C,G,T frequencies used by frequency-dependent models.
      logical, intent(in), optional :: pairwise_deletion !! Use pairwise missing-data deletion; default is global deletion.
      integer, allocatable :: work(:, :)
      logical, allocatable :: keep(:)
      real(dp) :: bf(4)
      logical :: pair_delete
      character(len=16) :: chosen

      chosen = uppercase(trim(adjustl(model)))
      if (trim(chosen) == 'BH87') then
         call dna_bh87_matrix(sequences, distance, info, pairwise_deletion)
         return
      end if
      allocate(distance(size(sequences, 1), size(sequences, 1)))
      distance = 0.0_dp
      info = 0
      if (size(sequences, 1) == 0) return
      if (present(base_frequency)) then
         bf = base_frequency
      else
         call dna_base_frequencies(sequences, bf)
      end if
      pair_delete = .false.
      if (present(pairwise_deletion)) pair_delete = pairwise_deletion
      if (trim(chosen) == 'INDEL' .or. trim(chosen) == 'INDELBLOCK') pair_delete = .true.
      if (pair_delete) then
         allocate(work(size(sequences, 1), size(sequences, 2)))
         work = sequences
      else
         call dna_global_deletion_mask(sequences, keep)
         allocate(work(size(sequences, 1), count(keep)))
         work = pack_columns(sequences, keep)
      end if
      if (trim(chosen) == 'GG95') then
         call gg95_matrix(work, distance, info)
      else
         call pairwise_distance_matrix(work, chosen, bf, distance, info, gamma_shape)
      end if
   end subroutine dna_distance_matrix

   pure subroutine dna_distance_with_variance(sequence_a, sequence_b, model, distance, variance, comparable_sites, info, &
      gamma_shape, base_frequency)
      !! Computes one ape DNA distance together with its analytical sampling variance when upstream defines one.
      integer, intent(in) :: sequence_a(:) !! First sequence encoded with the package DNA integer convention.
      integer, intent(in) :: sequence_b(:) !! Second sequence with matching site order and length.
      character(len=*), intent(in) :: model !! Distance model name; variance is supported where upstream defines it.
      real(dp), intent(out) :: distance !! Pairwise distance estimate, matching `dna_distance`.
      real(dp), intent(out) :: variance !! Upstream analytical variance; IEEE NaN when unavailable or mathematically undefined.
      integer, intent(out) :: comparable_sites !! Number of sites with unambiguous bases in both sequences.
      integer, intent(out) :: info !! Zero on success; distance errors propagate and 6 means variance is unavailable.
      real(dp), intent(in), optional :: gamma_shape !! Positive gamma shape for JC69, K80, F81, or TN93.
      real(dp), intent(in), optional :: base_frequency(4) !! A,C,G,T frequencies; pair-estimated when absent.
      real(dp) :: bf(4)
      integer :: variance_info

      if (present(base_frequency)) then
         bf = base_frequency
      else
         call pair_base_frequencies(sequence_a, sequence_b, bf)
      end if
      if (present(gamma_shape)) then
         call dna_distance_core(sequence_a, sequence_b, model, bf, distance, comparable_sites, info, gamma_shape)
      else
         call dna_distance_core(sequence_a, sequence_b, model, bf, distance, comparable_sites, info)
      end if
      variance = ieee_value(0.0_dp, ieee_quiet_nan)
      if (info /= 0) return
      if (present(gamma_shape)) then
         call dna_variance_core(sequence_a, sequence_b, model, bf, variance, variance_info, gamma_shape)
      else
         call dna_variance_core(sequence_a, sequence_b, model, bf, variance, variance_info)
      end if
      if (variance_info /= 0) info = variance_info
   end subroutine dna_distance_with_variance

   pure subroutine dna_distance_matrix_with_variance(sequences, model, distance, variance, info, gamma_shape, &
      base_frequency, pairwise_deletion)
      !! Computes an ape-style DNA distance matrix and analytical sampling-variance matrix for models with translated formulas.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site DNA matrix in the package integer encoding.
      character(len=*), intent(in) :: model !! Distance model name; unsupported analytical variances return status 6.
      real(dp), allocatable, intent(out) :: distance(:, :) !! Symmetric distance matrix with zero diagonal.
      real(dp), allocatable, intent(out) :: variance(:, :) !! Symmetric analytical variance matrix; diagonal is zero.
      integer, intent(out) :: info !! Zero on success; model/distance errors propagate and 6 means unavailable variance.
      real(dp), intent(in), optional :: gamma_shape !! Optional positive gamma shape for supported gamma-corrected models.
      real(dp), intent(in), optional :: base_frequency(4) !! Optional global A,C,G,T frequencies used by frequency-dependent models.
      logical, intent(in), optional :: pairwise_deletion !! Use pairwise deletion; default is ape-style global deletion.
      integer, allocatable :: work(:, :)
      logical, allocatable :: keep(:)
      real(dp) :: bf(4)
      logical :: pair_delete
      character(len=16) :: chosen
      integer :: i
      integer :: j
      integer :: comparable
      integer :: pair_info
      real(dp) :: d
      real(dp) :: v

      chosen = uppercase(trim(adjustl(model)))
      allocate(distance(size(sequences, 1), size(sequences, 1)))
      allocate(variance(size(sequences, 1), size(sequences, 1)))
      distance = 0.0_dp
      variance = 0.0_dp
      info = 0
      if (size(sequences, 1) == 0) return
      if (.not. variance_model_supported(chosen)) then
         call dna_distance_matrix(sequences, model, distance, info, gamma_shape, base_frequency, pairwise_deletion)
         if (info == 0) info = 6
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         do i = 1, size(variance, 1)
            variance(i, i) = 0.0_dp
         end do
         return
      end if
      if (present(base_frequency)) then
         bf = base_frequency
      else
         call dna_base_frequencies(sequences, bf)
      end if
      pair_delete = .false.
      if (present(pairwise_deletion)) pair_delete = pairwise_deletion
      if (pair_delete) then
         allocate(work(size(sequences, 1), size(sequences, 2)))
         work = sequences
      else
         call dna_global_deletion_mask(sequences, keep)
         allocate(work(size(sequences, 1), count(keep)))
         work = pack_columns(sequences, keep)
      end if
      if (trim(chosen) == 'GG95') then
         call gg95_matrix(work, distance, info, variance)
         return
      end if
      do i = 1, size(work, 1) - 1
         do j = i + 1, size(work, 1)
            if (present(gamma_shape)) then
               call dna_distance_core(work(i, :), work(j, :), chosen, bf, d, comparable, pair_info, gamma_shape)
            else
               call dna_distance_core(work(i, :), work(j, :), chosen, bf, d, comparable, pair_info)
            end if
            distance(i, j) = d
            distance(j, i) = d
            if (pair_info == 0) then
               if (present(gamma_shape)) then
                  call dna_variance_core(work(i, :), work(j, :), chosen, bf, v, pair_info, gamma_shape)
               else
                  call dna_variance_core(work(i, :), work(j, :), chosen, bf, v, pair_info)
               end if
            else
               v = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            variance(i, j) = v
            variance(j, i) = v
            if (info == 0 .and. pair_info /= 0) info = pair_info
         end do
      end do
   end subroutine dna_distance_matrix_with_variance

   pure subroutine dna_bh87_matrix(sequences, distance, info, pairwise_deletion)
      !! Computes the directional Barry-Hartigan 1987 distance matrix used by ape model BH87.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site integer DNA matrix; only unambiguous A/C/G/T states contribute.
      real(dp), allocatable, intent(out) :: distance(:, :) !! Directional BH87 distances from row taxon to column taxon.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if a required transition determinant is nonpositive.
      logical, intent(in), optional :: pairwise_deletion !! If false (default), globally delete incomplete sites first.
      integer, allocatable :: work(:, :)
      logical, allocatable :: keep(:)
      real(dp) :: table(4, 4)
      real(dp) :: p(4, 4)
      real(dp) :: row_sum
      real(dp) :: determinant
      logical :: pair_delete
      logical :: valid_transition
      integer :: i
      integer :: j
      integer :: k

      pair_delete = .false.
      if (present(pairwise_deletion)) pair_delete = pairwise_deletion
      if (pair_delete) then
         allocate(work(size(sequences, 1), size(sequences, 2)))
         work = sequences
      else
         call dna_global_deletion_mask(sequences, keep)
         allocate(work(size(sequences, 1), count(keep)))
         work = pack_columns(sequences, keep)
      end if
      allocate(distance(size(work, 1), size(work, 1)))
      distance = 0.0_dp
      info = 0
      do i = 1, size(work, 1) - 1
         do j = i + 1, size(work, 1)
            call dna_contingency_table(work(i, :), work(j, :), table)
            p = table
            valid_transition = .true.
            do k = 1, 4
               row_sum = sum(p(k, :))
               if (row_sum <= 0.0_dp) then
                  valid_transition = .false.
               else
                  p(k, :) = p(k, :) / row_sum
               end if
            end do
            if (.not. valid_transition) then
               distance(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
               info = max(info, 1)
            else
               determinant = det4(p)
               if (determinant <= 0.0_dp) then
                  distance(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
                  info = max(info, 2)
               else
                  distance(i, j) = -0.25_dp * log(determinant)
               end if
            end if

            p = transpose(table)
            valid_transition = .true.
            do k = 1, 4
               row_sum = sum(p(k, :))
               if (row_sum <= 0.0_dp) then
                  valid_transition = .false.
               else
                  p(k, :) = p(k, :) / row_sum
               end if
            end do
            if (.not. valid_transition) then
               distance(j, i) = ieee_value(0.0_dp, ieee_quiet_nan)
               info = max(info, 1)
            else
               determinant = det4(p)
               if (determinant <= 0.0_dp) then
                  distance(j, i) = ieee_value(0.0_dp, ieee_quiet_nan)
                  info = max(info, 2)
               else
                  distance(j, i) = -0.25_dp * log(determinant)
               end if
            end if
         end do
      end do
   end subroutine dna_bh87_matrix

   pure subroutine dna_base_frequencies(sequences, frequencies)
      !! Estimates A,C,G,T frequencies from all unambiguous nucleotides in an alignment.
      integer, intent(in) :: sequences(:, :) !! Integer DNA matrix in which values 1 through 4 denote A,C,G,T.
      real(dp), intent(out) :: frequencies(4) !! Normalized A,C,G,T frequencies; all zeros if no known bases occur.
      integer :: i
      integer :: j
      integer :: total

      frequencies = 0.0_dp
      total = 0
      do i = 1, size(sequences, 1)
         do j = 1, size(sequences, 2)
            if (.not. known_base(sequences(i, j))) cycle
            frequencies(sequences(i, j)) = frequencies(sequences(i, j)) + 1.0_dp
            total = total + 1
         end do
      end do
      if (total > 0) frequencies = frequencies / real(total, dp)
   end subroutine dna_base_frequencies

   pure subroutine dna_base_proportions(sequences, proportions, frequencies)
      !! Counts or normalizes all 17 ape DNAbin states in A,C,G,T,R,M,W,S,K,Y,V,H,D,B,N,gap,unknown order.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site alignment using the public DNA integer-state constants.
      real(dp), intent(out) :: proportions(17) !! State counts when `frequencies` is true, otherwise proportions of all cells.
      logical, intent(in), optional :: frequencies !! If true, return raw counts as ape `base.freq(..., freq=TRUE, all=TRUE)`.
      integer :: i
      integer :: j
      integer :: slot
      logical :: raw_counts

      proportions = 0.0_dp
      do i = 1, size(sequences, 1)
         do j = 1, size(sequences, 2)
            slot = dna_state_slot(sequences(i, j))
            if (slot > 0) proportions(slot) = proportions(slot) + 1.0_dp
         end do
      end do
      raw_counts = .false.
      if (present(frequencies)) raw_counts = frequencies
      if (.not. raw_counts .and. size(sequences) > 0) proportions = proportions / real(size(sequences), dp)
   end subroutine dna_base_proportions

   pure subroutine dna_leading_trailing_gaps_to_n(sequences, converted)
      !! Replaces each sequence's leading and trailing gaps by N, matching ape `latag2n`/`leading_trailing_gaps_to_N`.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site DNA alignment whose internal gaps are preserved.
      integer, allocatable, intent(out) :: converted(:, :) !! Copy with only terminal runs of `dna_gap` replaced by `dna_n`.
      integer :: i
      integer :: j

      allocate(converted(size(sequences, 1), size(sequences, 2)))
      converted = sequences
      do i = 1, size(sequences, 1)
         j = 1
         do while (j <= size(sequences, 2))
            if (sequences(i, j) /= dna_gap) exit
            converted(i, j) = dna_n
            j = j + 1
         end do
         j = size(sequences, 2)
         do while (j >= 1)
            if (sequences(i, j) /= dna_gap) exit
            converted(i, j) = dna_n
            j = j - 1
         end do
      end do
   end subroutine dna_leading_trailing_gaps_to_n

   pure subroutine dna_global_deletion_mask(sequences, keep)
      !! Marks alignment columns in which every taxon has an unambiguous A/C/G/T base.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site integer DNA alignment.
      logical, allocatable, intent(out) :: keep(:) !! Site mask equivalent to ape's `GlobalDeletionDNA` known-base test.
      integer :: i
      integer :: j

      allocate(keep(size(sequences, 2)))
      keep = .true.
      do j = 1, size(sequences, 2)
         do i = 1, size(sequences, 1)
            if (.not. known_base(sequences(i, j))) then
               keep(j) = .false.
               exit
            end if
         end do
      end do
   end subroutine dna_global_deletion_mask

   pure subroutine dna_segregating_sites(sequences, segregating, strict, trailing_gaps_as_n)
      !! Identifies ape `seg.sites`, including DNAbin ambiguity compatibility and terminal-gap conversion.
      integer, intent(in) :: sequences(:, :) !! Taxa-by-site DNA alignment using unambiguous, IUPAC, gap, or unknown states.
      logical, allocatable, intent(out) :: segregating(:) !! Logical mask marking columns that ape regards as segregating.
      logical, intent(in), optional :: strict !! If true, any unequal state codes differ; default is ape's ambiguity-aware mode.
      logical, intent(in), optional :: trailing_gaps_as_n !! Replace terminal gaps by N before testing; default is true as in R.
      integer, allocatable :: work(:, :)
      integer :: base
      integer :: i
      integer :: j
      logical :: done
      logical :: do_strict
      logical :: convert_terminal_gaps

      convert_terminal_gaps = .true.
      if (present(trailing_gaps_as_n)) convert_terminal_gaps = trailing_gaps_as_n
      if (convert_terminal_gaps) then
         call dna_leading_trailing_gaps_to_n(sequences, work)
      else
         allocate(work(size(sequences, 1), size(sequences, 2)))
         work = sequences
      end if

      allocate(segregating(size(work, 2)))
      segregating = .false.
      do_strict = .false.
      if (present(strict)) do_strict = strict
      do j = 1, size(work, 2)
         if (size(work, 1) <= 1) cycle
         if (do_strict) then
            base = work(1, j)
            do i = 2, size(work, 1)
               if (work(i, j) /= base) then
                  segregating(j) = .true.
                  exit
               end if
            end do
            cycle
         end if

         i = 1
         base = work(i, j)
         done = .false.
         do while (.not. known_base(base))
            i = i + 1
            if (i > size(work, 1)) then
               done = .true.
               exit
            end if
            if (base /= work(i, j)) then
               if (base /= dna_unknown .and. work(i, j) /= dna_unknown) then
                  if (dna_raw_code(base) > 4) then
                     if (work(i, j) == dna_gap) then
                        segregating(j) = .true.
                        done = .true.
                        exit
                     else if (dna_states_surely_different(base, work(i, j))) then
                        segregating(j) = .true.
                        done = .true.
                        exit
                     end if
                  else
                     segregating(j) = .true.
                     done = .true.
                     exit
                  end if
               end if
               base = work(i, j)
            end if
         end do
         if (done) cycle

         i = i + 1
         do while (i <= size(work, 1))
            if (work(i, j) /= base) then
               if (work(i, j) == dna_gap) then
                  segregating(j) = .true.
                  exit
               end if
               if (dna_states_surely_different(work(i, j), base)) then
                  segregating(j) = .true.
                  exit
               end if
            end if
            i = i + 1
         end do
      end do
   end subroutine dna_segregating_sites

   pure subroutine dna_contingency_table(sequence_a, sequence_b, table)
      !! Counts paired A,C,G,T states in a 4-by-4 contingency table.
      integer, intent(in) :: sequence_a(:) !! First integer DNA sequence.
      integer, intent(in) :: sequence_b(:) !! Second sequence; sites beyond the shorter input are ignored.
      real(dp), intent(out) :: table(4, 4) !! Counts with rows for sequence A states and columns for sequence B states.
      integer :: i
      integer :: n

      table = 0.0_dp
      n = min(size(sequence_a), size(sequence_b))
      do i = 1, n
         if (.not. known_base(sequence_a(i)) .or. .not. known_base(sequence_b(i))) cycle
         table(sequence_a(i), sequence_b(i)) = table(sequence_a(i), sequence_b(i)) + 1.0_dp
      end do
   end subroutine dna_contingency_table

   pure subroutine dna_pattern_positions(sequence, pattern, positions)
      !! Returns one-based starting positions of exact occurrences of an integer DNA pattern.
      integer, intent(in) :: sequence(:) !! Integer DNA sequence searched without ambiguity expansion.
      integer, intent(in) :: pattern(:) !! Nonempty integer pattern matched exactly.
      integer, allocatable, intent(out) :: positions(:) !! One-based starting positions of all overlapping matches.
      integer, allocatable :: buffer(:)
      integer :: i
      integer :: matches

      if (size(pattern) == 0 .or. size(sequence) < size(pattern)) then
         allocate(positions(0))
         return
      end if
      allocate(buffer(size(sequence) - size(pattern) + 1))
      matches = 0
      do i = 1, size(sequence) - size(pattern) + 1
         if (all(sequence(i:i + size(pattern) - 1) == pattern)) then
            matches = matches + 1
            buffer(matches) = i
         end if
      end do
      allocate(positions(matches))
      positions = buffer(1:matches)
   end subroutine dna_pattern_positions

   pure subroutine translate_dna(sequence, amino_acids, info, reading_frame, genetic_code)
      !! Translates DNA with IUPAC ambiguity resolution for ape genetic codes 1 through 6.
      integer, intent(in) :: sequence(:) !! DNA sequence using the public A/C/G/T, IUPAC ambiguity, gap, and unknown codes.
      character(len=1), allocatable, intent(out) :: amino_acids(:) !! One-letter amino acids; `*` is stop and `X` unresolved.
      integer, intent(out) :: info !! Status: zero on success, 1 for invalid frame, or 2 for unsupported genetic code.
      integer, intent(in), optional :: reading_frame !! One-based starting frame 1, 2, or 3; default is 1.
      integer, intent(in), optional :: genetic_code !! Ape translation-table number 1 through 6; default is standard code 1.
      integer :: code
      integer :: frame
      integer :: i
      integer :: k
      integer :: ncodon

      frame = 1
      if (present(reading_frame)) frame = reading_frame
      if (frame < 1 .or. frame > 3) then
         allocate(amino_acids(0))
         info = 1
         return
      end if
      code = 1
      if (present(genetic_code)) code = genetic_code
      if (code < 1 .or. code > 6) then
         allocate(amino_acids(0))
         info = 2
         return
      end if
      ncodon = max(0, (size(sequence) - frame + 1) / 3)
      allocate(amino_acids(ncodon))
      info = 0
      k = 0
      do i = frame, frame + 3 * ncodon - 1, 3
         k = k + 1
         amino_acids(k) = codon_amino(sequence(i), sequence(i + 1), sequence(i + 2), code)
      end do
   end subroutine translate_dna

   pure subroutine pairwise_distance_matrix(sequences, model, bf, distance, info, gamma_shape)
      integer, intent(in) :: sequences(:, :) !! Working taxa-by-site DNA alignment after deletion policy is applied.
      character(len=*), intent(in) :: model !! Uppercase distance-model name.
      real(dp), intent(in) :: bf(4) !! Global A,C,G,T frequencies for frequency-dependent models.
      real(dp), intent(inout) :: distance(:, :) !! Preallocated symmetric matrix receiving pairwise distances.
      integer, intent(out) :: info !! Zero if all pairs succeed, otherwise first pairwise error code.
      real(dp), intent(in), optional :: gamma_shape !! Optional gamma shape forwarded to supported models.
      integer :: comparable
      integer :: i
      integer :: j
      integer :: pair_info
      real(dp) :: d

      info = 0
      do i = 1, size(sequences, 1) - 1
         do j = i + 1, size(sequences, 1)
            if (present(gamma_shape)) then
               call dna_distance_core(sequences(i, :), sequences(j, :), model, bf, d, comparable, pair_info, gamma_shape)
            else
               call dna_distance_core(sequences(i, :), sequences(j, :), model, bf, d, comparable, pair_info)
            end if
            distance(i, j) = d
            distance(j, i) = d
            if (info == 0 .and. pair_info /= 0) info = pair_info
         end do
      end do
   end subroutine pairwise_distance_matrix

   pure subroutine dna_distance_core(sequence_a, sequence_b, model, bf, distance, comparable_sites, info, gamma_shape)
      integer, intent(in) :: sequence_a(:) !! First pairwise DNA sequence.
      integer, intent(in) :: sequence_b(:) !! Second pairwise DNA sequence of equal length.
      character(len=*), intent(in) :: model !! Uppercase ape distance model name.
      real(dp), intent(in) :: bf(4) !! A,C,G,T base frequencies for frequency-dependent models.
      real(dp), intent(out) :: distance !! Computed distance or IEEE NaN when undefined.
      integer, intent(out) :: comparable_sites !! Number of sites with two unambiguous bases.
      integer, intent(out) :: info !! Status code for the pairwise computation.
      real(dp), intent(in), optional :: gamma_shape !! Optional positive gamma shape for supported correction models.
      real(dp) :: f1(4)
      real(dp) :: f2(4)
      real(dp) :: table(4, 4)
      real(dp) :: p
      real(dp) :: p1
      real(dp) :: p2
      real(dp) :: q
      real(dp) :: r
      real(dp) :: a1
      real(dp) :: a2
      real(dp) :: a3
      real(dp) :: alpha
      real(dp) :: aa
      real(dp) :: bb
      real(dp) :: cc
      real(dp) :: e
      real(dp) :: g_r
      real(dp) :: g_y
      real(dp) :: k1
      real(dp) :: k2
      real(dp) :: k3
      real(dp) :: k4
      real(dp) :: w1
      real(dp) :: w2
      real(dp) :: w3
      real(dp) :: wg
      real(dp) :: determinant
      real(dp) :: denominator
      integer :: differences
      integer :: i
      integer :: transitions
      integer :: transitions_ag
      integer :: transitions_ct
      integer :: transv_one
      integer :: transv_two
      integer :: transversions
      logical :: use_gamma
      character(len=16) :: chosen

      distance = ieee_value(0.0_dp, ieee_quiet_nan)
      comparable_sites = 0
      info = 0
      if (size(sequence_a) /= size(sequence_b)) then
         info = 1
         return
      end if
      chosen = uppercase(trim(adjustl(model)))
      if (trim(chosen) == 'INDEL') then
         distance = real(indel_count(sequence_a, sequence_b), dp)
         comparable_sites = size(sequence_a)
         return
      end if
      if (trim(chosen) == 'INDELBLOCK') then
         distance = real(indel_block_count(sequence_a, sequence_b), dp)
         comparable_sites = size(sequence_a)
         return
      end if

      differences = 0
      transitions = 0
      transitions_ag = 0
      transitions_ct = 0
      transv_one = 0
      transv_two = 0
      transversions = 0
      do i = 1, size(sequence_a)
         if (.not. known_base(sequence_a(i)) .or. .not. known_base(sequence_b(i))) cycle
         comparable_sites = comparable_sites + 1
         if (sequence_a(i) == sequence_b(i)) cycle
         differences = differences + 1
         if (is_transition(sequence_a(i), sequence_b(i))) then
            transitions = transitions + 1
            if (is_ag(sequence_a(i), sequence_b(i))) then
               transitions_ag = transitions_ag + 1
            else
               transitions_ct = transitions_ct + 1
            end if
         else
            transversions = transversions + 1
            if (is_transversion_one(sequence_a(i), sequence_b(i))) then
               transv_one = transv_one + 1
            else
               transv_two = transv_two + 1
            end if
         end if
      end do
      if (comparable_sites == 0) then
         info = 2
         return
      end if
      use_gamma = present(gamma_shape)
      if (use_gamma) then
         alpha = gamma_shape
         if (alpha <= 0.0_dp) then
            info = 3
            return
         end if
      else
         alpha = 1.0_dp
      end if
      p = real(transitions, dp) / real(comparable_sites, dp)
      q = real(transversions, dp) / real(comparable_sites, dp)

      select case (trim(chosen))
      case ('RAW')
         distance = real(differences, dp) / real(comparable_sites, dp)
      case ('N')
         distance = real(differences, dp)
      case ('TS')
         distance = real(transitions, dp)
      case ('TV')
         distance = real(transversions, dp)
      case ('JC69')
         a1 = 1.0_dp - 4.0_dp * real(differences, dp) / (3.0_dp * real(comparable_sites, dp))
         if (a1 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            distance = 0.75_dp * alpha * (a1**(-1.0_dp / alpha) - 1.0_dp)
         else
            distance = -0.75_dp * log(a1)
         end if
      case ('K80')
         a1 = 1.0_dp - 2.0_dp * p - q
         a2 = 1.0_dp - 2.0_dp * q
         if (a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            distance = 0.5_dp * alpha * (a1**(-1.0_dp / alpha) &
               + 0.5_dp * a2**(-1.0_dp / alpha) - 1.5_dp)
         else
            distance = -0.5_dp * log(a1 * sqrt(a2))
         end if
      case ('F81')
         e = 1.0_dp - sum(bf * bf)
         p = real(differences, dp) / real(comparable_sites, dp)
         if (e <= 0.0_dp .or. 1.0_dp - p / e <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            distance = e * alpha * ((1.0_dp - p / e)**(-1.0_dp / alpha) - 1.0_dp)
         else
            distance = -e * log(1.0_dp - p / e)
         end if
      case ('K81')
         p = real(transitions, dp) / real(comparable_sites, dp)
         q = real(transv_one, dp) / real(comparable_sites, dp)
         r = real(transv_two, dp) / real(comparable_sites, dp)
         a1 = 1.0_dp - 2.0_dp * p - 2.0_dp * q
         a2 = 1.0_dp - 2.0_dp * p - 2.0_dp * r
         a3 = 1.0_dp - 2.0_dp * q - 2.0_dp * r
         if (a1 <= 0.0_dp .or. a2 <= 0.0_dp .or. a3 <= 0.0_dp) then
            info = 4
            return
         end if
         distance = -0.25_dp * log(a1 * a2 * a3)
      case ('F84')
         if (bf(dna_a) + bf(dna_g) <= 0.0_dp .or. bf(dna_c) + bf(dna_t) <= 0.0_dp) then
            info = 4
            return
         end if
         aa = bf(dna_a) * bf(dna_g) / (bf(dna_a) + bf(dna_g)) &
            + bf(dna_c) * bf(dna_t) / (bf(dna_c) + bf(dna_t))
         bb = bf(dna_a) * bf(dna_g) + bf(dna_c) * bf(dna_t)
         cc = (bf(dna_a) + bf(dna_g)) * (bf(dna_c) + bf(dna_t))
         if (aa <= 0.0_dp .or. cc <= 0.0_dp) then
            info = 4
            return
         end if
         a1 = 1.0_dp - p / (2.0_dp * aa) - (aa - bb) * q / (2.0_dp * aa * cc)
         a2 = 1.0_dp - q / (2.0_dp * cc)
         if (a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
            info = 4
            return
         end if
         distance = -2.0_dp * aa * log(a1) + 2.0_dp * (aa - bb - cc) * log(a2)
      case ('T92')
         wg = 2.0_dp * (bf(dna_c) + bf(dna_g)) * (1.0_dp - bf(dna_c) - bf(dna_g))
         if (wg <= 0.0_dp .or. wg >= 1.0_dp) then
            info = 4
            return
         end if
         a1 = 1.0_dp - p / wg - q
         a2 = 1.0_dp - 2.0_dp * q
         if (a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
            info = 4
            return
         end if
         distance = -wg * log(a1) - 0.5_dp * (1.0_dp - wg) * log(a2)
      case ('TN93')
         g_r = bf(dna_a) + bf(dna_g)
         g_y = bf(dna_c) + bf(dna_t)
         if (g_r <= 0.0_dp .or. g_y <= 0.0_dp) then
            info = 4
            return
         end if
         k1 = 2.0_dp * bf(dna_a) * bf(dna_g) / g_r
         k2 = 2.0_dp * bf(dna_c) * bf(dna_t) / g_y
         k3 = 2.0_dp * (g_r * g_y - bf(dna_a) * bf(dna_g) * g_y / g_r &
            - bf(dna_c) * bf(dna_t) * g_r / g_y)
         p1 = real(transitions_ag, dp) / real(comparable_sites, dp)
         p2 = real(transitions_ct, dp) / real(comparable_sites, dp)
         if (k1 <= 0.0_dp .or. k2 <= 0.0_dp) then
            info = 4
            return
         end if
         w1 = 1.0_dp - p1 / k1 - q / (2.0_dp * g_r)
         w2 = 1.0_dp - p2 / k2 - q / (2.0_dp * g_y)
         w3 = 1.0_dp - q / (2.0_dp * g_r * g_y)
         if (w1 <= 0.0_dp .or. w2 <= 0.0_dp .or. w3 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            k4 = 2.0_dp * (bf(dna_a) * bf(dna_g) + bf(dna_c) * bf(dna_t) + g_r * g_y)
            distance = alpha * (k1 * w1**(-1.0_dp / alpha) + k2 * w2**(-1.0_dp / alpha) &
               + k3 * w3**(-1.0_dp / alpha) - k4)
         else
            distance = -k1 * log(w1) - k2 * log(w2) - k3 * log(w3)
         end if
      case ('LOGDET')
         call dna_contingency_table(sequence_a, sequence_b, table)
         table = table / real(comparable_sites, dp)
         determinant = det4(table)
         if (determinant <= 0.0_dp) then
            info = 4
            return
         end if
         distance = -0.25_dp * log(determinant) - log(4.0_dp)
      case ('PARALIN')
         call dna_contingency_table(sequence_a, sequence_b, table)
         table = table / real(comparable_sites, dp)
         determinant = det4(table)
         call single_base_frequencies(sequence_a, f1)
         call single_base_frequencies(sequence_b, f2)
         denominator = sqrt(product(f1) * product(f2))
         if (determinant <= 0.0_dp .or. denominator <= 0.0_dp) then
            info = 4
            return
         end if
         distance = -0.25_dp * log(determinant / denominator)
      case default
         info = 5
      end select
   end subroutine dna_distance_core

   pure subroutine dna_variance_core(sequence_a, sequence_b, model, bf, variance, info, gamma_shape)
      !! Evaluates analytical distance-variance formulas translated from ape `src/dist_dna.c`.
      integer, intent(in) :: sequence_a(:) !! First sequence in the pair.
      integer, intent(in) :: sequence_b(:) !! Second sequence with the same site count.
      character(len=*), intent(in) :: model !! Uppercase distance model name.
      real(dp), intent(in) :: bf(4) !! A,C,G,T frequencies used by frequency-dependent variance formulas.
      real(dp), intent(out) :: variance !! Analytical sampling variance, or IEEE NaN if unsupported/undefined.
      integer, intent(out) :: info !! Zero on success, 4 for an undefined correction, or 6 for unavailable variance.
      real(dp), intent(in), optional :: gamma_shape !! Optional positive gamma shape for JC69, K80, F81, or TN93.
      real(dp) :: p
      real(dp) :: p1
      real(dp) :: p2
      real(dp) :: q
      real(dp) :: r
      real(dp) :: a1
      real(dp) :: a2
      real(dp) :: a3
      real(dp) :: aa
      real(dp) :: bb
      real(dp) :: cc
      real(dp) :: alpha
      real(dp) :: b
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: c3
      real(dp) :: e
      real(dp) :: g_a2
      real(dp) :: g_c2
      real(dp) :: g_g2
      real(dp) :: g_t2
      real(dp) :: g_ag
      real(dp) :: g_ct
      real(dp) :: g_r
      real(dp) :: g_r2
      real(dp) :: g_y
      real(dp) :: g_y2
      real(dp) :: k1
      real(dp) :: k2
      real(dp) :: k4
      real(dp) :: t1
      real(dp) :: t2
      real(dp) :: t3
      real(dp) :: w1
      real(dp) :: w2
      real(dp) :: w3
      real(dp) :: wg
      real(dp) :: table(4, 4)
      real(dp) :: inverse(4, 4)
      real(dp) :: f1(4)
      real(dp) :: f2(4)
      real(dp) :: accumulator
      real(dp) :: correction
      integer :: comparable
      integer :: differences
      integer :: i
      integer :: column
      integer :: row
      integer :: solve_info
      integer :: transitions
      integer :: transitions_ag
      integer :: transitions_ct
      integer :: transv_one
      integer :: transv_two
      integer :: transversions
      logical :: use_gamma
      character(len=16) :: chosen

      variance = ieee_value(0.0_dp, ieee_quiet_nan)
      info = 0
      if (size(sequence_a) /= size(sequence_b)) then
         info = 1
         return
      end if
      comparable = 0
      differences = 0
      transitions = 0
      transitions_ag = 0
      transitions_ct = 0
      transv_one = 0
      transv_two = 0
      transversions = 0
      do i = 1, size(sequence_a)
         if (.not. known_base(sequence_a(i)) .or. .not. known_base(sequence_b(i))) cycle
         comparable = comparable + 1
         if (sequence_a(i) == sequence_b(i)) cycle
         differences = differences + 1
         if (is_transition(sequence_a(i), sequence_b(i))) then
            transitions = transitions + 1
            if (is_ag(sequence_a(i), sequence_b(i))) then
               transitions_ag = transitions_ag + 1
            else
               transitions_ct = transitions_ct + 1
            end if
         else
            transversions = transversions + 1
            if (is_transversion_one(sequence_a(i), sequence_b(i))) then
               transv_one = transv_one + 1
            else
               transv_two = transv_two + 1
            end if
         end if
      end do
      if (comparable == 0) then
         info = 2
         return
      end if
      use_gamma = present(gamma_shape)
      if (use_gamma) then
         alpha = gamma_shape
         if (alpha <= 0.0_dp) then
            info = 3
            return
         end if
      else
         alpha = 1.0_dp
      end if
      chosen = uppercase(trim(adjustl(model)))
      p = real(transitions, dp) / real(comparable, dp)
      q = real(transversions, dp) / real(comparable, dp)

      select case (trim(chosen))
      case ('JC69')
         p = real(differences, dp) / real(comparable, dp)
         a1 = 1.0_dp - 4.0_dp * p / 3.0_dp
         if (a1 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            variance = p * (1.0_dp - p) / (a1**(-2.0_dp / (alpha + 1.0_dp)) * real(comparable, dp))
         else
            variance = p * (1.0_dp - p) / (a1 * a1 * real(comparable, dp))
         end if
      case ('K80')
         a1 = 1.0_dp - 2.0_dp * p - q
         a2 = 1.0_dp - 2.0_dp * q
         if (a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            b = -(1.0_dp / alpha + 1.0_dp)
            c1 = a1**b
            c2 = a2**b
            c3 = 0.5_dp * (c1 + c2)
         else
            c1 = 1.0_dp / a1
            c2 = 1.0_dp / a2
            c3 = 0.5_dp * (c1 + c2)
         end if
         variance = (c1 * c1 * p + c3 * c3 * q - (c1 * p + c3 * q)**2) / real(comparable, dp)
      case ('F81')
         p = real(differences, dp) / real(comparable, dp)
         e = 1.0_dp - sum(bf * bf)
         a1 = 1.0_dp - p / e
         if (e <= 0.0_dp .or. a1 <= 0.0_dp) then
            info = 4
            return
         end if
         if (use_gamma) then
            variance = p * (1.0_dp - p) / (a1**(-2.0_dp / (alpha + 1.0_dp)) * real(comparable, dp))
         else
            variance = p * (1.0_dp - p) / (a1 * a1 * real(comparable, dp))
         end if
      case ('K81')
         p = real(transitions, dp) / real(comparable, dp)
         q = real(transv_one, dp) / real(comparable, dp)
         r = real(transv_two, dp) / real(comparable, dp)
         a1 = 1.0_dp - 2.0_dp * p - 2.0_dp * q
         a2 = 1.0_dp - 2.0_dp * p - 2.0_dp * r
         a3 = 1.0_dp - 2.0_dp * q - 2.0_dp * r
         if (min(a1, a2, a3) <= 0.0_dp) then
            info = 4
            return
         end if
         aa = 0.5_dp * (1.0_dp / a1 + 1.0_dp / a2)
         bb = 0.5_dp * (1.0_dp / a1 + 1.0_dp / a3)
         cc = 0.5_dp * (1.0_dp / a2 + 1.0_dp / a3)
         variance = (aa * aa * p + bb * bb * q + cc * cc * r - (aa * p + bb * q + cc * r)**2) / 2.0_dp
      case ('F84')
         g_r = bf(dna_a) + bf(dna_g)
         g_y = bf(dna_c) + bf(dna_t)
         if (g_r <= 0.0_dp .or. g_y <= 0.0_dp) then
            info = 4
            return
         end if
         aa = bf(dna_a) * bf(dna_g) / g_r + bf(dna_c) * bf(dna_t) / g_y
         bb = bf(dna_a) * bf(dna_g) + bf(dna_c) * bf(dna_t)
         cc = g_r * g_y
         t1 = aa * cc
         t2 = cc * p / 2.0_dp
         t3 = (aa - bb) * q / 2.0_dp
         if (t1 - t2 - t3 <= 0.0_dp .or. cc - q / 2.0_dp <= 0.0_dp) then
            info = 4
            return
         end if
         a1 = t1 / (t1 - t2 - t3)
         a2 = aa * (aa - bb) / (t1 - t2 - t3) - (aa - bb - cc) / (cc - q / 2.0_dp)
         variance = (a1 * a1 * p + a2 * a2 * q - (a1 * p + a2 * q)**2) / real(comparable, dp)
      case ('T92')
         wg = 2.0_dp * (bf(dna_c) + bf(dna_g)) * (1.0_dp - bf(dna_c) - bf(dna_g))
         a1 = 1.0_dp - p / wg - q
         a2 = 1.0_dp - 2.0_dp * q
         if (wg <= 0.0_dp .or. a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
            info = 4
            return
         end if
         c1 = 1.0_dp / a1
         c2 = 1.0_dp / a2
         c3 = wg * (c1 - c2) + c2
         variance = (c1 * c1 * p + c3 * c3 * q - (c1 * p + c3 * q)**2) / real(comparable, dp)
      case ('TN93')
         g_r = bf(dna_a) + bf(dna_g)
         g_y = bf(dna_c) + bf(dna_t)
         if (g_r <= 0.0_dp .or. g_y <= 0.0_dp) then
            info = 4
            return
         end if
         k1 = 2.0_dp * bf(dna_a) * bf(dna_g) / g_r
         k2 = 2.0_dp * bf(dna_c) * bf(dna_t) / g_y
         p1 = real(transitions_ag, dp) / real(comparable, dp)
         p2 = real(transitions_ct, dp) / real(comparable, dp)
         w1 = 1.0_dp - p1 / k1 - q / (2.0_dp * g_r)
         w2 = 1.0_dp - p2 / k2 - q / (2.0_dp * g_y)
         w3 = 1.0_dp - q / (2.0_dp * g_r * g_y)
         if (k1 <= 0.0_dp .or. k2 <= 0.0_dp .or. min(w1, w2, w3) <= 0.0_dp) then
            info = 4
            return
         end if
         g_a2 = bf(dna_a) * bf(dna_a)
         g_c2 = bf(dna_c) * bf(dna_c)
         g_g2 = bf(dna_g) * bf(dna_g)
         g_t2 = bf(dna_t) * bf(dna_t)
         g_ag = bf(dna_a) * bf(dna_g)
         g_ct = bf(dna_c) * bf(dna_t)
         g_r2 = g_r * g_r
         g_y2 = g_y * g_y
         if (use_gamma) then
            b = -(1.0_dp + 1.0_dp / alpha)
            c1 = w1**b
            c2 = w2**b
            c3 = g_ag * c1 / g_r2 + g_ct * c2 / g_y2 + &
               ((g_a2 + g_g2) / (2.0_dp * g_r2) + (g_t2 + g_c2) / (2.0_dp * g_y2)) * w3**b
         else
            c1 = 1.0_dp / w1
            c2 = 1.0_dp / w2
            c3 = 2.0_dp * g_a2 * g_g2 / (g_r * (2.0_dp * g_ag * g_r - g_r2 * p1 - g_ag * q)) + &
               2.0_dp * g_c2 * g_t2 / (g_y * (2.0_dp * g_ct * g_y - g_y2 * p2 - g_ct * q)) + &
               (g_r2 * (g_t2 + g_c2) + g_y2 * (g_a2 + g_g2)) / (2.0_dp * g_r2 * g_y2 - g_r * g_y * q)
         end if
         k4 = c1 * p1 + c2 * p2 + c3 * q
         variance = (c1 * c1 * p1 + c2 * c2 * p2 + c3 * c3 * q - k4 * k4) / real(comparable, dp)
      case ('LOGDET')
         call dna_contingency_table(sequence_a, sequence_b, table)
         table = table / real(comparable, dp)
         call lapack_style_inverse4(table, inverse, solve_info)
         if (solve_info /= 0) then
            info = 4
            return
         end if
         accumulator = 0.0_dp
         do column = 1, 4
            do row = 1, 4
               accumulator = accumulator + inverse(row, column)**2 * table(column, row)
            end do
         end do
         variance = (accumulator - 16.0_dp) / (16.0_dp * real(comparable, dp))
      case ('PARALIN')
         call dna_contingency_table(sequence_a, sequence_b, table)
         table = table / real(comparable, dp)
         call lapack_style_inverse4(table, inverse, solve_info)
         call single_base_frequencies(sequence_a, f1)
         call single_base_frequencies(sequence_b, f2)
         if (solve_info /= 0 .or. any(f1 <= 0.0_dp) .or. any(f2 <= 0.0_dp)) then
            info = 4
            return
         end if
         accumulator = 0.0_dp
         do column = 1, 4
            do row = 1, 4
               accumulator = accumulator + inverse(row, column)**2 * table(column, row)
            end do
         end do
         correction = 4.0_dp * sum(1.0_dp / sqrt(f1 * f2))
         variance = (accumulator - correction) / (16.0_dp * real(comparable, dp))
      case default
         info = 6
      end select
   end subroutine dna_variance_core

   pure elemental logical function variance_model_supported(model) result(supported)
      !! Reports whether this translation currently provides ape's analytical sampling-variance formula for a model.
      character(len=*), intent(in) :: model !! Uppercase or mixed-case ape distance model name.
      character(len=16) :: chosen

      chosen = uppercase(trim(adjustl(model)))
      select case (trim(chosen))
      case ('JC69', 'K80', 'F81', 'K81', 'F84', 'T92', 'TN93', 'GG95', 'LOGDET', 'PARALIN')
         supported = .true.
      case default
         supported = .false.
      end select
   end function variance_model_supported

   pure subroutine gg95_matrix(sequences, distance, info, variance)
      integer, intent(in) :: sequences(:, :) !! Working taxa-by-site DNA alignment.
      real(dp), intent(inout) :: distance(:, :) !! Preallocated symmetric matrix receiving GG95 distances.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if alpha or a pairwise correction is undefined.
      real(dp), intent(inout), optional :: variance(:, :) !! Optional symmetric GG95 analytical sampling-variance matrix.
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: p(:)
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: ratio(:)
      integer, allocatable :: pair_i(:)
      integer, allocatable :: pair_j(:)
      integer, allocatable :: length_pair(:)
      integer :: comparable
      integer :: differences
      integer :: finite_count
      integer :: gc
      integer :: i
      integer :: j
      integer :: k
      integer :: npair
      integer :: transitions
      real(dp) :: a
      real(dp) :: k1
      real(dp) :: k2
      real(dp) :: mean_alpha
      real(dp) :: sum_ratio

      info = 0
      npair = size(sequences, 1) * (size(sequences, 1) - 1) / 2
      allocate(theta(size(sequences, 1)), p(npair), q(npair), ratio(npair))
      allocate(pair_i(npair), pair_j(npair), length_pair(npair))
      theta = 0.0_dp
      do i = 1, size(sequences, 1)
         comparable = 0
         gc = 0
         do j = 1, size(sequences, 2)
            if (.not. known_base(sequences(i, j))) cycle
            comparable = comparable + 1
            if (sequences(i, j) == dna_c .or. sequences(i, j) == dna_g) gc = gc + 1
         end do
         if (comparable == 0) then
            info = 1
            return
         end if
         theta(i) = real(gc, dp) / real(comparable, dp)
      end do

      k = 0
      do i = 1, size(sequences, 1) - 1
         do j = i + 1, size(sequences, 1)
            k = k + 1
            pair_i(k) = i
            pair_j(k) = j
            call pair_change_counts(sequences(i, :), sequences(j, :), differences, transitions, comparable)
            length_pair(k) = comparable
            if (comparable == 0) then
               ratio(k) = ieee_value(0.0_dp, ieee_quiet_nan)
               cycle
            end if
            p(k) = real(transitions, dp) / real(comparable, dp)
            q(k) = real(differences - transitions, dp) / real(comparable, dp)
            if (1.0_dp - 2.0_dp * q(k) <= 0.0_dp .or. 1.0_dp - 2.0_dp * p(k) - q(k) <= 0.0_dp) then
               ratio(k) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               a = log(1.0_dp - 2.0_dp * q(k))
               if (abs(a) <= tiny(1.0_dp)) then
                  ratio(k) = ieee_value(0.0_dp, ieee_quiet_nan)
               else
                  ratio(k) = 2.0_dp * (log(1.0_dp - 2.0_dp * p(k) - q(k)) - 0.5_dp * a) / a
               end if
            end if
         end do
      end do
      sum_ratio = 0.0_dp
      finite_count = 0
      do k = 1, npair
         if (ieee_is_finite(ratio(k))) then
            sum_ratio = sum_ratio + ratio(k)
            finite_count = finite_count + 1
         end if
      end do
      if (finite_count == 0) then
         info = 2
         return
      end if
      mean_alpha = sum_ratio / real(finite_count, dp)
      if (mean_alpha <= -1.0_dp) then
         info = 3
         return
      end if
      do k = 1, npair
         i = pair_i(k)
         j = pair_j(k)
         if (length_pair(k) == 0 .or. 1.0_dp - 2.0_dp * q(k) <= 0.0_dp) then
            distance(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
            distance(j, i) = distance(i, j)
            info = max(info, 4)
            cycle
         end if
         a = 1.0_dp - 2.0_dp * q(k)
         k1 = 1.0_dp + mean_alpha * (theta(i) * (1.0_dp - theta(i)) + theta(j) * (1.0_dp - theta(j)))
         k2 = mean_alpha * (theta(i) - theta(j))**2 / (mean_alpha + 1.0_dp)
         distance(i, j) = -0.5_dp * k1 * log(a) + k2 * (1.0_dp - a**(0.25_dp * (mean_alpha + 1.0_dp)))
         distance(j, i) = distance(i, j)
         if (present(variance)) then
            variance(i, j) = (k1 + 0.5_dp * k2 * (mean_alpha + 1.0_dp) * &
               a**(0.25_dp * (mean_alpha + 1.0_dp)))**2 * q(k) * (1.0_dp - q(k)) / &
               (a * a * real(length_pair(k), dp))
            variance(j, i) = variance(i, j)
         end if
      end do
   end subroutine gg95_matrix

   pure subroutine pair_change_counts(sequence_a, sequence_b, differences, transitions, comparable)
      integer, intent(in) :: sequence_a(:) !! First DNA sequence for pairwise event counting.
      integer, intent(in) :: sequence_b(:) !! Second DNA sequence of matching length.
      integer, intent(out) :: differences !! Number of unequal comparable sites.
      integer, intent(out) :: transitions !! Number of transition substitutions among comparable sites.
      integer, intent(out) :: comparable !! Number of sites with known bases in both sequences.
      integer :: i

      differences = 0
      transitions = 0
      comparable = 0
      do i = 1, min(size(sequence_a), size(sequence_b))
         if (.not. known_base(sequence_a(i)) .or. .not. known_base(sequence_b(i))) cycle
         comparable = comparable + 1
         if (sequence_a(i) == sequence_b(i)) cycle
         differences = differences + 1
         if (is_transition(sequence_a(i), sequence_b(i))) transitions = transitions + 1
      end do
   end subroutine pair_change_counts

   pure subroutine pair_base_frequencies(sequence_a, sequence_b, frequencies)
      integer, intent(in) :: sequence_a(:) !! First sequence contributing to pair-level default base frequencies.
      integer, intent(in) :: sequence_b(:) !! Second sequence contributing to pair-level default base frequencies.
      real(dp), intent(out) :: frequencies(4) !! Combined A,C,G,T frequencies over all known states in both sequences.
      integer :: i
      integer :: total

      frequencies = 0.0_dp
      total = 0
      do i = 1, size(sequence_a)
         if (known_base(sequence_a(i))) then
            frequencies(sequence_a(i)) = frequencies(sequence_a(i)) + 1.0_dp
            total = total + 1
         end if
      end do
      do i = 1, size(sequence_b)
         if (known_base(sequence_b(i))) then
            frequencies(sequence_b(i)) = frequencies(sequence_b(i)) + 1.0_dp
            total = total + 1
         end if
      end do
      if (total > 0) frequencies = frequencies / real(total, dp)
   end subroutine pair_base_frequencies

   pure subroutine single_base_frequencies(sequence, frequencies)
      integer, intent(in) :: sequence(:) !! DNA sequence whose known-state frequencies are requested.
      real(dp), intent(out) :: frequencies(4) !! A,C,G,T frequencies normalized over known states only.
      integer :: i
      integer :: total

      frequencies = 0.0_dp
      total = 0
      do i = 1, size(sequence)
         if (.not. known_base(sequence(i))) cycle
         frequencies(sequence(i)) = frequencies(sequence(i)) + 1.0_dp
         total = total + 1
      end do
      if (total > 0) frequencies = frequencies / real(total, dp)
   end subroutine single_base_frequencies

   pure function pack_columns(sequences, keep) result(work)
      integer, intent(in) :: sequences(:, :) !! Source taxa-by-site alignment.
      logical, intent(in) :: keep(:) !! Site mask with one entry per source column.
      integer, allocatable :: work(:, :)
      integer :: j
      integer :: k

      allocate(work(size(sequences, 1), count(keep)))
      k = 0
      do j = 1, min(size(sequences, 2), size(keep))
         if (.not. keep(j)) cycle
         k = k + 1
         work(:, k) = sequences(:, j)
      end do
   end function pack_columns

   pure integer function indel_count(sequence_a, sequence_b) result(count_value)
      integer, intent(in) :: sequence_a(:) !! First sequence in which gap status is compared sitewise.
      integer, intent(in) :: sequence_b(:) !! Second sequence of equal length.
      integer :: i

      count_value = 0
      do i = 1, min(size(sequence_a), size(sequence_b))
         if ((sequence_a(i) == dna_gap) .neqv. (sequence_b(i) == dna_gap)) count_value = count_value + 1
      end do
   end function indel_count

   pure integer function indel_block_count(sequence_a, sequence_b) result(count_value)
      integer, intent(in) :: sequence_a(:) !! First sequence whose gap runs are encoded by start position and length.
      integer, intent(in) :: sequence_b(:) !! Second sequence whose gap-run encoding is compared with the first.
      integer, allocatable :: block_a(:)
      integer, allocatable :: block_b(:)
      integer :: i

      call encode_indel_blocks(sequence_a, block_a)
      call encode_indel_blocks(sequence_b, block_b)
      count_value = 0
      do i = 1, min(size(block_a), size(block_b))
         if (block_a(i) /= block_b(i)) count_value = count_value + 1
      end do
   end function indel_block_count

   pure subroutine encode_indel_blocks(sequence, blocks)
      integer, intent(in) :: sequence(:) !! DNA sequence in which consecutive `dna_gap` codes define an indel block.
      integer, allocatable, intent(out) :: blocks(:) !! Zero except at each gap-run start, where the run length is stored.
      integer :: i
      integer :: run_length
      integer :: start

      allocate(blocks(size(sequence)))
      blocks = 0
      i = 1
      do while (i <= size(sequence))
         if (sequence(i) /= dna_gap) then
            i = i + 1
            cycle
         end if
         start = i
         run_length = 0
         do while (i <= size(sequence))
            if (sequence(i) /= dna_gap) exit
            run_length = run_length + 1
            i = i + 1
         end do
         blocks(start) = run_length
      end do
   end subroutine encode_indel_blocks

   pure subroutine lapack_style_inverse4(matrix_lu, inverse, info)
      !! Reproduces the 4-by-4 DGETRF/DGETRS state used by ape LogDet and ParaLin variance calculations.
      real(dp), intent(inout) :: matrix_lu(4, 4) !! Input matrix, overwritten by LAPACK-compatible combined L/U factors.
      real(dp), intent(out) :: inverse(4, 4) !! Matrix inverse obtained by solving against the identity matrix.
      integer, intent(out) :: info !! Zero on success, or the one-based singular pivot index.
      real(dp) :: row_buffer(4)
      real(dp) :: rhs_buffer(4)
      integer :: pivot(4)
      integer :: i
      integer :: j
      integer :: jp

      info = 0
      inverse = 0.0_dp
      do i = 1, 4
         inverse(i, i) = 1.0_dp
      end do
      do j = 1, 4
         jp = j - 1 + maxloc(abs(matrix_lu(j:4, j)), dim=1)
         pivot(j) = jp
         if (.not. (abs(matrix_lu(jp, j)) > 0.0_dp)) then
            info = j
            return
         end if
         if (jp /= j) then
            row_buffer = matrix_lu(j, :)
            matrix_lu(j, :) = matrix_lu(jp, :)
            matrix_lu(jp, :) = row_buffer
         end if
         if (j < 4) then
            matrix_lu(j + 1:4, j) = matrix_lu(j + 1:4, j) / matrix_lu(j, j)
            do i = j + 1, 4
               matrix_lu(i, j + 1:4) = matrix_lu(i, j + 1:4) - matrix_lu(i, j) * matrix_lu(j, j + 1:4)
            end do
         end if
      end do
      do j = 1, 4
         jp = pivot(j)
         if (jp /= j) then
            rhs_buffer = inverse(j, :)
            inverse(j, :) = inverse(jp, :)
            inverse(jp, :) = rhs_buffer
         end if
      end do
      do j = 1, 4
         if (j < 4) then
            do i = j + 1, 4
               inverse(i, :) = inverse(i, :) - matrix_lu(i, j) * inverse(j, :)
            end do
         end if
      end do
      do j = 4, 1, -1
         inverse(j, :) = inverse(j, :) / matrix_lu(j, j)
         if (j > 1) then
            do i = 1, j - 1
               inverse(i, :) = inverse(i, :) - matrix_lu(i, j) * inverse(j, :)
            end do
         end if
      end do
   end subroutine lapack_style_inverse4

   pure real(dp) function det4(matrix) result(det)
      real(dp), intent(in) :: matrix(4, 4) !! Four-by-four matrix whose determinant is evaluated explicitly.
      integer :: i
      real(dp) :: minor(3, 3)
      real(dp) :: sign

      det = 0.0_dp
      sign = 1.0_dp
      do i = 1, 4
         call first_row_minor(matrix, i, minor)
         det = det + sign * matrix(1, i) * det3(minor)
         sign = -sign
      end do
   end function det4

   pure subroutine first_row_minor(matrix, removed_column, minor)
      real(dp), intent(in) :: matrix(4, 4) !! Source matrix for a first-row cofactor minor.
      integer, intent(in) :: removed_column !! Column removed together with the first row, in the range 1 through 4.
      real(dp), intent(out) :: minor(3, 3) !! Resulting three-by-three cofactor matrix.
      integer :: column
      integer :: target

      target = 0
      do column = 1, 4
         if (column == removed_column) cycle
         target = target + 1
         minor(:, target) = matrix(2:4, column)
      end do
   end subroutine first_row_minor

   pure real(dp) function det3(matrix) result(det)
      real(dp), intent(in) :: matrix(3, 3) !! Three-by-three matrix whose determinant is evaluated explicitly.

      det = matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) - matrix(2, 3) * matrix(3, 2)) &
         - matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) - matrix(2, 3) * matrix(3, 1)) &
         + matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) - matrix(2, 2) * matrix(3, 1))
   end function det3

   pure elemental logical function known_base(base) result(known)
      integer, intent(in) :: base !! Integer nucleotide code tested for unambiguous A/C/G/T membership.

      known = base >= dna_a .and. base <= dna_t
   end function known_base

   pure elemental logical function is_transition(base_a, base_b) result(transition)
      integer, intent(in) :: base_a !! First unambiguous nucleotide code.
      integer, intent(in) :: base_b !! Second unambiguous nucleotide code.

      transition = is_ag(base_a, base_b) .or. is_ct(base_a, base_b)
   end function is_transition

   pure elemental logical function is_ag(base_a, base_b) result(value)
      integer, intent(in) :: base_a !! First unambiguous nucleotide code.
      integer, intent(in) :: base_b !! Second unambiguous nucleotide code.

      value = (base_a == dna_a .and. base_b == dna_g) .or. &
         (base_a == dna_g .and. base_b == dna_a)
   end function is_ag

   pure elemental logical function is_ct(base_a, base_b) result(value)
      integer, intent(in) :: base_a !! First unambiguous nucleotide code.
      integer, intent(in) :: base_b !! Second unambiguous nucleotide code.

      value = (base_a == dna_c .and. base_b == dna_t) .or. &
         (base_a == dna_t .and. base_b == dna_c)
   end function is_ct

   pure elemental logical function is_transversion_one(base_a, base_b) result(value)
      integer, intent(in) :: base_a !! First unambiguous nucleotide code.
      integer, intent(in) :: base_b !! Second unambiguous nucleotide code.

      value = ((base_a == dna_a .and. base_b == dna_t) .or. (base_a == dna_t .and. base_b == dna_a)) &
         .or. ((base_a == dna_g .and. base_b == dna_c) .or. (base_a == dna_c .and. base_b == dna_g))
   end function is_transversion_one

   pure elemental function uppercase(text) result(out)
      character(len=*), intent(in) :: text !! ASCII model name converted to uppercase for case-insensitive dispatch.
      character(len=len(text)) :: out
      integer :: code
      integer :: i

      out = text
      do i = 1, len(text)
         code = iachar(out(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) out(i:i) = achar(code - 32)
      end do
   end function uppercase

   pure elemental integer function dna_state_slot(base) result(slot)
      integer, intent(in) :: base !! Public DNA state code mapped to ape's 17-state base-frequency order.

      select case (base)
      case (dna_a)
         slot = 1
      case (dna_c)
         slot = 2
      case (dna_g)
         slot = 3
      case (dna_t)
         slot = 4
      case (dna_r)
         slot = 5
      case (dna_m)
         slot = 6
      case (dna_w)
         slot = 7
      case (dna_s)
         slot = 8
      case (dna_k)
         slot = 9
      case (dna_y)
         slot = 10
      case (dna_v)
         slot = 11
      case (dna_h)
         slot = 12
      case (dna_d)
         slot = 13
      case (dna_b)
         slot = 14
      case (dna_n)
         slot = 15
      case (dna_gap)
         slot = 16
      case (dna_unknown)
         slot = 17
      case default
         slot = 0
      end select
   end function dna_state_slot



   pure elemental integer function dna_raw_code(base) result(raw)
      integer, intent(in) :: base !! Public DNA state mapped to ape's historical DNAbin raw-byte encoding.

      select case (base)
      case (dna_a)
         raw = 136
      case (dna_c)
         raw = 40
      case (dna_g)
         raw = 72
      case (dna_t)
         raw = 24
      case (dna_r)
         raw = 192
      case (dna_m)
         raw = 160
      case (dna_w)
         raw = 144
      case (dna_s)
         raw = 96
      case (dna_k)
         raw = 80
      case (dna_y)
         raw = 48
      case (dna_v)
         raw = 224
      case (dna_h)
         raw = 176
      case (dna_d)
         raw = 208
      case (dna_b)
         raw = 112
      case (dna_n)
         raw = 240
      case (dna_gap)
         raw = 4
      case (dna_unknown)
         raw = 2
      case default
         raw = 0
      end select
   end function dna_raw_code

   pure elemental logical function dna_states_surely_different(base_a, base_b) result(different)
      integer, intent(in) :: base_a !! First public DNA state compared with ape's `DifferentBase` bit rule.
      integer, intent(in) :: base_b !! Second public DNA state compared with ape's `DifferentBase` bit rule.
      integer :: raw_a
      integer :: raw_b

      raw_a = dna_raw_code(base_a)
      raw_b = dna_raw_code(base_b)
      different = iand(raw_a, raw_b) < 16
   end function dna_states_surely_different

   pure elemental character(len=1) function codon_amino(first, second, third, genetic_code) result(amino)
      integer, intent(in) :: first !! First nucleotide state of a potentially ambiguous codon.
      integer, intent(in) :: second !! Second nucleotide state of a potentially ambiguous codon.
      integer, intent(in) :: third !! Third nucleotide state of a potentially ambiguous codon.
      integer, intent(in) :: genetic_code !! Ape genetic-code table number 1 through 6.
      character(len=4913), parameter :: code_1 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXXRSRSRRRRRSRRRRRSSIIMIIIIIIIIIIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXX*CWCXXXXXCXXXXXCCLFLFLLLLLFLLLLLFF*XXXXXXXXXXXXXXXX*XXXXXXX' // &
         'XXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*X' // &
         'XXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXRXRXRRRRRXRRRRRXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      character(len=4913), parameter :: code_2 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXX*S*S*****S*****SSMIMIMMMMMIMMMMMIIXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXXWCWCWWWWWCWWWWWCCLFLFLLLLLFLLLLLFFXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      character(len=4913), parameter :: code_3 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXXRSRSRRRRRSRRRRRSSMIMIMMMMMIMMMMMIIXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXXWCWCWWWWWCWWWWWCCLFLFLLLLLFLLLLLFF*XXXXXXXXXXXXXXXX*XXXXXXX' // &
         'XXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*X' // &
         'XXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXRXRXRRRRRXRRRRRXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      character(len=4913), parameter :: code_4 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXXRSRSRRRRRSRRRRRSSIIMIIIIIIIIIIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXXWCWCWWWWWCWWWWWCCLFLFLLLLLFLLLLLFF*XXXXXXXXXXXXXXXX*XXXXXXX' // &
         'XXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*X' // &
         'XXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXRXRXRRRRRXRRRRRXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      character(len=4913), parameter :: code_5 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXXSSSSSSSSSSSSSSSXXMIMIMMMMMIMMMMMIIXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXXWCWCWWWWWCWWWWWCCLFLFLLLLLFLLLLLFF*XXXXXXXXXXXXXXXX*XXXXXXX' // &
         'XXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*X' // &
         'XXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXRXRXRRRRRXRRRRRXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      character(len=4913), parameter :: code_6 = &
         'KNKNKKKKKNKKKKKNNTTTTTTTTTTTTTTTXXRSRSRRRRRSRRRRRSSIIMIIIIIIIIIIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XQHQHQQQQQHQQQQQHHPPPPPPPPPPPPPPPXXRRRRRRRRRRRRRRRXXLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXEDEDEEEEEDEEEEEDDAAAAAAAAAAAAAAAXXGGGGGGGGGGGGGGGXXVVVVVVVVVVVVVVVXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXX*Y*Y*****Y*****YYSSSSSSSSSSSSSSSXX*CWCXXXXXCXXXXXCCLFLFLLLLLFLLLLLFF*XXXXXXXXXXXXXXXX*XXXXXXX' // &
         'XXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*X' // &
         'XXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXX*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXRXRXRRRRRXRRRRRXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLXLXLLLLLXLLLLLXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' // &
         'XXXXXXXXXXXXXXXXX'
      integer :: index
      integer :: slot_first
      integer :: slot_second
      integer :: slot_third

      slot_first = dna_state_slot(first)
      slot_second = dna_state_slot(second)
      slot_third = dna_state_slot(third)
      if (slot_first == 0 .or. slot_second == 0 .or. slot_third == 0) then
         amino = 'X'
         return
      end if
      index = (slot_first - 1) * 17 * 17 + (slot_second - 1) * 17 + slot_third
      select case (genetic_code)
      case (1)
         amino = code_1(index:index)
      case (2)
         amino = code_2(index:index)
      case (3)
         amino = code_3(index:index)
      case (4)
         amino = code_4(index:index)
      case (5)
         amino = code_5(index:index)
      case (6)
         amino = code_6(index:index)
      case default
         amino = 'X'
      end select
   end function codon_amino

end module ape_dna
