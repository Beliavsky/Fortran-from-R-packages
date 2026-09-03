module grain_cpt
  use r_kinds, only : dp
  use grbase_types, only : table_t
  use grbase_tables, only : make_table, valid_table, table_normalize_first
  implicit none
  private

  public :: make_cpt
  public :: logical_and_cpt
  public :: logical_or_cpt
  public :: mendel_cpt
  public :: mendel_probability
  public :: replace_cpt_values

contains

  pure function make_cpt(var, dim, values, smooth) result(cpt)
    integer, intent(in) :: var(:) !! Node labels with the child first and parent labels following in conditioning order.
    integer, intent(in) :: dim(:) !! State counts corresponding positionally to `var`; all counts must be positive.
    real(dp), intent(in) :: values(:) !! Flattened conditional weights in Fortran/R column-major order.
    real(dp), value, optional :: smooth !! Nonnegative cell pseudocount used before normalization; default zero.
    type(table_t) :: cpt
    type(table_t) :: raw
    real(dp), allocatable :: work(:)
    real(dp) :: s

    if (size(var) == 0 .or. size(var) /= size(dim)) return
    if (size(values) /= product(dim)) return
    s = 0.0_dp
    if (present(smooth)) s = smooth
    if (s < 0.0_dp) return
    allocate(work(size(values)))
    work = values + s
    if (any(work < 0.0_dp)) return
    raw = make_table(var, dim, work)
    if (.not. valid_table(raw)) return
    cpt = table_normalize_first(raw)
  end function make_cpt

  pure function logical_and_cpt(child, parent1, parent2) result(cpt)
    integer, value :: child !! Binary child node label; state 1 denotes TRUE and state 2 denotes FALSE.
    integer, value :: parent1 !! First binary parent node label using TRUE/FALSE states 1/2.
    integer, value :: parent2 !! Second binary parent node label using TRUE/FALSE states 1/2.
    type(table_t) :: cpt
    real(dp), parameter :: values(8) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                                        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]

    cpt = make_cpt([child, parent1, parent2], [2, 2, 2], values)
  end function logical_and_cpt

  pure function logical_or_cpt(child, parent1, parent2) result(cpt)
    integer, value :: child !! Binary child node label; state 1 denotes TRUE and state 2 denotes FALSE.
    integer, value :: parent1 !! First binary parent node label using TRUE/FALSE states 1/2.
    integer, value :: parent2 !! Second binary parent node label using TRUE/FALSE states 1/2.
    type(table_t) :: cpt
    real(dp), parameter :: values(8) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]

    cpt = make_cpt([child, parent1, parent2], [2, 2, 2], values)
  end function logical_or_cpt

  pure function mendel_cpt(n_alleles, child, father, mother) result(cpt)
    integer, value :: n_alleles !! Number of distinct alleles; must be at least one.
    integer, value :: child !! Child genotype node label.
    integer, value :: father !! Father genotype node label.
    integer, value :: mother !! Mother genotype node label.
    type(table_t) :: cpt
    real(dp), allocatable :: values(:)
    integer :: child_state
    integer :: father_state
    integer :: mother_state
    integer :: n_genotypes
    integer :: entry

    if (n_alleles < 1) return
    n_genotypes = n_alleles * (n_alleles + 1) / 2
    allocate(values(n_genotypes**3))
    entry = 0
    do mother_state = 1, n_genotypes
      do father_state = 1, n_genotypes
        do child_state = 1, n_genotypes
          entry = entry + 1
          values(entry) = mendel_probability(child_state, mother_state, father_state, n_alleles)
        end do
      end do
    end do
    cpt = make_cpt([child, father, mother], [n_genotypes, n_genotypes, n_genotypes], values)
  end function mendel_cpt

  pure real(dp) function mendel_probability(child_state, mother_state, father_state, n_alleles) result(probability)
    integer, value :: child_state !! One-based unordered child-genotype state index.
    integer, value :: mother_state !! One-based unordered mother-genotype state index.
    integer, value :: father_state !! One-based unordered father-genotype state index.
    integer, value :: n_alleles !! Number of possible alleles defining the genotype state space.
    integer :: child_pair(2)
    integer :: father_pair(2)
    integer :: mother_pair(2)
    real(dp) :: p1
    real(dp) :: p2

    probability = 0.0_dp
    child_pair = genotype_pair(child_state, n_alleles)
    mother_pair = genotype_pair(mother_state, n_alleles)
    father_pair = genotype_pair(father_state, n_alleles)
    if (any(child_pair == 0) .or. any(mother_pair == 0) .or. any(father_pair == 0)) return

    p1 = transmission_probability(child_pair(1), mother_pair) * &
         transmission_probability(child_pair(2), father_pair)
    if (child_pair(1) == child_pair(2)) then
      probability = p1
    else
      p2 = transmission_probability(child_pair(1), father_pair) * &
           transmission_probability(child_pair(2), mother_pair)
      probability = p1 + p2
    end if
  end function mendel_probability

  pure subroutine replace_cpt_values(cpt, values, ok)
    type(table_t), intent(inout) :: cpt !! Existing CPT whose domain is retained while its numerical cell values are replaced.
    real(dp), intent(in) :: values(:) !! Replacement conditional weights; length must equal the current CPT cell count.
    logical, intent(out) :: ok !! True when replacement values are nonnegative, correctly sized, and normalize successfully.
    type(table_t) :: raw

    ok = .false.
    if (.not. valid_table(cpt)) return
    if (size(values) /= size(cpt%value)) return
    if (any(values < 0.0_dp)) return
    raw = make_table(cpt%var, cpt%dim, values)
    cpt = table_normalize_first(raw)
    ok = valid_table(cpt)
  end subroutine replace_cpt_values

  pure function genotype_pair(state, n_alleles) result(pair)
    integer, value :: state !! One-based unordered genotype index to decode.
    integer, value :: n_alleles !! Number of allele labels in the genotype state space.
    integer :: pair(2)
    integer :: a
    integer :: b
    integer :: k

    pair = 0
    if (n_alleles < 1 .or. state < 1) return
    k = 0
    do a = 1, n_alleles
      do b = a, n_alleles
        k = k + 1
        if (k == state) then
          pair = [a, b]
          return
        end if
      end do
    end do
  end function genotype_pair

  pure real(dp) function transmission_probability(allele, genotype) result(probability)
    integer, value :: allele !! Allele label whose transmission probability is requested.
    integer, intent(in) :: genotype(2) !! Two allele labels defining an unordered diploid genotype.

    probability = 0.5_dp * real(count(genotype == allele), dp)
  end function transmission_probability

end module grain_cpt
