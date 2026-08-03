# API overview

All integer values use `integer(ik)`, where `ik` is `int64`.

## Magic squares

- `magic_square_of_order(n)` selects the appropriate odd, doubly-even, or
  singly-even construction.
- `magic_2np1(m)`, `magic_4n(m)`, `magic_4np2(m)`, `strachey_square(m)`, and
  `lozenge_square(m)` expose individual construction families.
- `hudson_square(n)`, `magic_prime(n)`, `panmagic_4n(m)`,
  `panmagic_6np1(m)`, and `panmagic_6nm1(m)` construct special squares.
- `magic_product_fast(a,b)` forms the compound product of two squares.
- `is_magic`, `is_semimagic`, `is_panmagic`, `is_normal_square`,
  `is_mostperfect`, `is_associative`, and symmetry predicates inspect them.
- `is_antimagic`, `is_totally_antimagic`, `is_heterosquare`, `is_sam`, and
  `is_stam` implement the corresponding line-sum tests.
- `is_multiplicative_magic` is the product analogue of `is_magic`.

## Tensor operations

`integer_tensor` stores `shape(:)` and column-major `values(:)`. Its methods
are `rank`, `size`, `get`, `set`, and `valid`.

Constructors and transformations include `make_tensor`, `sequence_tensor`,
`tensor_permute`, `tensor_reverse`, `tensor_shift`, `tensor_rotate`,
`tensor_block_diag`, `tensor_overlay_add`, `tensor_subsums`, `tensor_pad`,
`tensor_take`, and `tensor_drop`.

## Hypercubes

- `magiccube_2np1(m)` constructs odd-order magic cubes.
- `magichypercube_4n(m,dimension)` constructs order `4m` hypercubes.
- `is_semimagichypercube`, `is_diagonally_correct`, `is_magichypercube`,
  `is_latinhypercube`, `is_perfect_hypercube`, and `is_alicehypercube`
  implement the upstream invariant tests.
- `diagonal_subhypercubes` exposes the subhypercubes used by the perfectness
  test.

## Combinatorial constructions

- `latin_square`, `incidence`, `unincidence`, `incidence_move`,
  `another_incidence`, `another_latin`, and `random_latin_squares`
- `sylvester_hadamard` and `is_hadamard`
- `sam_square`
- `cilleruelo_square`
- `bernhardsson_a`, `bernhardsson_b`, and `bernhardsson_matrix`

Procedures that can reject an input accept an optional `magic_error` argument.
