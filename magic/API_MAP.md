# Upstream-to-Fortran API map

## Direct or close equivalents

| R function | Fortran interface |
|---|---|
| `adiag()` | repeated `tensor_block_diag()` |
| `aplus()` | `tensor_overlay_add()` |
| `arev()` | `tensor_reverse()` / `reverse_square()` |
| `arot()` | `tensor_rotate()` / `rotate_square()` |
| `ashift()` | `tensor_shift()` |
| `arow()` | `tensor_axis_coordinates()` |
| `circulant()`, `latin()` | `circulant()`, `latin_square()` |
| `diag.off()`, `allsums()` | `diag_off()`, `all_square_sums()` |
| `magic()` | `magic_square_of_order()` |
| `magic.2np1()` | `magic_2np1()` |
| `magic.4n()` | `magic_4n()` |
| `magic.4np2()`, `strachey()` | `magic_4np2()`, `strachey_square()` |
| `lozenge()` | `lozenge_square()` |
| `hudson()` | `hudson_square()` |
| `magic.prime()` | `magic_prime()` |
| `magic.product.fast()` | `magic_product_fast()` |
| `panmagic.4n()` | `panmagic_4n()` |
| `panmagic.6npm1()` | `panmagic_6npm1()` |
| `magiccube.2np1()` | `magiccube_2np1()` |
| `magichypercube.4n()` | `magichypercube_4n()` |
| `subsums()` | `tensor_subsums()` |
| `apad()` | `tensor_pad()` |
| `apltake()`, `apldrop()` | `tensor_take()`, `tensor_drop()` |
| `fnsd()` | `first_nonsingleton_dimensions()` |
| `transf()` | `transform_square()` |
| `as.standard()` | `standardize_square()` |
| `is.magic()` | `is_magic()` |
| `is.semimagic()` | `is_semimagic()` |
| `is.panmagic()` | `is_panmagic()` |
| `is.normal()` | `is_normal_square()` |
| `is.associative()` | `is_associative()` |
| `is.mostperfect()` | `is_mostperfect()` |
| `is.magichypercube()` | `is_magichypercube()` |
| `is.perfect()` | `is_perfect_hypercube()` |
| `allsubhypercubes()` | `diagonal_subhypercubes()` |
| `incidence()`, `unincidence()` | same names |
| `inc_to_inc()` | `incidence_move()` |
| `another_latin()`, `rlatin()` | `another_latin()`, `random_latin_squares()` |
| `sylvester()` | `sylvester_hadamard()` |
| `sam()` | `sam_square()` |
| `cilleruelo()` | `cilleruelo_square()` |
| `bernhardsson*()` | `bernhardsson_*()` |

## Adapted interfaces

R's optional `func` callback in magicness tests is represented by separate
sum-based and product-based predicates. `is_multiplicative_magic` covers the
upstream `func=prod` use case. Arbitrary user callbacks are not stored inside
Fortran tensor values.

`standardize_square` chooses the lexicographically smallest dihedral image.
This is a deterministic canonical form and agrees with Frenicle-style
standardization for the normal squares exercised by the test suite, although
its definition is slightly stronger than the upstream corner-and-adjacency
rule.

`lozenge_square` uses the equivalent odd-order Siamese construction rather
than reproducing the upstream coordinate placement order exactly.

## Intentionally omitted

- `magicplot()` and all graphics
- R S3/S4 dispatch, lists, `dimnames`, replacement syntax, and operators such
  as `%eq%`
- `do.index()` with an arbitrary R closure
- packaged `.rda` data as compiled Fortran constants
- `magic.8()` exhaustive enumeration and the customizable binary-carpet
  forms of `panmagic.4()` and `panmagic.8()`
- textual APL replacement functions and R-specific coercion helpers
