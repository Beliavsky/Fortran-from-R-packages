# Validation

Validation was performed with GNU Fortran 14.2 and GCC 14.2 against the C
sources supplied in `rngWELL-master.zip`.

## Strict Fortran build

The library, permanent tests, and example compile with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O2
```

Both permanent FPM-style tests pass:

- `test_sequences.f90`: reference streams for all 17 variants, state
  round-trip, compatibility option mapping, and matrix fill order.
- `test_properties.f90`: all variants, output range, loose distribution smoke
  checks, and MT2002 seed expansion.

## Bit-for-bit C differential test

A separate development harness compiled the original WELL C files and the
Fortran translation independently. For each of the 17 variants, 4,000 raw
32-bit outputs were generated from each of four seeds:

```text
0
1
0xffffffff
123456789
```

That is 272,000 generated 32-bit words. The concatenated C and Fortran binary
files were byte-identical, with SHA-256:

```text
b85d7a8d646d3d1ce026f48ee5c7fe3c357e5e4f6bd0306937279f57044083cd
```

This length crosses every circular-state boundary multiple times, including
the 1,391-word WELL44497 state.

## State differential test

For every variant, the upstream C generator and Fortran generator were seeded
with 42, advanced by 1,733 outputs, and their canonical exported states were
compared byte-for-byte. All state images matched. The concatenated state image
SHA-256 was:

```text
0055a11fd9c81486e85bb445cc9ac0c347ee1e5a7076cfaf71b3d304779c1867
```
