# Translation coverage

The port covers the package's central computational themes:

- arbitrary-rank integer-array manipulation
- normal, panmagic, associative, most-perfect, antimagic, sparse antimagic,
  and multiplicative magic squares
- odd, doubly-even, singly-even, Hudson, prime, compound, and panmagic
  construction families
- odd magic cubes and order-4n magic hypercubes
- semimagic, diagonal, perfect, Latin, and Alice hypercube tests
- Latin-square incidence tensors and random incidence moves
- Sylvester Hadamard and Bernhardsson constructions

The exclusions are presentation or R-language infrastructure, except for the
explicitly documented exhaustive/customizable generators in `API_MAP.md`.
The retained upstream snapshot allows omitted interfaces to be audited or
ported later.
