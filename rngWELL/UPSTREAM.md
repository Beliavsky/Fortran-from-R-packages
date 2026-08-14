# Upstream provenance

Translated from the supplied `rngWELL` source package:

- Package: rngWELL
- Version: 0.10-10
- Package authors: Christophe Dutang and Petr Savicky
- WELL C code contributors: Francois Panneton, Pierre L'Ecuyer, Makoto Matsumoto
- CRAN publication date in supplied DESCRIPTION: 2024-10-17

The WELL C files identify their algorithmic source as the WELL generators of
Panneton, L'Ecuyer, and Matsumoto and reference:

F. Panneton, P. L'Ecuyer, and M. Matsumoto (2006),
"Improved Long-Period Generators Based on Linear Recurrences Modulo 2",
ACM Transactions on Mathematical Software.

R registration, `.Call`/`.C` plumbing, configure scripts, locale wrappers, and
package-version reporting are not computational RNG algorithms and were not
ported. There is no plotting code in this package.
