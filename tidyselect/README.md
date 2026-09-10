# tidyselect

Restricted modern Fortran translation of **tidyselect 1.2.1**, *Select from a
Set of Strings*. Selectors accept available names and return one-based integer
positions. Implemented operations include `everything`, `all_of`, `any_of`,
`one_of`, `starts_with`, `ends_with`, `contains`, fixed-pattern `matches`,
`last_col`, `num_range`, selection combination, relocation, pulling, and
renaming. R quosures, data masks, and nonstandard evaluation are excluded.

The returned indices can be passed directly to the translated `dplyr`
`select`, `rename`, and `relocate` overloads and to `tidyr::pivot_longer`.

```text
fpm build
fpm test
```
