# Licensing

## ewens-derived translation code

The files in the top-level `src/`, `test/`, and `example/` directories are a
modern Fortran translation of **ewens 0.1.0** by Chris Hanretty. The upstream
package is distributed under the MIT license (`MIT + file LICENSE`). The
original copyright notice is:

> Copyright (c) 2026 Chris Hanretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Vendored copula-fortran dependency

`vendor/copula-fortran` is the supplied Fortran translation of the R package
`copula` and is licensed **GPL-3.0-or-later**. Its own `LICENSE`, `NOTICE`, and
provenance files are retained unchanged.

Because `ewens-fortran` links against that GPL-3.0-or-later dependency, a
linked/distributed combined binary must be distributed under terms compatible
with GPL-3.0-or-later. The MIT terms above continue to apply to the
`ewens`-derived source files themselves.
