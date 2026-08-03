# Building Clarabel Fortran

## Why the Rust bridge is runtime-loaded

FPM can declare an external library with `build.link`, but it does not build a
Rust crate or add a package-relative library search directory. The original
manifest therefore emitted `-lclarabel_fortran_bridge` and failed during
executable linking unless users manually configured `-L`/`LIBRARY_PATH`.

The package now compiles `src/clarabel_dynamic_loader.c` as part of the FPM
library. It provides the C ABI symbols expected by the Fortran modules and
loads the real Clarabel bridge shared library only when needed. Consequently:

- `fpm build` works without Cargo and without linker search-path variables.
- `scripts/build_backend.bat` or `.sh` builds the production Rust backend.
- `CLARABEL_FORTRAN_BRIDGE` can select an absolute backend library path.
- A missing backend is reported at solver creation instead of link time.
- The default demo gives setup instructions and exits normally when no backend is installed.

## Windows

```bat
scripts\build_with_backend.bat run
```

The wrapper calls `build_backend.bat`, then runs the requested FPM command with
`CLARABEL_FORTRAN_BRIDGE` set to the absolute DLL path. The backend script builds
`clarabel_fortran_bridge.dll`, copies it to `rust_bridge\bin`, copies common
MinGW runtime DLLs when required, and verifies a real QP solve through the
Fortran interface.

Because the bridge is runtime-loaded, a MinGW-built Fortran executable can use
a DLL produced by either the GNU or MSVC Rust host toolchain; no MinGW import
library is required. After the first successful backend build, plain `fpm run`
from the package root finds `rust_bridge\bin\clarabel_fortran_bridge.dll`.

One-command wrappers are also available:

```bat
scripts\build_with_backend.bat build
scripts\build_with_backend.bat run
scripts\build_with_backend.bat test
```

## GNU/Linux and macOS

```sh
./scripts/build_backend.sh
fpm build
fpm run
```

The loader uses `dlopen` and checks `rust_bridge/bin` in addition to the normal
shared-library search path.

## Windows loader troubleshooting

If the DLL file exists but Windows returns error 126 ("The specified module could not be found"), the missing module is usually a dependency of the DLL rather than the bridge itself. Re-run:

```bat
scripts\build_backend.bat
```

The script now copies common GNU runtime DLLs and verifies the bridge. To list remaining dependencies manually:

```bat
objdump -p rust_bridge\bin\clarabel_fortran_bridge.dll | findstr "DLL Name"
```
