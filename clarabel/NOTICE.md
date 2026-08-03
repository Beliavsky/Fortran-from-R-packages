# Notices and provenance

- Source package translated: `clarabel` R package 0.11.2.
- Solver dependency represented by the bundled source: Clarabel.rs 0.11.1.
- Original uploaded archive: `upstream/clarabel-r-main.zip`.
- Original Rust dependency archive: `rust_bridge/vendor.tar.xz`.
- The Fortran source and the new Rust C-ABI bridge are licensed under Apache-2.0.
- The original R package and Clarabel.rs notices are retained in `LICENSE`, `LICENSE.note`, `docs/UPSTREAM_AUTHORS`, and `docs/UPSTREAM_CITATION`.
- Third-party crates inside the vendor archive retain their individual license files and Cargo checksum metadata.

No R plotting, S3/S4 class machinery, `Matrix` dispatch, package-site files, or R printing methods are translated into Fortran.

## FPM linker-path correction

This revision replaces the manifest-level `-lclarabel_fortran_bridge` link with
a small runtime loader compiled by FPM. Plain `fpm build` therefore succeeds
without a prebuilt Rust library or manual `-L` configuration. Actual solves
still use the unmodified Clarabel.rs backend built by the supplied scripts.
