# Third-party licenses

`rust_bridge/vendor.tar.xz` is the dependency archive distributed with the attached R package. It contains Clarabel.rs and all Cargo dependencies needed for an offline build. Each vendored crate retains its own `LICENSE*`, `COPYING*`, or package metadata. Extract the archive and inspect `rust_bridge/vendor/<crate>/` for the controlling text.

The top-level project does not relicense those dependencies. The Apache-2.0 top-level license applies only to the new Fortran interface and Rust bridge, while Clarabel.rs and each dependency remain under their stated licenses.
