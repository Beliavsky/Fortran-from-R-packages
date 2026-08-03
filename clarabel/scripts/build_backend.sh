#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
bridge="$root/rust_bridge"
bindir="$bridge/bin"

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: Cargo/Rust is required to build the Clarabel backend" >&2
  exit 1
fi
if [ ! -d "$bridge/vendor/clarabel" ]; then
  echo "Extracting the vendored Rust dependency tree..."
  tar -xJf "$bridge/vendor.tar.xz" -C "$bridge"
fi
mkdir -p "$bridge/.cargo" "$bindir"
cat > "$bridge/.cargo/config.toml" <<'CFG'
[source.crates-io]
replace-with = "vendored-sources"
[source.vendored-sources]
directory = "vendor"
[net]
offline = true
CFG
cargo build --manifest-path "$bridge/Cargo.toml" --release --offline --locked

case "$(uname -s)" in
  Darwin)
    library="$bridge/target/release/libclarabel_fortran_bridge.dylib"
    ;;
  Linux|FreeBSD|NetBSD|OpenBSD)
    library="$bridge/target/release/libclarabel_fortran_bridge.so"
    ;;
  *)
    echo "error: unsupported Unix platform $(uname -s)" >&2
    exit 1
    ;;
esac
if [ ! -f "$library" ]; then
  echo "error: Cargo completed but the shared bridge library was not found: $library" >&2
  exit 1
fi
cp -f "$library" "$bindir/"
echo "Built Clarabel backend in $bindir"
echo "A normal fpm build no longer needs LIBRARY_PATH, LD_LIBRARY_PATH, or -L flags."
