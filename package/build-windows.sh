#!/bin/bash
# Build portable Windows exe (no installer required)
#
# Prerequisites:
#   rustup target add x86_64-pc-windows-gnu
#   apt install gcc-mingw-w64-x86-64   # on Debian/Ubuntu
#
# Output: target/x86_64-pc-windows-gnu/release/fastqc-compliant-rs-gui.exe
#         package/dist/fastqc-compliant-rs-windows/

set -euo pipefail
cd "$(dirname "$0")/.."

TARGET="x86_64-pc-windows-gnu"
DIST_DIR="package/dist/fastqc-compliant-rs-windows"

echo "=== Building fastqc-compliant-rs for Windows (portable exe) ==="

# Build CLI
RUSTFLAGS="-C target-cpu=native" \
    cargo build --release --target "$TARGET" --bin fastqc-compliant-rs

# Build GUI
RUSTFLAGS="-C target-cpu=native" \
    cargo build --release --target "$TARGET" --features gui --bin fastqc-compliant-rs-gui

echo "=== Packaging ==="
rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR"

cp "target/$TARGET/release/fastqc-compliant-rs.exe" "$DIST_DIR/"
cp "target/$TARGET/release/fastqc-compliant-rs-gui.exe" "$DIST_DIR/"
cp README.md "$DIST_DIR/"
cp LICENSE "$DIST_DIR/" 2>/dev/null || cp fastqc/LICENSE "$DIST_DIR/"

echo "=== Done ==="
echo "Portable distribution: $DIST_DIR/"
echo "  fastqc-compliant-rs.exe      — CLI tool"
echo "  fastqc-compliant-rs-gui.exe  — GUI tool (double-click to run)"
ls -lh "$DIST_DIR/"
