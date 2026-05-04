#!/bin/bash
# Build macOS .app bundle
#
# Prerequisites:
#   On macOS: just run this script
#   Cross from Linux: rustup target add x86_64-apple-darwin / aarch64-apple-darwin
#
# Output: package/dist/FastQC.app/

set -euo pipefail
cd "$(dirname "$0")/.."

# Detect architecture
if [[ "$(uname -m)" == "arm64" ]]; then
    TARGET="aarch64-apple-darwin"
else
    TARGET="x86_64-apple-darwin"
fi

APP_DIR="package/dist/FastQC.app"
CONTENTS="$APP_DIR/Contents"
MACOS_DIR="$CONTENTS/MacOS"
RESOURCES="$CONTENTS/Resources"

echo "=== Building fastqc-compliant-rs for macOS ($TARGET) ==="

# Build CLI
RUSTFLAGS="-C target-cpu=native" \
    cargo build --release --target "$TARGET" --bin fastqc-compliant-rs 2>/dev/null || \
    cargo build --release --bin fastqc-compliant-rs

# Build GUI
RUSTFLAGS="-C target-cpu=native" \
    cargo build --release --target "$TARGET" --features gui --bin fastqc-compliant-rs-gui 2>/dev/null || \
    cargo build --release --features gui --bin fastqc-compliant-rs-gui

echo "=== Creating .app bundle ==="
rm -rf "$APP_DIR"
mkdir -p "$MACOS_DIR" "$RESOURCES"

# Copy binaries
if [ -f "target/$TARGET/release/fastqc-compliant-rs-gui" ]; then
    cp "target/$TARGET/release/fastqc-compliant-rs-gui" "$MACOS_DIR/FastQC"
    cp "target/$TARGET/release/fastqc-compliant-rs" "$MACOS_DIR/fastqc-compliant-rs"
else
    cp "target/release/fastqc-compliant-rs-gui" "$MACOS_DIR/FastQC"
    cp "target/release/fastqc-compliant-rs" "$MACOS_DIR/fastqc-compliant-rs"
fi

# Create Info.plist
cat > "$CONTENTS/Info.plist" << 'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key>
    <string>FastQC</string>
    <key>CFBundleDisplayName</key>
    <string>FastQC</string>
    <key>CFBundleIdentifier</key>
    <string>uk.ac.babraham.fastqc-compliant-rs</string>
    <key>CFBundleVersion</key>
    <string>0.4.2</string>
    <key>CFBundleShortVersionString</key>
    <string>0.4.2</string>
    <key>CFBundleExecutable</key>
    <string>FastQC</string>
    <key>CFBundleIconFile</key>
    <string>fastqc</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>NSHighResolutionCapable</key>
    <true/>
    <key>CFBundleDocumentTypes</key>
    <array>
        <dict>
            <key>CFBundleTypeExtensions</key>
            <array>
                <string>fastq</string>
                <string>fq</string>
                <string>fastq.gz</string>
                <string>bam</string>
                <string>sam</string>
            </array>
            <key>CFBundleTypeName</key>
            <string>Sequence File</string>
            <key>CFBundleTypeRole</key>
            <string>Viewer</string>
        </dict>
    </array>
</dict>
</plist>
PLIST

# Convert icon to icns if sips is available (macOS only)
if command -v sips &>/dev/null && command -v iconutil &>/dev/null; then
    ICONSET="$RESOURCES/fastqc.iconset"
    mkdir -p "$ICONSET"
    # Use the PNG icon as source
    SRC_ICON="src/resources/icons/fastqc_icon.png"
    for SIZE in 16 32 64 128 256 512; do
        sips -z $SIZE $SIZE "$SRC_ICON" --out "$ICONSET/icon_${SIZE}x${SIZE}.png" 2>/dev/null || true
        DOUBLE=$((SIZE * 2))
        sips -z $DOUBLE $DOUBLE "$SRC_ICON" --out "$ICONSET/icon_${SIZE}x${SIZE}@2x.png" 2>/dev/null || true
    done
    iconutil -c icns "$ICONSET" -o "$RESOURCES/fastqc.icns" 2>/dev/null || true
    rm -rf "$ICONSET"
else
    echo "Note: sips/iconutil not available — .app will use default icon"
    # Copy PNG as fallback
    cp src/resources/icons/fastqc_icon.png "$RESOURCES/fastqc.png"
fi

echo "=== Done ==="
echo "App bundle: $APP_DIR"
echo "  Double-click to open, or: open $APP_DIR"
du -sh "$APP_DIR"
