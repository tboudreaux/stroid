#!/bin/bash

set -e

# Confirm that the script is being run on a macOS aarch64 host
ARCHITECTURE=$(uname -m)
OS_NAME=$(uname -s)
if [[ "$OS_NAME" != "Darwin" || "$ARCHITECTURE" != "arm64" ]]; then
    echo "This release script must be run on a macOS aarch64 host."
    exit 1
fi

PROJECT_NAME="stroid"
MACOS_MIN_VERSION="15.0"
OUTPUT_DIR="releases"

# Confirm that Meson and Ninja are installed
if ! command -v meson &> /dev/null || ! command -v ninja &> /dev/null; then
    echo "Meson and Ninja are required to build the project."
    echo "Please install them via pip:"
    echo "  pip install meson ninja"
    exit 1
fi

# Confirm that the current macOS version meets the minimum requirement
CURRENT_MACOS_VERSION=$(sw_vers -productVersion)
if [[ "$(printf '%s\n%s\n' "$MACOS_MIN_VERSION" "$CURRENT_MACOS_VERSION" | sort -V | head -n1)" != "$MACOS_MIN_VERSION" ]]; then
    echo "Current macOS version ($CURRENT_MACOS_VERSION) does not meet the minimum requirement ($MACOS_MIN_VERSION)."
    exit 1
fi


X86_IMAGE="quay.io/pypa/manylinux_2_28_x86_64"
ARM_IMAGE="quay.io/pypa/manylinux_2_28_aarch64"

mkdir -p "$OUTPUT_DIR"

echo "Starting release build for $PROJECT_NAME..."

echo "Building for macOS (ARM64, Target >= $MACOS_MIN_VERSION)..."
rm -rf build-mac
export MACOSX_DEPLOYMENT_TARGET=$MACOS_MIN_VERSION

meson setup build-mac \
    --buildtype release \
    -Dcpp_args="-mmacosx-version-min=$MACOS_MIN_VERSION" \
    -Dcpp_link_args="-mmacosx-version-min=$MACOS_MIN_VERSION"

meson compile -C build-mac
cp build-mac/tools/stroid "$OUTPUT_DIR/stroid-macos-arm64"

echo "Building for Linux (x86_64) via Manylinux 2.28..."
docker run --rm --platform linux/amd64 \
    -v "$(pwd)":/src -w /src \
    $X86_IMAGE /bin/bash -c "
        /opt/python/cp311-cp311/bin/pip install meson ninja
        export PATH=/opt/python/cp311-cp311/bin:\$PATH

        rm -rf build-linux-x86
        meson setup build-linux-x86 --buildtype release -Dbuild_tests=false
        meson compile -C build-linux-x86
        cp build-linux-x86/tools/stroid /src/$OUTPUT_DIR/stroid-linux-x86_64
    "

echo "Building for Linux (ARM64) via Manylinux 2.28..."
docker run --rm --platform linux/arm64 \
    -v "$(pwd)":/src -w /src \
    $ARM_IMAGE /bin/bash -c "
        /opt/python/cp311-cp311/bin/pip install meson ninja
        export PATH=/opt/python/cp311-cp311/bin:\$PATH

        rm -rf build-linux-arm
        meson setup build-linux-arm --buildtype release -Dbuild_tests=false
        meson compile -C build-linux-arm
        cp build-linux-arm/tools/stroid /src/$OUTPUT_DIR/stroid-linux-arm64
    "

echo "---"
echo "All binaries generated in /$OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"