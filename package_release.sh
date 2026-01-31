#!/bin/bash

set -e

RELEASE_DIR="releases"
STAGING_DIR="stroid-dist"

REQUIRED_BINARIES=(
    "stroid-macos-arm64"
    "stroid-linux-x86_64"
    "stroid-linux-arm64"
)

for bin in "${REQUIRED_BINARIES[@]}"; do
    if [ ! -f "$RELEASE_DIR/$bin" ]; then
        echo "Error: Missing binary $RELEASE_DIR/$bin. Run release.sh first."
        exit 1
    fi
done

chmod +x "$RELEASE_DIR/stroid-macos-arm64"
VERSION_OUTPUT=$("$RELEASE_DIR/stroid-macos-arm64" info -v)
VERSION=$(echo "$VERSION_OUTPUT" | awk '{print $3}')

if [ -z "$VERSION" ]; then
    echo "Error: Could not extract version from binary output."
    exit 1
fi

TARBALL_NAME="stroid-${VERSION}.tar.gz"

echo "Packaging Stroid version ${VERSION}..."

rm -rf "$STAGING_DIR"
mkdir -p "$STAGING_DIR/bin"

for bin in "${REQUIRED_BINARIES[@]}"; do
    cp "$RELEASE_DIR/$bin" "$STAGING_DIR/bin/"
done

cat << 'EOF' > "$STAGING_DIR/install.sh"
#!/bin/bash
set -e

OS=$(uname -s)
ARCH=$(uname -m)
INSTALL_PATH=$1

echo "Detecting system: $OS ($ARCH)"

case "$OS" in
    Darwin)
        if [ "$ARCH" = "arm64" ]; then
            SELECTED_BIN="bin/stroid-macos-arm64"
        else
            echo "Error: Intel Macs (x86_64) are not supported by this package."
            exit 1
        fi
        ;;
    Linux)
        if [ "$ARCH" = "x86_64" ]; then
            SELECTED_BIN="bin/stroid-linux-x86_64"
        elif [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then
            SELECTED_BIN="bin/stroid-linux-arm64"
        else
            echo "Error: Linux architecture $ARCH is not supported."
            exit 1
        fi
        ;;
    *)
        echo "Error: Operating system $OS is not supported."
        exit 1
        ;;
esac

echo "Verifying binary integrity..."
chmod +x "$SELECTED_BIN"
if ! ./"$SELECTED_BIN" info -v > /dev/null 2>&1; then
    echo "Error: The selected binary failed the version check. It may be incompatible with this system's libraries."
    exit 1
fi

VERSION_INFO=$(./"$SELECTED_BIN" info -v)
echo "Verification successful: $VERSION_INFO"

if [ -z "$INSTALL_PATH" ]; then
    if [ "$EUID" -eq 0 ]; then
        INSTALL_PATH="/usr/local/bin"
    else
        INSTALL_PATH="$HOME/.local/bin"
    fi
fi

mkdir -p "$INSTALL_PATH"

echo "Installing to: $INSTALL_PATH/stroid"
cp "$SELECTED_BIN" "$INSTALL_PATH/stroid"
chmod +x "$INSTALL_PATH/stroid"

echo "Installation complete."
if [[ ":$PATH:" != *":$INSTALL_PATH:"* ]]; then
    echo "Note: $INSTALL_PATH is not in your PATH. Add it to your shell profile to run 'stroid' from anywhere."
fi
EOF

chmod +x "$STAGING_DIR/install.sh"

export COPYFILE_DISABLE=1
tar -cz --no-xattrs -f "$TARBALL_NAME" -C "$STAGING_DIR" .

echo "Distribution package created: $TARBALL_NAME"

rm -rf "$STAGING_DIR"