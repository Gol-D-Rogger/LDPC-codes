#!/bin/bash
# LSF Shim Installation Script

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== LSF Shim Installation ==="
echo ""
echo "Installation directory: $SCRIPT_DIR"
echo ""

# Make scripts executable
echo "Making scripts executable..."
chmod +x "$SCRIPT_DIR/bsub"
chmod +x "$SCRIPT_DIR/bjobs"
chmod +x "$SCRIPT_DIR/bkill"
chmod +x "$SCRIPT_DIR/lsf_runner.py"

echo "✓ Scripts are now executable"
echo ""

# Check Python version
echo "Checking Python version..."
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "Found Python $PYTHON_VERSION"

# Verify Python 3.6+
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 6 ]); then
    echo "⚠ Warning: Python 3.6+ is required. Found Python $PYTHON_VERSION"
    exit 1
fi

echo "✓ Python version is compatible"
echo ""

# Create default database directory
DEFAULT_DB_DIR="$HOME/.lsf_shim"
if [ ! -d "$DEFAULT_DB_DIR" ]; then
    echo "Creating default database directory: $DEFAULT_DB_DIR"
    mkdir -p "$DEFAULT_DB_DIR"
    echo "✓ Database directory created"
else
    echo "✓ Database directory already exists: $DEFAULT_DB_DIR"
fi
echo ""

# Installation options
echo "=== Installation Complete ==="
echo ""
echo "To use the LSF shim, add it to your PATH:"
echo ""
echo "  export PATH=\"$SCRIPT_DIR:\$PATH\""
echo ""
echo "Add this line to your ~/.bashrc or ~/.zshrc to make it permanent."
echo ""
echo "Optional: Create symlinks to ~/bin (if ~/bin is in your PATH):"
echo ""
echo "  mkdir -p ~/bin"
echo "  ln -sf $SCRIPT_DIR/bsub ~/bin/bsub"
echo "  ln -sf $SCRIPT_DIR/bjobs ~/bin/bjobs"
echo "  ln -sf $SCRIPT_DIR/bkill ~/bin/bkill"
echo ""
echo "=== Configuration ==="
echo ""
echo "Environment variables (optional):"
echo "  LSF_SHIM_DB        - Path to job database (default: ~/.lsf_shim/jobs.json)"
echo "  LSF_SHIM_PEND_SEC  - PEND duration in seconds (default: 10)"
echo ""
echo "Example:"
echo "  export LSF_SHIM_DB=/tmp/my_project/lsf_jobs.json"
echo "  export LSF_SHIM_PEND_SEC=5"
echo ""
