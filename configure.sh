#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# configure.sh  —  Bootstrap and configure dropflow2D
#
# Usage:
#   ./configure.sh [configure options]
#
# Examples:
#   ./configure.sh                         standard optimised build
#   ./configure.sh --enable-debug          debug build (-g -O0)
#   ./configure.sh --prefix=$HOME/.local   install to home directory
#
# What this script does:
#   1. Checks that the required autotools are present (autoconf, automake).
#   2. Runs 'autoreconf -fi' to regenerate configure, Makefile.in, etc.
#   3. Runs './configure' with any arguments you passed in.
#   4. If autotools are missing, falls back to the standalone Makefile.
#
# After this script completes, just run:
#   make
#   make install   (optional)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Colours (disabled if not a terminal) ─────────────────────────────────────
if [ -t 1 ]; then
    RED='\033[0;31m'; YELLOW='\033[1;33m'
    GREEN='\033[0;32m'; BOLD='\033[1m'; RESET='\033[0m'
else
    RED=''; YELLOW=''; GREEN=''; BOLD=''; RESET=''
fi

info()    { echo -e "${BOLD}[configure.sh]${RESET} $*"; }
success() { echo -e "${GREEN}[configure.sh] $*${RESET}"; }
warn()    { echo -e "${YELLOW}[configure.sh] WARNING: $*${RESET}"; }
error()   { echo -e "${RED}[configure.sh] ERROR: $*${RESET}" >&2; exit 1; }

# ── Must be run from the project root ────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ ! -f configure.ac ]; then
    error "configure.ac not found. Run this script from the dropflow2D source root."
fi

# ── Check for autotools ───────────────────────────────────────────────────────
HAVE_AUTOTOOLS=true
for tool in autoconf automake autoreconf aclocal; do
    if ! command -v "$tool" &>/dev/null; then
        warn "$tool not found."
        HAVE_AUTOTOOLS=false
    fi
done

# ─────────────────────────────────────────────────────────────────────────────
# Autotools path
# ─────────────────────────────────────────────────────────────────────────────
if $HAVE_AUTOTOOLS; then

    info "Autotools found. Bootstrapping..."

    # Remove stale generated files so autoreconf starts clean.
    info "Cleaning stale autotools artefacts..."
    rm -f configure Makefile.in aclocal.m4
    rm -rf autom4te.cache

    # Generate configure, Makefile.in, config.h.in, etc.
    info "Running autoreconf -fi ..."
    autoreconf -fi

    success "Bootstrap complete."

    # Run configure, forwarding all command-line arguments.
    info "Running ./configure $* ..."
    echo ""
    ./configure "$@"

    echo ""
    success "Configuration complete."
    echo ""
    echo -e "  ${BOLD}Next steps:${RESET}"
    echo "    make              — build the executable"
    echo "    make install      — install to prefix  [optional]"
    echo "    make vtkclean     — remove VTK output directories"
    echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Fallback: use the standalone Makefile directly
# ─────────────────────────────────────────────────────────────────────────────
else
    warn "One or more autotools programs were not found."
    warn "Falling back to the standalone Makefile."
    echo ""

    if [ ! -f Makefile ]; then
        error "Makefile not found either. Cannot build."
    fi

    # Parse --enable-debug from the forwarded arguments so the fallback
    # path still honours the same flag.
    DEBUG_FLAG=""
    for arg in "$@"; do
        if [ "$arg" = "--enable-debug" ]; then
            DEBUG_FLAG='OPT="-g -O0"'
            warn "--enable-debug detected: setting OPT=-g -O0"
        fi
    done

    info "Standalone Makefile is ready."
    echo ""
    echo -e "  ${BOLD}To build, run:${RESET}"
    if [ -n "$DEBUG_FLAG" ]; then
        echo "    make $DEBUG_FLAG"
    else
        echo "    make"
    fi
    echo ""
    echo "  To install autotools on common platforms:"
    echo "    Ubuntu/Debian:  sudo apt-get install autoconf automake"
    echo "    Fedora/RHEL:    sudo dnf install autoconf automake"
    echo "    macOS (brew):   brew install autoconf automake"
    echo ""
fi
