# ─────────────────────────────────────────────────────────────────────────────
# Makefile.am  —  Automake input for dropflow2D
#
# Source files live in src/.  The 'subdir-objects' option in configure.ac
# tells automake to place .o files alongside their .cpp (src/*.o), which
# avoids name collisions and keeps the build root clean.
# ─────────────────────────────────────────────────────────────────────────────

# ── Executable ───────────────────────────────────────────────────────────────
bin_PROGRAMS = dropflow2D

# ── Sources ───────────────────────────────────────────────────────────────────
dropflow2D_SOURCES = \
    src/dropflow.cpp \
    src/field.cpp    \
    src/boundary.cpp \
    src/init.cpp     \
    src/twophase.cpp \
    src/velocity.cpp \
    src/timeadv.cpp  \
    src/output.cpp

# Headers are listed so 'make dist' includes them in release tarballs.
noinst_HEADERS = \
    src/field.h    \
    src/boundary.h \
    src/init.h     \
    src/twophase.h \
    src/velocity.h \
    src/timeadv.h  \
    src/output.h

# ── Extra files included in 'make dist' tarballs ─────────────────────────────
EXTRA_DIST = \
    configure.ac \
    Makefile.am  \
    configure.sh

# ── Remove generated VTK output directories ──────────────────────────────────
.PHONY: vtkclean
vtkclean:
	rm -rf time_step-* *.pvtr
	@echo "VTK output files removed."
