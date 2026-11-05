#!/bin/bash
# build_patched_staar.sh -- build the root-checked STAAR and STAARpipeline into the NEW library
#             ${R_LIBS_ROOT}/staar_rootcheck_4.4 (the default library is never touched).
#
# Steps: (1) fresh copies of the two untrusted source trees (exact upstream commits) are made under
# build/; (2) patches/STAAR_saddlepoint_root_check.patch and patches/STAARpipeline_saddlepoint_root_check.patch are applied
# (patch -p1; the build stops if a hunk fails); (3) R CMD INSTALL --library=<new lib>; (4) the
# installed DESCRIPTIONs are checked for the RootcheckPatch marker and the record is written to
# build/BUILD_RECORD.txt (versions, upstream SHAs, tarball sha256, patch sha256, sessionInfo).
# Idempotent: re-running rebuilds from scratch. Usage: bash build_patched_staar.sh   (login shell or jobs/01_build.sbatch)
set -euo pipefail
PKG_DIR=${ADSP_ROOT}/revision/ukb_rootcheck_20261019
source "$PKG_DIR/env_rootcheck.sh"
cd "$PKG_DIR"

STAAR_SHA=9db9dd504905b9f469146f670e5f6dbe3e08d01a
PIPE_SHA=7e868f23c85ba0713dd00c8512ad23cca7fbc8f9
STAAR_SRC="$PKG_DIR/src_untrusted/STAAR-$STAAR_SHA"
PIPE_SRC="$PKG_DIR/src_untrusted/STAARpipeline-$PIPE_SHA"
[ -d "$STAAR_SRC" ] || { echo "missing $STAAR_SRC (run the curl commands in README.md first)"; exit 2; }
[ -d "$PIPE_SRC" ]  || { echo "missing $PIPE_SRC"; exit 2; }

# the two installed packages must still be the ones the patch was written against
for p in STAAR STAARpipeline; do
  want=$([ $p = STAAR ] && echo $STAAR_SHA || echo $PIPE_SHA)
  have=$(grep -E "^RemoteSha:" "$DEFAULT_LIB/$p/DESCRIPTION" | awk '{print $2}')
  [ "$have" = "$want" ] || { echo "default-library $p RemoteSha $have != $want"; exit 3; }
done

rm -rf "$PKG_DIR/build" && mkdir -p "$PKG_DIR/build" && cd "$PKG_DIR/build"
cp -r "$STAAR_SRC" STAAR
cp -r "$PIPE_SRC" STAARpipeline
patch -p1 -d STAAR         < "$PKG_DIR/patches/STAAR_saddlepoint_root_check.patch"
patch -p1 -d STAARpipeline < "$PKG_DIR/patches/STAARpipeline_saddlepoint_root_check.patch"
grep -q "^RootcheckPatch:" STAAR/DESCRIPTION && grep -q "^RootcheckPatch:" STAARpipeline/DESCRIPTION
# the exported C++ signatures are unchanged, so RcppExports.* stay valid; make sure of it
grep -q "double NR_Binary_SPA(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter)" STAAR/src/NR_Binary_SPA.cpp

mkdir -p "$ROOTCHECK_LIB"
export MAKEFLAGS="${MAKEFLAGS:--j2}"
echo "[$(date)] installing STAAR into $ROOTCHECK_LIB"
R CMD INSTALL --no-docs --no-multiarch --library="$ROOTCHECK_LIB" STAAR 2>&1 | tee STAAR_install.log | tail -3
echo "[$(date)] installing STAARpipeline into $ROOTCHECK_LIB"
R CMD INSTALL --no-docs --no-multiarch --library="$ROOTCHECK_LIB" STAARpipeline 2>&1 | tee STAARpipeline_install.log | tail -3

# verify: the patched packages load from the new library, the originals are untouched
Rscript -e '
lib <- Sys.getenv("ROOTCHECK_LIB"); def <- Sys.getenv("DEFAULT_LIB")
for (p in c("STAAR","STAARpipeline")) {
  d <- packageDescription(p, lib.loc = lib); stopifnot(!is.null(d$RootcheckPatch))
  o <- packageDescription(p, lib.loc = def); stopifnot(is.null(o$RootcheckPatch))
  cat(sprintf("%-14s rootcheck lib: version %s | %s\n", p, d$Version, d$RootcheckPatch))
  cat(sprintf("%-14s default lib  : version %s RemoteSha %s (unchanged)\n", p, o$Version, o$RemoteSha))
}
suppressPackageStartupMessages({ library(STAAR); library(STAARpipeline) })
stopifnot(normalizePath(dirname(system.file(package = "STAAR"))) == normalizePath(lib))
stopifnot(normalizePath(dirname(system.file(package = "STAARpipeline"))) == normalizePath(lib))
cat("loaded from:", system.file(package = "STAAR"), "\n")
'
{
  echo "BUILD_RECORD  $(date -Iseconds)  host $(hostname)  user $USER"
  echo "R: $(R --version | head -1)"
  echo "compiler: $(R CMD config CXX) / $(gcc --version | head -1)"
  echo "library: $ROOTCHECK_LIB"
  echo "STAAR          upstream commit $STAAR_SHA   tarball sha256 $(sha256sum "$PKG_DIR/src_untrusted/STAAR-$STAAR_SHA.tar.gz" | cut -d' ' -f1)"
  echo "STAARpipeline  upstream commit $PIPE_SHA   tarball sha256 $(sha256sum "$PKG_DIR/src_untrusted/STAARpipeline-$PIPE_SHA.tar.gz" | cut -d' ' -f1)"
  echo "patch sha256: $(sha256sum "$PKG_DIR/patches/STAAR_saddlepoint_root_check.patch" | cut -d' ' -f1)  STAAR_saddlepoint_root_check.patch"
  echo "patch sha256: $(sha256sum "$PKG_DIR/patches/STAARpipeline_saddlepoint_root_check.patch" | cut -d' ' -f1)  STAARpipeline_saddlepoint_root_check.patch"
  for p in STAAR STAARpipeline; do
    echo "installed $p: $(grep -E '^(Version|RootcheckPatch):' "$ROOTCHECK_LIB/$p/DESCRIPTION" | tr '\n' ' ')"
    echo "installed $p .so sha256: $(sha256sum "$ROOTCHECK_LIB/$p/libs/$p.so" | cut -d' ' -f1)"
  done
  echo "default library untouched: STAAR $(grep -E '^RemoteSha' "$DEFAULT_LIB/STAAR/DESCRIPTION") ; STAARpipeline $(grep -E '^RemoteSha' "$DEFAULT_LIB/STAARpipeline/DESCRIPTION")"
} > "$PKG_DIR/build/BUILD_RECORD.txt"
cat "$PKG_DIR/build/BUILD_RECORD.txt"
echo "BUILD_OK $(date -Iseconds)" > "$PKG_DIR/build/BUILD_OK"
echo "[$(date)] BUILD_OK"
