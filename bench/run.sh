#!/bin/sh
# Run the multivector benchmarks on every implementation and save the results.
#
#   bench/run.sh                      # from the repository root
#   CURRY=/path/to/curry bench/run.sh # use a specific Curry binary
#
# Writes bench/results/<date>-<host>-<commit>.txt and prints it.  The commit
# is the code measured; "-dirty" means there were uncommitted changes.

set -eu

cd "$(dirname "$0")/.."

CHEZ=${CHEZ:-scheme}
CURRY=${CURRY:-$HOME/src/curry/build-release/curry}

commit=$(git rev-parse --short HEAD)
if ! git diff --quiet HEAD -- curry-r6rs bench; then
    commit="$commit-dirty"
fi

mkdir -p bench/results
out="bench/results/$(date +%Y-%m-%d)-$(hostname -s)-$commit.txt"

{
    echo "# $(date -Iseconds)"
    echo "# commit: $commit"
    echo "# cpu:   $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //' || uname -m)"
    echo "# os:    $(uname -sr)"
    echo "# chez:  $($CHEZ --version 2>&1)"
    if [ -x "$CURRY" ]; then
        echo "# curry: $($CURRY -v 2>&1 | head -1) ($CURRY)"
    else
        echo "# curry: not found at $CURRY"
    fi
    echo "#"
    echo "# columns: implementation workload algebra operations best-ms median-ms ns-per-op checksum"
    echo "#   chez     = (curry-r6rs multivector), Chez default (safe) optimize-level"
    echo "#   chez-o3  = the same at optimize-level 3 (unsafe)"
    echo "#   curry-c  = Curry's built-in C multivectors"

    "$CHEZ" --program bench/multivector.sps
    "$CHEZ" --optimize-level 3 --program bench/multivector.sps | sed 's/^chez /chez-o3 /'
    if [ -x "$CURRY" ]; then
        "$CURRY" --clear-cache bench/multivector-curry-c.scm
    fi
} | tee "$out"
