#!/bin/sh
# Check that (curry-r6rs multivector) on Chez computes exactly what Curry's
# built-in C multivectors compute, operation by operation.
#
#   tests/differential.sh                      # from the repository root
#   CURRY=/path/to/curry tests/differential.sh
#
# Expect an exact match against a Curry built without FMA contraction (GCC's
# x86-64 default).  On arm64, or with -march=native, the C compiler fuses
# multiply-adds and the last bit of some products differs.

set -eu

cd "$(dirname "$0")/.."

CHEZ=${CHEZ:-scheme}
CURRY=${CURRY:-$HOME/src/curry/build-release/curry}

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$CHEZ" --program tests/differential.sps > "$tmp/chez.txt"
"$CHEZ" --optimize-level 3 --program tests/differential.sps > "$tmp/chez-o3.txt"
"$CURRY" --clear-cache tests/differential-curry-c.scm > "$tmp/curry.txt"

echo "Chez optimize-level 3 vs Chez default:"
"$CHEZ" --program tests/compare-results.sps "$tmp/chez.txt" "$tmp/chez-o3.txt"
echo "Curry's C code vs Chez:"
"$CHEZ" --program tests/compare-results.sps "$tmp/curry.txt" "$tmp/chez.txt"
