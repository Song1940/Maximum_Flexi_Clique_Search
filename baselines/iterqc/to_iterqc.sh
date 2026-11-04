#!/bin/bash
# Convert an edge list of datasets/graphs to the text input of IterQC: a header "n m", then
# every undirected edge once as "u v" (u < v), without duplicates and self-loops, with node
# ids shifted to start at 0 (n = largest id + 1).
#
#   bash baselines/iterqc/to_iterqc.sh <edge list> <output> [<smallest node id>]
#
# The smallest node id of the input (0 or 1) is detected when it is not given. SORT_MEM
# (default 8G) and TMPDIR set the memory and the temporary directory of sort.
set -eo pipefail
src=$1; out=$2; base=$3
[ -n "$src" ] && [ -n "$out" ] || { echo "usage: $0 <edge list> <output> [<smallest node id>]" >&2; exit 1; }
if [ -z "$base" ]; then
  base=$(awk '$1 !~ /^[#%]/ && NF >= 2 { if (m == "" || $1 < m) m = $1; if ($2 < m) m = $2 } END { print m + 0 }' "$src")
fi
tmp=$(mktemp "${TMPDIR:-/tmp}/iterqc.XXXXXX")
trap 'rm -f "$tmp" "$out.partial"' EXIT
awk -v b="$base" '$1 !~ /^[#%]/ && NF >= 2 && $1 ~ /^[0-9]/ {
       u = $1 - b; v = $2 - b; if (u != v) { if (u < v) print u, v; else print v, u } }' "$src" \
  | LC_ALL=C sort -u -S "${SORT_MEM:-8G}" -T "${TMPDIR:-/tmp}" > "$tmp"
m=$(wc -l < "$tmp")
n=$(awk 'BEGIN { mx = 0 } { if ($1 > mx) mx = $1; if ($2 > mx) mx = $2 } END { print mx + 1 }' "$tmp")
{ echo "$n $m"; cat "$tmp"; } > "$out.partial"
mv "$out.partial" "$out"
echo "$src: n=$n m=$m -> $out"
