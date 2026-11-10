#!/usr/bin/env bash
# Download SK-2005 from the Laboratory for Web Algorithmics (LAW) and convert it to the input format
# of `flexi`: graphs/SK.txt (any other LAW graph can be given by its basename, e.g. it-2004,
# which is written to graphs/IT.txt).
#
#   bash graphs/convert_law.sh                  # SK-2005
#   bash graphs/convert_law.sh cnr-2000         # another LAW graph (basename of its files)
#
# The graphs are distributed in the compressed WebGraph format. They are symmetrised and their
# self-loops removed with WebGraph 3.6.12 (Transform simplifyOffline), written as an arc list
# (ArcListASCIIGraph), and every undirected edge is kept once (u < v), with the LAW node ids.
#
# Requirements: Java 11+, curl, md5sum, awk. The WebGraph jars are fetched from Maven Central.
# Disk: the final edge list of SK takes 32 GB; the intermediate arc list, which holds every edge in
# both directions, takes about twice as much. Memory: set JAVA_HEAP (default 64g;
# the paper used 200g).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-$HERE}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
LIB="$OUT/_webgraph"
WORK="$OUT/_download"
HEAP="${JAVA_HEAP:-64g}"
LAW_URL="http://data.law.di.unimi.it/webdata"
MAVEN="https://repo1.maven.org/maven2"

# WebGraph 3.6.12 and its runtime dependencies (group path / artifact / version).
JARS="
it/unimi/dsi webgraph 3.6.12
it/unimi/dsi dsiutils 2.7.4
it/unimi/dsi fastutil 8.5.16
it/unimi/dsi sux4j 5.4.1
it/unimi/di jsap 20210129
com/google/guava guava 33.4.0-jre
com/google/guava failureaccess 1.0.2
com/google/guava listenablefuture 9999.0-empty-to-avoid-conflict-with-guava
com/google/code/findbugs jsr305 3.0.2
com/google/errorprone error_prone_annotations 2.36.0
com/google/j2objc j2objc-annotations 3.0.0
org/checkerframework checker-qual 3.43.0
org/apache/commons commons-collections4 4.4
org/apache/commons commons-configuration2 2.11.0
org/apache/commons commons-lang3 3.12.0
org/apache/commons commons-math3 3.6.1
org/apache/commons commons-text 1.12.0
commons-logging commons-logging 1.3.2
org/slf4j slf4j-api 2.0.3
ch/qos/logback logback-classic 1.3.4
ch/qos/logback logback-core 1.3.4
net/sf/jung jung-algorithms 2.1.1
net/sf/jung jung-api 2.1.1
net/sf/jung jung-io 2.1.1
com/sun/mail javax.mail 1.6.2
javax/activation activation 1.1
"

mkdir -p "$LIB" "$WORK"
CP=""
while read -r group artifact version; do
    [ -z "$group" ] && continue
    jar="$LIB/$artifact-$version.jar"
    if [ ! -s "$jar" ]; then
        echo "fetching $artifact-$version.jar"
        curl -fsSL -o "$jar.part" "$MAVEN/$group/$artifact/$version/$artifact-$version.jar"
        mv "$jar.part" "$jar"
    fi
    CP="${CP:+$CP:}$jar"
done <<< "$JARS"
JAVA=(java "-Xmx$HEAP" -cp "$CP")

abbr() { case "$1" in it-2004) echo IT ;; sk-2005) echo SK ;; *) echo "$1" ;; esac; }

[ $# -eq 0 ] && set -- sk-2005
for n in "$@"; do
    dst="$OUT/$(abbr "$n").txt"
    if [ -s "$dst" ]; then echo "[$n] $dst exists"; continue; fi
    W="$WORK/$n"; mkdir -p "$W"; cd "$W"
    for ext in graph properties md5sums; do
        if [ ! -s "$n.$ext" ]; then
            echo "[$n] downloading $n.$ext"
            curl -fSL -o "$n.$ext.part" "$LAW_URL/$n/$n.$ext"
            mv "$n.$ext.part" "$n.$ext"
        fi
    done
    for ext in graph properties; do
        grep " $n\.$ext\$" "$n.md5sums" | md5sum -c - || { echo "[$n] MD5 check failed for $n.$ext"; exit 1; }
    done
    echo "[$n] offsets"
    "${JAVA[@]}" it.unimi.dsi.webgraph.BVGraph -o -O -L "$n"
    echo "[$n] symmetrise and remove self-loops"
    mkdir -p tmp
    "${JAVA[@]}" it.unimi.dsi.webgraph.Transform simplifyOffline "$n" "$n-sym" 1000000 tmp
    echo "[$n] arc list"
    "${JAVA[@]}" it.unimi.dsi.webgraph.ArcListASCIIGraph "$n-sym" "$n-sym.arcs"
    echo "[$n] keep each edge once"
    { printf '# %s from LAW (symmetrised, self-loops removed): %s/%s/\n' "$n" "$LAW_URL" "$n"
      awk -F'\t' '$1 < $2' "$n-sym.arcs"; } > "$dst.part"
    mv "$dst.part" "$dst"
    rm -rf tmp "$n-sym.arcs"
    echo "[$n] wrote $dst ($(grep -vc '^#' "$dst") edges)"
done
