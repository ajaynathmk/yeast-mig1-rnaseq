#!/usr/bin/env bash
# Download every FASTQ listed in metadata/fastq_urls.txt into data/raw.
# ENA serves some files only over https/http and others only over ftp, so each file
# is tried over https, then http, then ftp. Partial files are resumed (wget -c).
# A file counts as done only when its size equals the size ENA lists (metadata/fastq_sizes.tsv).
set -u
dest=data/raw
mkdir -p "$dest"

want_size() { awk -F'\t' -v f="$1" '$1==f {print $2}' metadata/fastq_sizes.tsv; }
have_size() { if [ -f "$dest/$1" ]; then stat -c %s "$dest/$1"; else echo 0; fi; }

while read -r url; do
  path=${url#*://}                 # drop the scheme, keep host/path
  f=$(basename "$path")
  want=$(want_size "$f")
  if [ -z "$want" ]; then echo "[FAIL] $f: no expected size in metadata/fastq_sizes.tsv"; continue; fi
  if [ "$(have_size "$f")" = "$want" ]; then echo "[ok]   $f already complete"; continue; fi

  for scheme in https http ftp; do
    echo "[try]  $f via $scheme"
    wget -c -q -t 3 --waitretry=20 -P "$dest" "$scheme://$path"
    if [ "$(have_size "$f")" = "$want" ]; then echo "[done] $f"; break; fi
    echo "[warn] $f incomplete after $scheme ($(have_size "$f") of $want bytes)"
  done

  [ "$(have_size "$f")" = "$want" ] || echo "[FAIL] $f"
  sleep 5
done < metadata/fastq_urls.txt
echo "[end]  pass finished"
