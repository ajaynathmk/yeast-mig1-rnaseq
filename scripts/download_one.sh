#!/usr/bin/env bash
# Usage: scripts/download_one.sh <url>
# Download ONE FASTQ into data/raw. The scheme in the url is ignored: https, http and ftp are
# tried in that order (ENA serves some files on only some of them). Partial files are resumed.
# Success means the file size equals the size ENA lists in metadata/fastq_sizes.tsv.
set -u
dest=data/raw
mkdir -p "$dest"
url=$1
path=${url#*://}
f=$(basename "$path")
want=$(awk -F'\t' -v f="$f" '$1==f {print $2}' metadata/fastq_sizes.tsv)
have() { if [ -f "$dest/$f" ]; then stat -c %s "$dest/$f"; else echo 0; fi; }

[ -z "$want" ] && { echo "[FAIL] $f: no expected size in metadata/fastq_sizes.tsv"; exit 1; }
[ "$(have)" = "$want" ] && { echo "[ok]   $f already complete"; exit 0; }

for scheme in https http ftp; do
  echo "[try]  $f via $scheme"
  wget -c -q -t 3 --waitretry=20 -P "$dest" "$scheme://$path"
  [ "$(have)" = "$want" ] && { echo "[done] $f"; exit 0; }
  echo "[warn] $f incomplete after $scheme ($(have) of $want bytes)"
done
echo "[FAIL] $f"
exit 1
