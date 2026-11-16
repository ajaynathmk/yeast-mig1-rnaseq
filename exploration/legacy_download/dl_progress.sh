#!/usr/bin/env bash
# Show download progress for every FASTQ in metadata/fastq_urls.txt.
# Expected sizes come from metadata/fastq_sizes.tsv (columns: file name, bytes).
printf '%-24s %9s %9s %6s\n' FILE DOWNLOADED TOTAL DONE
sum_have=0; sum_total=0
while read -r url; do
  f=$(basename "$url")
  total=$(awk -F'\t' -v f="$f" '$1==f {print $2}' metadata/fastq_sizes.tsv)
  [ -z "$total" ] && { printf '%-24s  (no expected size found)\n' "$f"; continue; }
  have=0; [ -f "data/raw/$f" ] && have=$(stat -c %s "data/raw/$f")
  printf '%-24s %9d %9d %5d%%\n' "$f" $((have/1048576)) $((total/1048576)) $((100*have/total))
  sum_have=$((sum_have+have)); sum_total=$((sum_total+total))
done < metadata/fastq_urls.txt
echo "--------------------------------------------------------"
printf '%-24s %9d %9d %5d%%\n' TOTAL $((sum_have/1048576)) $((sum_total/1048576)) $((100*sum_have/sum_total))
echo "wget processes running: $(pgrep -c wget)"
