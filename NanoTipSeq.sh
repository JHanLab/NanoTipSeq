#!/bin/bash

# 09/10/23 v4
# author: Xuanming Zhang, Sam Stuckert
# usage: bash NanoTipSeq.sh <alignment.bam> <polyA count>
# input a nanopore alignment file and output potential insertions region based on soft-clipped reads
# L1 index
set -e
set -o pipefail
L1MdTf_btw2_index=#insert path to bowtie index

# Packages and environments
# Make sure these 2 packages are installed
# bowtie2
# ripgrep https://github.com/BurntSushi/ripgrep/blob/master/GUIDE.md

samextractclip= #insert path to jar file
echo "Index Path: $L1MdTf_btw2_index"
echo "ClipJarPath: $samextractclip"


# Input
input_bam=$1
echo "Bam file: $input_bam"

# Checks for user input of polyA count, if nothing defaults to 19 As, otherwise uses user input to determine A count
echo "Checking for polyA"
if [ -z "$2" ]; then
  polya="AAAAAAAAAAAAAAAAAA"
  ACount="19"
else
  polya=$(printf '%0.sA' $(seq 1 "$2"))
  ACount=$2
fi

# Creates output directory for all pipeline files
echo "Creating pipeline directory"
output_dir="${1%.bam}"
mkdir -p "$output_dir"
echo "Pipeline directory created: $output_dir"
echo "Using a polyA count of $ACount"

# Check if pre-required file is existing
echo "Checking for prereq files"
if ! ls "$output_dir"/*L1mdtf1.sorted.bam &>/dev/null; then
  echo "Prereq files not found. Extracting soft clip reads."
  for i in $input_bam; do
    echo "$output_dir/${i%.bam}soft_clip.fastq"
    java -jar "$samextractclip" "$i" > "$output_dir/${i%.bam}soft_clip.fastq"
  done

  # Trim query name
  echo "Trimming Query Name"
  for i in "$output_dir"/*soft_clip.fastq; do
    cat "$i" | paste - - - - | cut -f 1 | cut -f 1-2 -d ";" | paste - <(cat "$i" | paste - - - - | cut -f 2-4) | sed 's/\t/\n/g' > "${i%.fastq}.short_name.fastq"
  done

  # Map to L1MdTf1.fa bowtie2
  echo "Mapping reference genome using bowtie2"
  for i in "$output_dir"/*short_name.fastq; do
    date >> "$output_dir/stats.file.txt"
    echo "$i" >> "$output_dir/stats.file.txt"
    bowtie2 --sensitive -p 20 -N 1 --mp 1,1 --rdg 5,2 --rfg 5,2 -x "$L1MdTf_btw2_index" -U "$i" 2>> "$output_dir/stats.file.txt" | samtools view -hb - | samtools sort -o "${i%.fastq}L1mdtf1.sorted.bam" -
  done
fi

# To filter out useful mapped reads:
echo "Filtering reads"
for i in "$output_dir"/*L1mdtf1.sorted.bam; do
  echo $i
  samtools view -F 4 "$i" | grep  "$polya" | cut -f 1 | sed "s/;/_/g" > "${i%L1mdtf1.sorted.bam}L1.mapped.reads_contain-$ACount""pA"".name.txt"
done

# Extract genomic location of the reads mapped to L1 that contain poly A tails
echo "Extracting L1 reads"
for i in "$output_dir"/*soft_clip.fastq; do
  cat "$i" | awk 'NR%4==1 {print substr($1,2)}' | sed "s/;/\t/g" | sed "s/'/ /g" | awk '{OFS="\t";$1=$1"_"$2; print}' > "${i%.fastq}.name.txt"
  time rg -f "${i%fastq}short_nameL1.mapped.reads_contain-$ACount""pA"".name.txt" "${i%.fastq}.name.txt" > "${i%fastq}temp.test_name.txt"
  cat "${i%fastq}temp.test_name.txt" | awk -v OFS="\t" '$7==5 {print $3, $4-150, $4+20, $4 , "plus" , $1}' > "${i%fastq}temp.plus.txt"
  cat "${i%fastq}temp.test_name.txt" | awk -v OFS="\t" '$7==3 { chr = $3; start = $4; cigar = $6; len = 0;
  while (match(cigar, /[0-9]+[MINX=]/)) {
      len += substr(cigar, RSTART, RLENGTH - 1);
      cigar = substr(cigar, RSTART + RLENGTH);
  }
  $8 = start + len;
  $9 = $8 + 150;
  print $3, $8, $9, $4, "minus", $1;}' > "${i%fastq}temp.minus.txt"
  cat "${i%fastq}temp.plus.txt" "${i%fastq}temp.minus.txt" | sort -k1,1 -k2,2n - > "${i%fastq}$ACount""pA"".merged.bed"
  rm "${i%fastq}temp"*
done

echo "Checking for merged.bed"
if ls "$output_dir"/*"$ACount""pA"".merged.bed" &>/dev/null; then
  echo "merged.bed was created"
else
  echo "ERROR: merged.bed was not created"
fi

# Step 3: Cluster
# Make sure to sort the bed file
echo "Clustering region"

for i in "$output_dir"/*"$ACount""pA"".merged.bed"; do
  bedtools cluster -i "$i" > "${i%bed}cluster.bed"
  cat "${i%bed}cluster.bed" | cut -f 7 | sort | uniq -c | awk '$1>3 {print $2}' > "${i%bed}tempID"
  while read -r line; do
    awk -v r="$line" '$7==r' "${i%bed}cluster.bed" >> "${i%bed}cluster.3hits.2.bed"
  done <"${i%bed}tempID"
done

# Checking outputs
echo "Checking outputs"
if ls "$output_dir"/*"$ACount""pA"".merged.cluster.bed" &>/dev/null; then
  echo "cluster.bed was created"
else
  echo "ERROR: cluster.bed was not created"
fi

if ls "$output_dir"/*"$ACount""pA"".merged.tempID" &>/dev/null; then
  echo "tempID was created"
else
  echo "ERROR: tempID was not created"
fi

if ls "$output_dir"/*"$ACount""pA"".merged.cluster.3hits.2.bed" &>/dev/null; then
  echo "cluster.3hits.2.bed was created"
else
  echo "ERROR: cluster.3hits.2.bed was not created"
fi

for i in "$output_dir"/*"$ACount""pA"".merged.cluster.3hits.2.bed"; do
  bedtools merge -i <(sort -k1,1 -k2,2n "$i") > "${i%cluster.3hits.2.bed}new_inst-3_sup.txt"
done

if ls "$output_dir"/*"$ACount""pA"".merged.new_inst-3_sup.txt" &>/dev/null; then
  echo "Predicted insertions list created"
else
  echo "ERROR: Predicted insertions list was not created"
fi
