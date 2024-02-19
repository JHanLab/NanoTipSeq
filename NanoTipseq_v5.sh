#!/bin/bash

# 09/10/23 v4 
# author: Xuanming Zhang, Sam Stucker 
# usage: bash this_script.sh bamfile number_of_polyA
# input a nanopore alignment file and output potential insertions region base on soft-clipped reads

# L1 index 
L1MdTf_btw2_index=/Volumes/SanDisk_SSD/L1MdTf1_bt2_index/L1MdTf1


#packages and enviroments 

#make sure these 2 pacakages are installed 
#bowtie2 
#ripgrep https://github.com/BurntSushi/ripgrep/blob/master/GUIDE.md

filterbyname=/Users/xuanmingzhang/bbmap/filterbyname.sh
samextractclip=/Users/xuanmingzhang/jvarkit/dist/samextractclip.jar
bbmap_reformat=/Users/xuanmingzhang/bbmap/reformat.sh
eval "$(conda shell.bash hook)"
conda activate bioinfo


# input 
input_bam=$1

#Checks for user input of polyA count, if nothing defaults to 19 As, otherwise uses user input to determine A count
if [ -z "$2" ]
  then
    polya="AAAAAAAAAAAAAAAAAA"
    ACount="19"
else
    polya=$(printf '%0.sA' $(seq 1 $2))
    ACount=$2
fi

#Creates output directory for all pipeline files
mkdir ${1%".bam"}
output_dir=${1%".bam"}

echo "Using a polyA count of $ACount"

# check if pre-required file is existing 
if [ -f $output_dir/*L1mdtf1.sorted.bam ]
then
  echo "sorted.bam file already exists, skipping creation"
else
# extract soft clip reads : 
# will use samextractclip.jar from javarkit/dist folder , link : http://lindenb.github.io/jvarkit/SamExtractClip.html
  for i in $input_bam ; do
    java -jar $samextractclip $i > $output_dir/${i%bam}"soft_clip.temp.fastq"
    $bbmap_reformat in=$output_dir/${i%bam}"soft_clip.temp.fastq" out=$output_dir/${i%bam}"soft_clip.fastq" qin=33 minlength=80 maxlength=200
    rm $output_dir/${i%bam}"soft_clip.temp.fastq"
  done

# trim query name 

  for i in $output_dir/*soft_clip.fastq ; do 
    cat $i | paste - - - - | cut -f 1  | cut -f 1-2 -d ";" | paste - <( cat $i | paste - - - - | cut -f 2-4) | sed 's/\t/\n/g' > ${i%fastq}short_name.fastq
  done

# map to L1MdTf1.fa bowtie2 
  for i in $output_dir/*short_name.fastq ; do 
    date >> $output_dir/stats.file.txt 
    echo $i >> $output_dir/stats.file.txt
    bowtie2 --sensitive -p 20 -N 1 --mp 1,1 --rdg 5,2 --rfg 5,2 -x $L1MdTf_btw2_index -U $i 2>> $output_dir/stats.file.txt |  samtools view -hb - |  samtools sort -o ${i%fastq}"L1mdtf1.sorted.bam" -
  done
fi

# To filter out useful mapped reads: 
for i in $output_dir/*L1mdtf1.sorted.bam ; do
  samtools view -F 4 $i | grep  $polya |  cut -f 1 | sed "s/;/_/g" > ${i%L1mdtf1.sorted.bam}L1.mapped.reads_contain-$ACount"pA".name.txt
done 
# will test bbduk to search for primer sequence instead of poly-As 
# bbduk.sh in=reads.fq out=unmatched.fq outm=matched.fq literal=ACGTACGTACGTACGTAC k=18 mm=f hdist=2 
# Make sure "k" is set to the exact length of the sequence. "hdist" controls the number of substitutions allowed. "outm" gets the reads that match. By default this also looks for the reverse-complement; you can disable that with "rcomp=f". 
# or 
# I can extract by mapped location : 

# for i in `find . -name *L1mdtf1.sorted.bam | grep 10pA`; do 
#     samtools view -F 4 $i  |  awk '$4>6800 && length($10)>20' | cut -f 4 > $i.map.position.txt 
#     samtools view -F 4 $i  |  awk '$4>6800 && length($10)>20' | awk '{print length($10)}' > $i.map.length.txt 
# done

# extract genomic location of the reads map to L1 that contain poly A tails 
  ###step 1 do not use Grep or awk, too slow for reliterately grepping. use bbmap / filterbyname.sh . 100x faster 



#maybe we can simply this by use grep 
#1 : print fastq name (since this contain all infor we need)
for i in $output_dir/*soft_clip.fastq ; do 
  cat $i | awk 'NR%4==1 {print substr($1,2)}' | sed "s/;/\t/g" | sed "s/'/ /g" | awk '{OFS="\t";$1=$1"_"$2; print}' > ${i%.fastq}.name.txt
    time rg -f ${i%fastq}"short_name.L1.mapped.reads_contain-"$ACount"pA.name.txt" ${i%.fastq}.name.txt > ${i%fastq}temp.test_name.txt 
    cat ${i%fastq}temp.test_name.txt  | awk -v OFS="\t" '$7==5 {print $3, $4-150, $4+20, $4 , "plus" , $1}' > ${i%fastq}temp.plus.txt
    cat ${i%fastq}temp.test_name.txt | awk -v OFS="\t" '$7==3 { chr = $3; start = $4; cigar = $6; len = 0;
    while (match(cigar, /[0-9]+[MINX=]/)) {
        len += substr(cigar, RSTART, RLENGTH - 1);
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    $8 = start + len;
    $9 = $8 + 150;
    print $3, $8, $9, $4, "minus", $1;}' > ${i%fastq}temp.minus.txt
    cat ${i%fastq}temp.plus.txt ${i%fastq}temp.minus.txt  | sort -k1,1 -k2,2n - > ${i%fastq}$ACount"pA.merged.bed"
    rm ${i%fastq}temp*
done

for i in $output_dir/*soft_clip.fastq ; do 
  cat $i | awk 'NR%4==1 {print substr($1,2)}' | sed "s/;/\t/g" | sed "s/'/ /g" | awk '{OFS="\t";$1=$1"_"$2; print}' > ${i%.fastq}.name.txt
    time rg -f ${i%fastq}"short_name.L1.mapped.reads_contain-"$ACount"pA.name.txt" ${i%.fastq}.name.txt > ${i%fastq}temp.test_name.txt 
    cat ${i%fastq}temp.test_name.txt  | awk -v OFS="\t" '$7==5 {print $3, $4-150, $4+20, $4 , "plus" , $1}' > ${i%fastq}temp.plus.txt
    cat ${i%fastq}temp.test_name.txt | awk -v OFS="\t" '$7==3 { chr = $3; start = $4; cigar = $6; len = 0;
    while (match(cigar, /[0-9]+[MINX=]/)) {
        len += substr(cigar, RSTART, RLENGTH - 1);
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    $8 = start + len;
    $9 = $8 + 150;
    print $3, $8, $9, $4, "minus", $1;}' > ${i%fastq}temp.minus.txt
    cat ${i%fastq}temp.plus.txt ${i%fastq}temp.minus.txt  | sort -k1,1 -k2,2n - > ${i%fastq}$ACount"pA.merged.bed"
    rm ${i%fastq}temp*
done

 
if [ -f $output_dir/*$ACount"pA.merged.bed" ]
then
  echo "merged.bed was created"
else
  echo "ERROR: merged.bed was not created"
fi



  ###step 3 cluster the region 
#makesure sort the bed file
for i in $output_dir/*$ACount"pA.merged.bed" ;do 
  bedtools cluster -i $i > ${i%bed}"cluster.bed" 
  cat ${i%bed}"cluster.bed" | cut -f 7 | sort | uniq -c | awk '$1>3 {print $2}' > ${i%bed}"tempID" 
    while read -r line ;do 
      awk -v r=$line '$7==r' ${i%bed}"cluster.bed"  >>   ${i%bed}"cluster.3hits.2.bed"
    done <${i%bed}"tempID"
done 


#Each following if statement checks if the correct files were created
#1. *merged.bed
#2. *cluster.bed
#3. *tempID
#4. *cluster.3hits.2.bed
#5. *samtools.merged.bed

if [ -f $output_dir/*$ACount"pA.merged.cluster.bed" ]
then
  echo "cluster.bed was created"
else
  echo "ERROR: cluster.bed was not created"
fi

if [ -f $output_dir/*$ACount"pA.merged.tempID" ]
then
  echo "tempID was created"
else
  echo "ERROR: tempID was not created"
fi

if [ -f $output_dir/*$ACount"pA.merged.cluster.3hits.2.bed" ]
then
  echo "cluster.3hits.2.bed was created"
else
  echo "ERROR: cluster.3hits.2.bed was not created"
fi


for i in $output_dir/*$ACount"pA.merged.cluster.3hits.2.bed" ; do 
  bedtools merge -i <(sort -k1,1 -k2,2n $i) > ${i%cluster.3hits.2.bed}"new_inst-3_sup.txt"
done

if [ -f $output_dir/*$ACount"pA.merged.new_inst-3_sup.txt" ]
then
  echo "predicited insertions list created"
else
  echo "ERROR: predicited insertions list was not created"
fi
rm -f $output_dir/*soft_clip.fastq
rm -f $output_dir/*short_name.fastq
rm -f $output_dir/*L1mdtf1.sorted.bam
rm -f $output_dir/*.name.txt
rm -f $output_dir/*$ACount"pA.merged.bed"
rm -f $output_dir/*$ACount"pA.merged.cluster.bed"
rm -f $output_dir/*$ACount"pA.merged.tempID"

# write a code to check if in each input, there are the following output files : 

#1. *merged.bed
#2. *cluster.bed
#3. *tempID
#4. *cluster.3hits.2.bed
#5. *samtools.merged.bed

#  for i in `ls */*merged.sorted.soft_clip.merged.insertions-sup-3reads.txt` ; do wc -l $i ;done
#      259 barcode01.all_data.merged.sorted/barcode01.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      157 barcode02.all_data.merged.sorted/barcode02.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      113 barcode03.all_data.merged.sorted/barcode03.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      117 barcode04.all_data.merged.sorted/barcode04.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      203 barcode05.all_data.merged.sorted/barcode05.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      104 barcode06.all_data.merged.sorted/barcode06.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#       92 barcode07.all_data.merged.sorted/barcode07.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      215 barcode08.all_data.merged.sorted/barcode08.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#      176 barcode09.all_data.merged.sorted/barcode09.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt
#       87 barcode10.all_data.merged.sorted/barcode10.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt

# intervene upset --names=1,2,3,4,5,6,7,8,9,10 --save-overlaps --project highMM_oyctes  -i \
# barcode01.all_data.merged.sorted/barcode01.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode02.all_data.merged.sorted/barcode02.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode03.all_data.merged.sorted/barcode03.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode04.all_data.merged.sorted/barcode04.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode05.all_data.merged.sorted/barcode05.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode06.all_data.merged.sorted/barcode06.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode07.all_data.merged.sorted/barcode07.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode08.all_data.merged.sorted/barcode08.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode09.all_data.merged.sorted/barcode09.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt \
# barcode10.all_data.merged.sorted/barcode10.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.txt 

# for i in `ls */*merged.insertions-sup-3reads.txt`;
# do bedtools intersect -v -b /Volumes/SanDisk_SSD/L1MdTf_1-3_add500to3End.bed -a $i > ${i%txt}exclude_L1Tf.txt 
# done 

# intervene upset --names=2,3,4,6,7,9 -o highMM_excludeL1Tf --save-overlaps -i barcode02.all_data.merged.sorted/barcode02.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt barcode03.all_data.merged.sorted/barcode03.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt barcode04.all_data.merged.sorted/barcode04.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt barcode06.all_data.merged.sorted/barcode06.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt barcode07.all_data.merged.sorted/barcode07.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt barcode09.all_data.merged.sorted/barcode09.all_data.merged.sorted.soft_clip.merged.insertions-sup-3reads.exclude_L1Tf.txt


# chr12:10,981,205-10,981,341 could be de novo insertion wiht 3' transduction 

