#!/bin/bash
# Example script to extract a region around SNPs for 
# genotyping primer design.  This will select one SNP from
# each of the $N longest contigs.  We also need:
#
# longest-contigs.txt - file with contig names sorted by length
# $BAMS               - the bam files to use
# $N                  - the number of contigs to pick SNPs from
# 
# usage: ./gensnps.sh inputfile > outfile
#
#  where inputfile is a file containing the potential SNP sites
#  sorted according to importance. (i.e. by sorting the output from 
#  varan on the appropriate column).

N=30
BAMS=$(echo bams/{ATL,ANT,PAC}_sorted.bam.bam bams/Dwf.bam)

CONTIGS=$(head -$N longest-contigs.txt | cut -f1)

# replace wildcard in char position 101 with a regexp-like expression
fixsnp(){ sed \
  -e 's/\(^[A-Za-z]\{100\}\)R/\1[A\/G]/g' \
  -e 's/\(^[A-Za-z]\{100\}\)Y/\1[C\/T]/g' \
  -e 's/\(^[A-Za-z]\{100\}\)S/\1[C\/G]/g' \
  -e 's/\(^[A-Za-z]\{100\}\)W/\1[A\/T]/g' \
  -e 's/\(^[A-Za-z]\{100\}\)K/\1[G\/T]/g' \
  -e 's/\(^[A-Za-z]\{100\}\)M/\1[A\/C]/g' 
}

for a in $CONTIGS; do
	HIT=$(grep "^$a	" $1 | head -1)
	LOC=$(echo "$HIT" | cut -f2)
	echo \>$HIT
	START=$((LOC-100))
	END=$((LOC+100))
	samtools mpileup -B -r $a:$START-$END $BAMS 2> /dev/null | vextr | fixsnp
	samtools mpileup -B -r $a:$START-$END $BAMS 2> /dev/null | varan -s -v > "$a:$LOC.out"
	echo
done

