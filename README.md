# Brief user guide for varan

This is a small guide for an even smaller program that analyses gene
variant data from pooled samples, with the intention of identifying
the variants that have different underlying allele frequencies in the
pools.

## Usage

First you need to align your pooled reads into one BAM-file each, say
`sample1.bam`, `sample2.bam`, etc.
Basically, you should be able to get away with something like:

~~~~~~
samtools mpileup -f ref.fasta sample1.bam sample2.bam sample3.bam | varan [options] > OUT
~~~~~~

The output from `samtools pileup` looks something like this:

~~~~~~
                                        [-- sample #1 -->               [-- sample #2 --->    
   contig       pos    ref    #reads    char            qual            #reads .... etc
LSalAtl2s1      6       A       12      ,,,,,,,,,,.,    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      7       A       12      ,,,,,,,,,,.,    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      8       A       12      ,,,,,,,,,,.,    ::::::::::::    9       ,,,,,,,,^F,     :::::::::
LSalAtl2s1      9       A       12      ,,,,,,,,,,.,    ::::::::::::    9       ,,,,,,,,,       :::::::::
LSalAtl2s1      10      T       13      ,,,,,,,,,,.,^F. :::::::::::::   9       ,,,,,,,,,       :::::::::
LSalAtl2s1      11      C       13      ,,,,,,,,,,.,.   :::::::::::::   9       ,,,,,,,,,       :::::::::
LSalAtl2s1      12      A       13      ,,,,,,,,,,.,.   :::::::::::::   9       ,,,,,,,,,       :::::::::
LSalAtl2s1      13      T       13      ,,,,,,,,,,.,.   :::::::::::::   9       ,$,,,,,,,,      :::::::::
LSalAtl2s1      14      C       13      ,$,,,,,,,,,.,.  :::::::::::::   8       ,,,,,,,,        ::::::::
LSalAtl2s1      15      A       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      16      T       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      17      A       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      18      A       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      19      T       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
LSalAtl2s1      20      A       12      ,,,,,,,,,.,.    ::::::::::::    8       ,,,,,,,,        ::::::::
~~~~~~

## Output

After filtering this through `varan` as above, you get

~~~~~~
                                 <-lib 1->       <-lib 2->              conf    angle   Fst  PIk  dCI  Indels
                                 A  C G T        A  C G T               ACGT    
LSalAtl2s1      153     A        18 0 0 0        20 0 0 0       -       ....    1.000   0.000   0.000   0.04    []              
LSalAtl2s1      154     T        0 0 0 17        0 0 0 20       -       ....    1.000   0.000   0.000   0.06    []              
LSalAtl2s1      155     A        17 0 0 0        20 0 0 0       -       ....    1.000   0.000   0.000   0.06    []              
LSalAtl2s1      156     T        0 0 0 17        0 0 0 20       -       ....    1.000   0.000   0.000   0.06    []              
LSalAtl2s1      157     A        16 0 0 0        20 0 0 0       -       ....    1.000   0.000   0.000   0.08    []              
LSalAtl2s1      158     C        1 15 0 0        10 10 0 0      -       ++..    0.753   0.223   0.500   2.31    []              
LSalAtl2s1      159     A        14 2 0 0        19 1 0 0       -       ....    0.996   0.018   0.163   0.53    []              
LSalAtl2s1      160     A        17 0 0 0        20 0 0 0       -       ....    1.000   0.000   0.000   0.06    []              
LSalAtl2s1      161     T        0 0 0 18        0 0 0 20       -       ....    1.000   0.000   0.000   0.04    []              
LSalAtl2s1      162     T        0 0 0 18        0 0 0 19       -       ....    1.000   0.000   0.000   0.02    []              
LSalAtl2s1      163     G        0 0 18 0        0 0 19 0       -       ....    1.000   0.000   0.000   0.02    []              

~~~~~~

*NOTE*: This is a bit out of date, and describes deprecated or
unfinished features.  Use as an example of output format only, not as
a recommendation of statistics to use!

Here we see that each input library (i.e. BAM file) has had the allele frequencies counted
up as the number of As, Cs, Gs, and Ts for each position.  Note that a
space is used to separate these, while TAB is used between the columns
proper.  

The `conf` column for position 158 displays a `+` 
indicating that the 95% confidence intervals for that allele's
frequency don't overlap, a `*` would indicate 99% confidence.  It is
perhaps useful to `grep` on this column to select a subset of the
output.

Then follows the angle between allele frequency spectra, the $F_ST$,
$\pi_k$, and $\Delta_{CI}$, the latter is the distance between
the confidence intervals measured in standard deviations. 
Specifically I use Agresti-Coull's confidence interval approximation,
which according to Agresti and Coull is a better way to do it than to
calculate the binomials directly.

The final column lists indels and their counts, if any.

# Examining specific regions

It is also easy to extract a specific region (e.g. a specific gene)
through `samtools`, if you index the fasta-file first.  For instance:

    samtools mpileup -f ref.fasta -r CHR8:4500-6000 sample1.bam sample2.bam | varan > out

