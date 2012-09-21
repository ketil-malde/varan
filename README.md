# Brief user guide for deltasnp

This is a small guide for an even smaller program that analyses gene
variant data from pooled samples, with the intention of identifying
the variants that have different underlying allele frequencies in the
pools.

## Usage

First you need to align your pooled reads into one BAM-file each, say
`sample1.bam`, `sample2.bam`, etc.
Basically, you should be able to get away with something like:

~~~~~~
samtools mpileup -f ref.fasta sample1.bam sample2.bam sample3.bam | deltasnp | grep '[+*]' > OUT
~~~~~~

The output from `samtools pileup` looks something like this:

~~~~~~
                                <---sample #1 -->       <--sample #2 --->       <--sample#3----->  
   contig       pos    ref    #reads   char    qual   #reads .... etc
ATLCOD1Ac00001  30      A       1       ^>.     @       0       *       *       0       *       *
ATLCOD1Ac00001  31      C       1       .       @       0       *       *       0       *       *
ATLCOD1Ac00001  32      A       1       .       @       0       *       *       0       *       *
ATLCOD1Ac00001  33      C       1       .       D       0       *       *       0       *       *
ATLCOD1Ac00001  34      A       1       .       F       0       *       *       0       *       *
ATLCOD1Ac00001  35      C       1       .       D       0       *       *       0       *       *
ATLCOD1Ac00001  36      A       1       .       D       0       *       *       0       *       *
ATLCOD1Ac00001  37      C       1       .       D       0       *       *       0       *       *
ATLCOD1Ac00001  38      A       1       .       B       0       *       *       0       *       *
ATLCOD1Ac00001  39      C       1       .       A       0       *       *       0       *       *
~~~~~~

## Output

After filtering this through `deltansp` (and `grep` as above,
otherwise you'd see quite a larger number of lines), we should get
something like:

~~~~~~
                                <-lib 1-->      <-lib2-->       <-lib3-->               conf2   conf3
                                 A  C G T        A C G T         A C G T                ACGT    ACGT             
ATLCOD1Ac00006  227     C       (0,11,0,0)      (0,0,0,0)       (0,2,0,4)       -       ....    .+.+
ATLCOD1Ac00014  309     C       (3,16,0,0)      (0,7,0,0)       (0,72,0,0)      -       ....    ++..
ATLCOD1Ac00014  466     T       (0,0,0,5)       (0,0,0,1)       (10,2,4,8)      -       ....    ...+
ATLCOD1Ac00014  471     G       (0,0,4,0)       (0,0,0,0)       (8,0,2,0)       -       ....    +.+.
ATLCOD1Ac00014  482     G       (0,0,4,0)       (0,0,0,0)       (12,0,12,12)    -       ....    ..+.
ATLCOD1Ac00021  350     C       (0,7,0,0)       (0,1,0,0)       (0,3,0,10)      -       ....    .+.+
ATLCOD1Ac00023  87      C       (0,33,0,0)      (0,1,1,0)       (0,14,0,0)      -       .++.    ....
ATLCOD1Ac00027  101     G       (8,0,1,0)       (0,0,0,0)       (0,0,3,0)       -       ....    +.+.
ATLCOD1Ac00027  314     T       (1,0,0,9)       (0,0,0,0)       (6,0,0,2)       -       ....    +..+
ATLCOD1Ac00027  408     C       (0,13,0,1)      (0,0,0,0)       (0,5,0,6)       -       ....    .+.+
~~~~~~

Here we see that each library has had the allele frequencies counted
up, and the tuples number As, Cs, Gs, and Ts for each position.  The
final columns compare the samples to sample number 1, and a `+`
indicates that the 95% confidence intervals for that allele's
frequency don't overlap, a `*` is the same, but for 99% confidence.

So in position 227 in contig 6 we see a significant difference
in allele frequency for `C` and `T`.  This is a bit surprising, since
`C` is the major allele (in fact the only one) in the first sample.
But the confidence intervals are

~~~~~~
> confidenceInterval 1.65 11 0           -- 11 "successful" observasjoner, no "failures"
(0.7676537322704904,1.0339494741423352)  -- 95% interval above 0.76
> confidenceInterval 1.65 2 4            -- 2 "successes", 4 "failures"
(0.11345603146594596,0.657251907771658)  -- 95% interval less than 0.66
~~~~~~

Specifically I use Agresti-Coull's confidence interval approximation,
which according to the same is a better way to do it than to calculate
the binomials directly.
