# Brief user guide for varan

This is a small guide for a set of not-so-small-anymore programs that
analyse gene variant data from pooled samples, with the intention of
identifying the variants that have different underlying allele
frequencies in the pools.

## Usage for `varan`

First you need to align your pooled reads into one BAM-file each, say
`sample1.bam`, `sample2.bam`, etc.
Basically, you should then be able to get away with something like:

~~~~~~
samtools mpileup -f ref.fasta sample1.bam sample2.bam sample3.bam | varan [options] > OUT
~~~~~~

The output from `samtools mpileup` looks something like this:

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

You can of course also feel this file directly to `varan`, e.g.

~~~~~~
varan [options] sample1.mpileup > OUT
~~~~~~

## Output

After filtering this through `varan` as above, you get something that
looks like:

~~~~~~
                                 <-lib 1->       <-lib 2->              conf    angle   Fst  PIk  dCI  Indels
                                 A  C G T        A  C G T               ACGT    
LSalAtl2s1      153     A        18:0:0:0        20:0:0:0       -       ....    1.000   0.000   0.000   0.04    []
LSalAtl2s1      154     T        0:0:0:17        0:0:0:20       -       ....    1.000   0.000   0.000   0.06    []
LSalAtl2s1      155     A        17:0:0:0        20:0:0:0       -       ....    1.000   0.000   0.000   0.06    []
LSalAtl2s1      156     T        0:0:0:17        0:0:0:20       -       ....    1.000   0.000   0.000   0.06    []
LSalAtl2s1      157     A        16:0:0:0        20:0:0:0       -       ....    1.000   0.000   0.000   0.08    []
LSalAtl2s1      158     C        1:15:0:0        10:10:0:0      -       ++..    0.753   0.223   0.500   2.31    []
LSalAtl2s1      159     A        14:2:0:0        19:1:0:0       -       ....    0.996   0.018   0.163   0.53    []
LSalAtl2s1      160     A        17:0:0:0        20:0:0:0       -       ....    1.000   0.000   0.000   0.06    []
LSalAtl2s1      161     T        0:0:0:18        0:0:0:20       -       ....    1.000   0.000   0.000   0.04    []
LSalAtl2s1      162     T        0:0:0:18        0:0:0:19       -       ....    1.000   0.000   0.000   0.02    []
LSalAtl2s1      163     G        0:0:18:0        0:0:19:0       -       ....    1.000   0.000   0.000   0.02    []

~~~~~~

*NOTE*: The actual output will depend on the options you give `varan`
Use as an example of output format only, not as a recommendation of
statistics to use!  Generally, each output option will result in a
tab-seperated column, and if the output option causes multiple values
to be output, they are separated with a space.  The exception is the
allele counts (columns four an on), which are colon-separated, similar
to the `sync`-files output by the `popoolation` tool.

Here we see that each input library (i.e. BAM file) has had the allele frequencies counted
up as the number of As, Cs, Gs, and Ts for each position.

# Examining specific regions

It is also easy to extract a specific region (e.g. a specific gene)
through `samtools`, if you index the fasta-file first.  For instance:

    samtools mpileup -f ref.fasta -r CHR8:4500-6000 sample1.bam sample2.bam | varan > out

# Specific output options for `varan`

As usual, you can specify the `--help` option to get an overview of
the available options.

~~~~~~
  -s --suppress     omit non-variant lines from output
  -v --variants     output list of non-SNP variants
  -f --fst --f-st   estimate fixation index, F_st
     --pi-k         estimate nucleotide diversity, Pi_k
  -c --conf         check if major allele frequency confidence intervals
                    overlap
     --ds           distance between major allele frequency confidence
                    intervals, using Agresti-Coull
     --dsw          lower bound for distance between major allele
                    frequencies, using Wald
  -e --esi          output conservative expected site information for SNPs
                    using Agresti-Coull intervals
     --pconf        pairwise major allele confidence
  -n --nd-all       nucleotide diversity (unadjusted), per sample and overall
  -o --output=FILE  output file name
  -g --global       calculate global statistics
  -t --threads=INT  queue lenght for parallel execution
     --min-cov=INT  minimum coverage to include
     --max-cov=INT  maximum coverage to include
  -? --help         Display help message
  -V --version      Print version information
~~~~~~

Specifying `-f` will calculate the overall F_ST.  This uses a simple
formula based on allele frequency estimates.  Another standard measure
is the expected nucleotide diversity, `--pi-k`.  The formula adjusts
for coverage, which I'm not convinced is the right thing for
sequencing data.  Interpret with caution!  On a similar note, `-n`
will calculate a simple nucleotide diversity, based on allele
frequency estimates.

The `-c` option will output a column with a '.' if the confidence
intervals for major allele frequency overlap (i.e. we can't say for
sure that they are different), a '+' if 95%-intervals are separate,
and a '*' if 99%-intervals are separate.  This can be useful to
quickly identify or count (e.g. using `fgrep -c`) candidates for
informative SNPs.  The `--pconf` option does the same, but between all
pairs of samples.

Similarly, `--ds` ouptuts the distance between confidence intervals,
and `--dsw` does the same, but using Wald's formula with pseudocounts,
which has been shown to perform well.  The latter is probably the best
choice here.

The expected information from `-e` is probably the most relevant for
selecting diagnostic SNPs.

Finally, `-g` outputs some global statistics, including a matrix of
global F_ST and nucleotide diversity between all sample pairs.

# Extracting regions with `vextr`

After filtering and sorting `varan`'s output with regular Unix tools
or what have you, you now want to extract the sequence for a
particular region, perhaps to use it for primer design.  This is one
way to do this:

    samtools mpileup -f ref.fasta -r CHR8:4500-4700 sample1.bam sample2.bam | vextr

This will simply output the sequence on a single line.  By default,
`vextr` will use IUPAC codes for ambigous bases.  If there are alternative
bases present but only at low frequency (by default singletons or less
than five percent, can be adjusted with the `--mincount` and
`--minfreq` options), `vextr` outputs the major allele, but in lower
case.

Using the `--fasta` option generates a FASTA header (using chromosome
name and starting position to construct the sequence name), and breaks
the sequence into 60 character lines.

If IUPAC is not what you want, you can specify `--format=reg` to get
regular expressions (e.g., `[A/C]`), or `--format=xs` to have ambigous
bases output as `X` characters.

# Quick and dirty visualization with `sparks`

Another way to look at mpileup data is as sparklines.  The `sparks`
executable reads from standard input.  Using `sparks disp` generates
one line of consensus sequence per sample, but highlighting the
different bases with colors, and ambigous bases with Unicode
characters so that the amount of color corresponds to base frequencies
in that position.  Alternatively `sparks info` draws a single
sparkline, showing a histogram of the information value for each
allele in each position (using the same colors as `sparks disp`).

This requires a terminal able to display Unicode graphical characters,
(e.g., `xterm` won't work, and `uxterm` requires the right font). You
can use `sparks test` to check this - it should display a series of
gradual changes.  If you see solid color blocks, you need to change
terminal and/or font.



