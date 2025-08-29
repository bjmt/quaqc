[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/quaqc/README.html)
# quaqc: *QU*ick *A*TAC-seq *Q*uality *C*ontrol

_Also works with any unspliced DNA-seq experiment! Compatible with plant genomes!!_

## Installation

Requires gcc/clang and GNU Make, tested with macOS and Linux. Basic install:

```sh
git clone https://github.com/bjmt/quaqc  # Or download latest release
cd quaqc
make release-full
make test  # Make sure quaqc produces the expected outputs (optional)
make install  # Copy quaqc + manual to /usr/local (optional, may require sudo)
```

See [INSTALL](./INSTALL) for configuration options.

Alternatively for conda users, grab it from [Bioconda](https://bioconda.github.io/recipes/quaqc/README.html):

```sh
# Example environment setup:
conda create -y --name quaqc -c conda-forge -c bioconda quaqc
```

## Quick start

To see a basic description of all available parameters, run `./quaqc -h`. For
more detailed help, see the [manual page](./doc/quaqc.1.md).

quaqc comes built-in with sensible defaults, meaning few parameters must be
specified for regular usage. (The hardcoded defaults can be changed before
compilation by editing [quaqc.h](src/quaqc.h).) To run quaqc using the test
data:

```sh
./quaqc -v --output-dir . --peaks test/peak.bed --tss test/tss.bed \
    --target-list test/target.bed test/reads.bam
cat ./reads.quaqc.txt
```

This command changes the location of the output QC report to the current directory,
and uses the peak and TSS files to generate additional stats. Since this example
BAM is just a subset of reads contained within a small region of the Arabidopsis
genome, the `target.bed`
file restricts quaqc to only consider that region of the genome (and thus adjust
how it calculates the stats). See a copy of the output
report [here](./test/reads.quaqc.txt). 

Note that quaqc can also save the final reads passing all filters to a new BAM
file via the `-S/--keep` flag. (This will increase the runtime.)

## quaqcr

For those who prefer to analyze their data with R, a companion package
[`quaqcr`](https://github.com/bjmt/quaqcr) is available. This package contains
functions for both executing quaqc from within R as well as analyzing all of
its outputs made available when outputting the results in JSON format. See
the README for a detailed tutorial.

## Help

Please open an [issue](https://github.com/bjmt/quaqc/issues) on GitHub.

## Citation

Please cite the following journal article if you find this software useful
in your research:

Tremblay, B.J.M. and Qüesta, J.I. (2024). quaqc: efficient and quick ATAC-seq quality control and filtering. _Bioinformatics_ **40**, btae649. [https://doi.org/10.1093/bioinformatics/btae649](https://doi.org/10.1093/bioinformatics/btae649)

## Additional use cases

### Creating insertion BED files for use by MACS2/3

A common method for calling ATAC-seq peaks is to provide the reads as a
BED file in order to trick MACS2/3 into using individual reads during pileup
instead of the fragment coordinates in paired-end data. Using the `--bed`
option, this can be done using quaqc, generating identical output as
`bedtools bamtobed`. For example:

```sh
quaqc -0 --bed Sample.bam

macs3 callpeak [...] \
    --treatment Sample.bed.gz --format BED \
    --call-summits --keep-dup all \
    --shift -75 --extsize 150 \
    --nomodel --nolambda
```

However since we only care about the 5-prime ends of each read (i.e. the
insertion sites) when calling peaks this way, we can also add the `--bed-ins`
flag to tell quaqc to only output the 5-prime coordinates of the reads with
no additional BED fields, saving a decent amount of disk space in a way
which will not affect MACS2/3.

### Creating bedGraph files centered around the 5-prime ends of reads

Typically read pileups of ATAC-seq data (such as bedGraph files) are made
in such a way as to visualize the sites of transposition, which are the
5-prime ends of the reads. This means needing to resize each read before
creating a bedGraph file. Both of these steps can be performed simultaneously
using quaqc:

```sh
quaqc --bedGraph --bedGraph-qlen 150 Sample.bam
```

In this example, a bedGraph file is created where each read is resized
to 150 bp, centered around the 5-prime end of each read. (To use the original read sizes without resizing, simply set `--bedGraph-qlen` to 0.) This functionality can also be used to visualize the exact insertions, adjusting for the +4/-5 transposition shift:

```sh
quaqc --bedGraph --bedGraph-tn5 --bedGraph-qlen 1 Sample.bam
```

The +4/-5 shift can also be adjusted to use custom values via `--tn5-fwd` and `--tn5-rev`.

### Filtering a single-cell ATAC-seq BAM by cell type

Using a set of barcodes, quaqc can subset a BAM by making use of the `--rg-list` and
`--keep` options. For example, using the common `CB` SAM flag typically used for storing
barcodes:

```sh
quaqc -0 -A --rg-list barcodes.txt --rg-tag CB --keep --keep-ext .MyCellType.bam Sample.bam
```

The `-0` flag prevents quaqc from creating a QC report, thus making the output BAM
the only file created by quaqc. The `-A` flag can be used if this BAM has already been
filtered, and no additional filtering is desired.

### Fast iteration of the effects of different read quality thresholds

quaqc is fairly performant, but it can still require some patience for decently
sized BAMs (e.g., about one minute for 80 million reads from Arabidopsis ATAC-seq).
To substantially speed up quaqc to iterate through different filtering
thresholds, make use of the `-n/--target-names` to restrict quaqc to only
look at single chromosomes. This is often enough to get an idea of the quality of
the data while getting a several fold speed up (e.g., restricting quaqc to the
first chromosome of Arabidopsis takes just over 10 seconds for the previously
mentioned example).

For example, to try out several MAPQ thresholds for an Arabidopsis ATAC-seq sample
by only scanning chromosome 1 (which is just '1' for this version of TAIR10):

```sh
for i in 5 10 15 20 25 30 ; do
    quaqc -n 1 -q $i -i "MAPQ >= ${i}" -O .quaqc.mapq${i}.txt Sample.bam
done
```

Now, six different MAPQ thresholds have been tested in the time it would take
to process the entire BAM file once across all chromosomes. (The `-i/--title`
flag is used to add a unique title to each run.)

### Quick motif footprints

The TSS pileup functionality of quaqc can be used to instead generate motif
footprints. Note that these are only useful for testing, as no correction
of base composition is peformed to reduce Tn5 transposase insertion biases. The output
must be saved in JSON format to recover the footprint data, which can then
be plotted with an external program (such as [`quaqcr`](https://github.com/bjmt/quaqcr)).

In this example code, the `--footprint` preset will adjust the TSS pileup
to produce single base resolution data of Tn5 transposase insertion
frequency. The `--nfr` preset is also used to only consider nucleosome-free
reads.

```sh
quaqc --nfr --footprint --tss motif.bed --json sample_motif.json Sample.bam
```

## The metrics

See this [annotated report](doc/metrics.md) for a description of the metrics.

### A note on additional metrics

The set of metrics output by quaqc are those that I personally have
found useful. If you can think of additional important ones,
please create an issue on GitHub and I would be happy to implement
them if it seems feasible to do so.

## Extra

### A note on the default Tn5 shift values (+4/-5)

It is standard practice to adjust the start coordinates of ATAC-seq
reads to account for the ~9 bp footprint of the Tn5 transposase by
adding 4 bp to the start and substracting 5 bp from the end (+4/-5).
However, some articles such as [this one](https://www.biorxiv.org/content/10.1101/2024.12.25.630221v2.full)
instead suggest +4/-4 might better capture the true Tn5 offset and allow
for better estimation of the sequence composition bias of transposition.
To use this alternative offset in quaqc, simply add `--tn5-rev 4`.

### Information on input BAMs

quaqc is meant to work with unfiltered BAMs produced by sequence
alignment tools such as `bowtie2`. However, some processing
steps are recommended and/or required. These include:

- `samtools fixmate` so that quaqc can compute accurate
  mate and fragment size metrics for PE experiments (optional,
  though not having proper read mate data may require using
  the `--use-nomate` flag to stop quaqc from filtering all PE reads)
- `samtools sort` if the reads are not already
  coordinate sorted (required)
- `samtools markdup` (or an equivalent command) so that
  quaqc can compute read duplication metrics (optional)
- `samtools index` to generate an accompanying index file
  (optional; quaqc will create one if missing)

My personal recommendation is to run all of these commands simultaneously
during sequence alignment. For example:

```sh
[bowtie2/etc ... ] \
  | samtools fixmate -m -u - - \
  | samtools sort -u - \
  | samtools markdup --write-index - file.bam
```

This requires a sequence aligner which can output its results to `stdout`,
as well as being grouped by read name. For further information check out
[this article](http://www.htslib.org/algorithms/duplicate.html).

quaqc only works with indexed BAMs, though actually indexing them before
using quaqc is optional. If no index can be found then quaqc will
attempt to create the index.

### CRAM

As far as I can tell CRAM files work fine as well, as long as the reference
is available locally. See the notes in the installation instructions
about compiling quaqc with libcurl if you need to work with CRAMs depending
on remote references. The Bioconda version is dynamically linked to the full version of HTSlib, making it fully compatible with any type of CRAM.

### Long reads

Some quick testing shows that quaqc seems to work with DNA-seq long read
BAMs too. However if you try it out and see obvious errors in the output,
please open an issue and I will do my best to fix the problem. Overall
this is probably the wrong tool to use to run QC on long reads. Keep in
mind quaqc will not properly handle supplementary alignments.

### Spliced reads

I have not tested quaqc with BAMs containing spliced reads (e.g. RNA-seq),
but I expect most of the statistics (such as alignment and fragment lengths)
would be incorrectly calculated. Even the number of reads may be wrongly
counted, since quaqc is not set up to properly handle supplementary alignments.

