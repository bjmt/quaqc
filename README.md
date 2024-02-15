# quaqc: *QU*ick *A*TAC-seq *Q*uality *C*ontrol

_Also works with any unspliced DNA-seq experiment! Compatible with plant genomes!!_

## Installation

Requires gcc/clang and GNU Make, tested with macOS and Linux. Basic install:

```sh
git clone https://github.com/bjmt/quaqc  # Or download latest release
cd quaqc
make release-full
make install  # Copy quaqc + manual to /usr/local (optional, may require sudo)
make test  # Make sure quaqc produces the expected outputs (optional)
```

See [INSTALL](./INSTALL) for configuration options.

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
BAM is just a subset of reads contained within a small region of the Arabidpsis
genome, the `target.bed`
file restricts quaqc to only consider that region of the genome (and thus adjust
how it calculates the stats). 

Note that quaqc can also save the final reads passing all filters to a new BAM
file via the `-S/--keep` flag. (This will increase the runtime.)

## Help

Please open an [issue](https://github.com/bjmt/quaqc/issues) on GitHub or email
me at `benjmtremblay` _at_ `gmail.com`.

## Citation

No article about quaqc is yet published. Until then, please cite this repository
if you find this software useful in your research.

Tremblay, B.J.M. (2024). quaqc: Quick ATAC-seq QC (Version 0.1) [Computer software]. https://github.com/bjmt/quaqc

## Additional use cases

### Fast iteration of the effects of different read quality thresholds

quaqc is fairly performant, but it can still require some patience for decently
sized BAMs (e.g., about one minute for 40 million reads from Arabidopsis ATAC-seq).
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
to process the entire BAM file once across all chromosomes.

### Quick motif footprints

The TSS pileup functionality of quaqc can be used to instead generate motif
footprints. Note that these are only useful for testing, as no correction
of base composition is peformed to reduce Tn5 insertion biases. The output
must be saved in JSON format to recover the footprint data, which can then
be plotted with an external program.

In this example code, the `--footprint` preset will adjust the TSS pileup
to produce single base resolution data of Tn5 transposase insertion
frequency. The `--nfr` present is also used to only consider nucleosome-free
reads.

```sh
quaqc --nfr --footprint --tss motif.bed --json sample_motif.json Sample.bam
```

## The metrics

* Make sure to explain "effective".
* Emphasize the difference between "read" and "alignment".

### A note on additional metrics

The set of metrics output by quaqc are those that I personally have
found useful. If you can think of additional important ones,
please create an issue on GitHub and I would be happy to implement
them if it seems feasible to do so.

## Extra

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

### CRAM

As far as I can tell CRAM files work fine as well, as long as the reference
is available locally. See the notes in the installation instructions
about compiling quaqc with libcurl if you need to work with CRAMs depending
on remote references.

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

