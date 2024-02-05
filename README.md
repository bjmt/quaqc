# quaqc: *QU*ick *A*TAC-seq *Q*uality *C*ontrol

_Also works with any DNA-seq experiment! Compatible with plant genomes!!_

TODO:
* Create a separate quaqc-data repo that contains small test BAMs as well as Arabidopsis TSS bed.
* Create AUTHORS, NEWS files
* Double check fragment size average is correct, shouldn't it be close to the read
  size average in the ATAC-seq? I think in general I just need a suite of tests.
* Clean up the human-readable report. Center the text and use fancy wchars.
* Create a set of scripts which use jq and awk to create usable TSVs from JSON

Also benchmarks. Make sure to run these on the cluster w/ gcc as well.
* Using `ATAC_cold_MO_30.fix.bam`:
  + quaqc with --peaks,--tss and default settings: ~55 seconds, ~5.5 MB
  + samtools view -q60: ~40 seconds, ~3 MB
  + ataqv w/o tss pileup: ~1 min 45 seconds, ~62 MB
  + ataqv w/ tss pileup: ~5 min 15 seconds, ~77 MB
  + ATACseqQC? QualiMap? AIAP?

Global enrichment? (`https://ivanek.github.io/analysisOfGenomicsDataWithR/12_ATACSeq_html.html`)
Essentially just sort the read depth hist and normalize a cumsum to 0-1?

## Installation

Requires gcc/clang and GNU Make. Basic install:

TODO: Move release-full rule to top of Makefile

```sh
git clone https://github.com/bjmt/quaqc  # Or download lastest release
cd quaqc
make          # Creates quaqc binary in project folder
make install  # Optional, may require sudo
```

See [INSTALL](./INSTALL) for configuration options.

TODO: Move this to INSTALL.
### Configuration options

The following options can be added when invoking `make`.

`with_lzma=1`, `with_bz2=1`: According to the HTSlib
[installation instructions](https://github.com/samtools/htslib/blob/develop/INSTALL)
bipz2 and liblzma are only required for some CRAM files.
These options tell HTSlib to build with these libraries
when invoking `make libhts`. (Note that they must be
system-wide shared libraries, or made visible to gcc/clang
in some way.)

`with_curl=1`: By default HTSlib is built without curl
support. This causes HTSlib functions to print errors
when reading CRAM files, but as long as the accompanying
references are available locally these can be ignored.
If you need to process CRAM files which require a remote
connection, add this flag when invoking `make libhts`.
(Note that curl must be a system-wide shared library,
or made visible to gcc/clang in some way.)

`hts_dir=/path/to/HTSlib`: quaqc comes bundled with its
own version of HTSlib. To build and/or use a different
static HTSlib location, set this variable.

`z_dir=/path/to/Zlib`: quaqc comes bundled with its
own version of Zlib (specifically, the Cloudfare fork
of Zlib). To build and/or use a different static Zlib
location, set this variable.

`hts_dyn=1`: Set this to use a system-wide shared HTSlib
installation instead of a static HTSlib.

`z_dyn=1`: Set this to use a system-wide shared zlib
installation instead of a static zlib.

`native=1`: I have left out the gcc/clang
`-march=native` flag by default since older versions of clang do
not support it for Apple M1 CPUs. Use this
option to get a potential performance boost if it is available
on your system.

## Information on input BAMs

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

### Read groups

Sadly I did not consider BAMs containing multiple read groups during
the design of quaqc. This means that if a BAM contains multiple
samples, quaqc will treat them all as being one.

### CRAM

As far as I can tell CRAM files work fine as well, as long as the reference
is available locally. See the notes in the installation instructions section
about compiling quaqc with libcurl if you need to work with CRAMs depending
on remote references.

### Long reads

Some quick testing shows that quaqc seems to work with DNA-seq long read
BAMs too. However if you try it out and see obvious errors in the output,
please open an issue and I will do my best to fix the problem. Overall
this is probably the wrong tool to use to run QC on long reads.

### Spliced reads

I have not tested quaqc with BAMs containing spliced reads (e.g. RNA-seq),
but I expect most of the statistics (such as alignment and fragment lengths)
would be incorrectly calculated. Even the number of reads may be wrongly
counted, since quaqc is not set up to properly handle supplementary alignments.

## Examples usage of quaqc

TODO.

### A note on runtime and memory usage

quaqc relies entirely on HTSlib functions to access and manipulate reads
in the BAM, which are extremely efficient. In addition, quaqc is
careful to allocate the minimum required amount of memory for calculation
of the various metrics. It also makes sure to collect all necessary
data in a single pass. As a result, quaqc generally requires less than
10 MB of memory and one minute of runtime for 100 million 150 bp reads
(as tested on my MacBook Pro M1). Multi-threading is implemented very
simply, meaning it repeatedly allocates the required memory for all
independent threads, thus linearly increasing the memory usage with
every additional thread.

### Comprehensive NGS QC



### Fast iteration of read quality threshold testing



### Motif footprints



## The metrics

* Make sure to explain "effective".
* Emphasize the difference between "read" and "alignment".

### A note on additional metrics

The set of metrics output by quaqc are those that I personally have
found useful. If you can think of additional important ones,
please create an issue on GitHub and I would be happy to implement
them if it seems feasible to do so.

## Program options

TODO: Move this to doc/quaqc.md + doc/quaqc.1

Run `quaqc -h` to see an abbreviated description of all the options.

### `-m, --mitochondria [chrM,ChrM,Mt]`

Comma-separated mitochondrial reference names (no spaces). This causes quaqc
to skip reads mapping to these references for certain metrics. The remaining
metrics will be merged for all matching references under a single set
of metrics for "mitochondrial" reads. By default quaqc
checks for references with names matching to chrM, ChrM and Mt. To skip
this check entirely, simply provide an empty string (e.g. `-m''`).

### `-p, --plastids [chrC,ChrC,Pt]`

This option behaves the exact same as the `--mitochondria` option. Metrics
for all matching references will be merged under a single set of metrics
for "plastidic" reads. By default quaqc checks for references with names
matching to chrC, ChrC and Pt. To skip
this check entirely, simply provide an empty string (e.g. `-p''`).

### `--peaks [FILE.bed[.gz]]`

Filename of a BED file containing ranges that will be use to calculate
the fraction of reads in peaks (FRIP). The impact on the runtime of
quaqc when using this option is generally negligible.

### `--tss [FILE.bed[.gz]]`

Filename of a BED file containing ranges that will be used to construct
a vector of read densities around the center of the ranges.
Afterwards a final TSS enrichment score (TES) is calculated by dividing
the average read depth around the middle 10% of the ranges by the
average read depth in the first 25% of the ranges. This function is
strand-aware. If ranges are missing strand information, they are assumed
to be on the forward/positive strand. The impact on the runtime of
quaqc when using this option is generally negligible.

### `-n, --target-names [STR[,...]]`

Comma-separated reference names that quaqc will be restricted to (no spaces).
This can include reference names also provided to `--mitochondria`
and `--plastids`. The effective genome size for metric calculations
will also be adjusted. Using this option will also proportionally
scale the runtime down to the fraction of the total genome being
scanned, since quaqc will simply skip over the entire section of the
BAM not containing matching reads.

### `-t, --target-list [FILE.bed[.gz]]`

Filename of a BED file containing ranges that quaqc will be restricted
to. The effective genome size for metric calculations will also be
adjusted. Cannot be used simultaneously with `--blacklist`.

### `-b, --blacklist [FILE.bed[.gz]]`

Filename of a BED file containing ranges that quaqc will ignore. The
effective genome size for metric calculations will also be adjusted.
Cannot be used simultaneously with `--target-list`.

### `-2, --use-secondary`

Do not filter out alignments marked with the secondary alignment flag
for use in the calculation of metrics done after read filtering. Useful
when you wish to include multi-mapped reads (such as when `bowtie2` is run
with the `-k` flag).

### `--use-chimeric`

Do not filter out alignments marked with the supplementary alignment flag
for use in the calculation of metrics done after read filtering. Generally
this argument should never be set outside of special cases.

### `--use-nomate`

Do not filter out paired-end reads when it is properly mapped but
its mate read is not (for use in the calculation of metrics done
after read filtering). Note that quaqc gets fragment and mate
information from the individual read metadata and SAM flags
(i.e. it does not wait to see the corresponding read to know if it
will be properly mapped or calculate the fragment size). Depending on the
sequence aligner which created the BAM, this information may be
absent without using `samtools fixmate` or equivalent. If this is
the case, using this argument may be required to prevent quaqc from
classifying all paired-ends as filtered.

### `--use-dups`

Do not filter out duplicate reads for use in the calculation of metrics
done after read filtering. Note that quaqc can only know if reads are
duplicates based on their SAM flags, which may require use of
`samtools markdup` or equivalent.

### `-q, --mapq [30]`

Filter out reads with a MAPQ below this value for use in the calculation
of metrics done after read filtering. The default value is 30. Note that
reads with MAPQ scores of 255 are classified as unscored and
will ignore this setting.

### `--min-qlen [1]`

Filter out reads which cover fewer bases on the reference than this value
for use in the calculation of metrics done after read filtering. The
default value is 1.

### `--min-flen [1]`

Filter out paired-end reads with a fragment size below this value
for use in the calculation of metrics done after read filtering. The
default value is 1. Note that this argument is effectively ignored
when `--use-nomate` is set.

### `--max-qlen [250]`

Filter out reads which cover more bases on the reference than this value
for use in the calculation of metrics done after read filtering. The
default value is 250.

### `--max-flen [2000]`

Filter out paired-end reads with a fragment size higher than this value
for use in the calculation of metrics done after read filtering. The
default value is 2000. Note that this argument is effectively ignored
when `--use-nomate` is set.

### `--max-depth [100000]`

The maximum allowed value for the read depth histogram. Depths with
a higher value than this will be counted as this value. The default
is 100,000. Be aware that this value is used directly during memory
allocation. To estimate how much memory will be used for the depths
histogram, multiply the value by 4 bytes (e.g. the default uses 400 kB).
If quaqc detects depth values equal to the set maximum during QC
report generation, it will add a note warning that some depth values
may have been truncated (which will result in incorrect stats). The
average depth can be calculated without this data so is unaffected
by this issue.

### `--max-qhist [250]`

The maximum allowed value for the covered bases per read histogram.
Density values greater than this value will be reported as this value.
The default is the value of `--max-qlen`. Be aware that this value is
used directly during memory
allocation. To estimate how much memory will be used for the read size
histogram, multiply the value by 4 bytes (e.g. the default uses 1 kB).
If quaqc detects size values equal to the set maximum during QC
report generation, it will add a note warning that some size values
may have been truncated (which will result in incorrect stats). The
average size can be calculated without this data so is unaffected
by this issue.

### `--max-fhist [2000]`

The maximum allowed value for the fragment size histogram.
Density values greater than this value will be reported as this value.
The default is the value of `--max-flen`. Be aware that this value is
used directly during memory
allocation. To estimate how much memory will be used for the fragment size
histogram, multiply the value by 4 bytes (e.g. the default uses 8 kB).
If quaqc detects size values equal to the set maximum during QC
report generation, it will add a note warning that some size values
may have been truncated (which will result in incorrect stats). The
average size can be calculated without this data so is unaffected
by this issue.

### `--tss-size [2000]`

The size of the density vector range generated when `--tss` is set, in bases.
Ranges are first centered at their midpoints, then resized in both
directions to a final width of the set value. The default is 2000.
Be aware that this value is used directly during memory
allocation. To estimate how much memory will be used for the TSS density
values, multiply the value by 4 bytes (e.g. the default uses 8 kB).

### `--tss-qlen [100]`

The final size of adjusted read coordinates when generating the read
density values when `--tss` is set. Reads are first set to size 1 (anchored
from their five-prime ends), then resized in both directions to a final
width of the set value. The default is 100. To prevent read resizing and
instead use the actual coordinates of the reads, set this value to 0.

### `--tss-tn5`

When resizing the reads as described in the `--tss-qlen` option, adjust
the read five-prime coordinates forward 4 bases (to center the coordinate
in the middle of the Tn5 transposase binding site). This option is ignored
when `--tss-qlen` is set to 0.

### `-k, --keep`

Save the nuclear reads passing all filters in a new BAM. This will
significantly slow down quaqc.

### `--keep-dir [DIR/]`

By default, when `--keep-filtered` is set a new filtered BAM is
created in the same directory as the input BAM. Setting this will
change the final directory where the new BAM will be written.

### `--keep-ext [.filt.bam]`

By default, when `--keep-filtered` is set a new filtered BAM is
created with the text `.filt.bam` appended to the file name. Use this
argument to change it. If an existing `.bam` or `.cram` extension
exists, it will be stripped.

### `--omit-gc`

Skip GC content metrics. This can shave off a small percentage of the
runtime for regular short read experiments (<10%). The savings may be
more substantial for long read experiments, as quaqc has to iterate
over every base in the alignments to count GC bases.

### `--omit-depth`

Skip generation of the read depths histogram. This can shave off a
small percentage of the runtime for regular short read experiments (<10%).
The savings may be more substantial for long read experiments, as quaqc
has to iterate over the entire alignment length to count per-base depths.

### `-f, --fast`

Set both `--omit-gc` and `--omit-depth`, thus skipping the two metric which
require iterating over the entire read lengths. Together this can shave off
about 15% of the runtime for regular short read experiments. The savings
may be more substantial for long read experiments.

### `-j, --threads [1]`

Set the number of child threads used to process input BAMs. At minimum,
one child thread is launched (meaning quaqc technically uses two threads,
though not simultaneously), and at maximum, one child thread per sample
is launched (in addition to the main parent thread). All of the data
structures are duplicated for each new thread, meaning memory usage will
increase linearly with increasing thread count. When using default settings,
the `--max-depth` option has the biggest impact on memory growth. Set this
to a lower value to mitigate this.

### `-c, --continue`

If set when processing more than one input file, quaqc will keep running if
it encounters errors processing individual files (e.g. one file is unsorted).

## Misc

Note that all BED inputs can also be gzipped.

It is possible to get an odd number of passing PE alignments even when not using --use-nomate.
This is because quaqc only can look at the mate info encoded into the R1 of each
fragment, meaning that if it suggests the fragment is ok it will pass, even if
the R2 does not pass filters.

The fragment size information in SAM/BAM is implementation-defined when reads
overlap as such:
      |-------------->
<---------|
For example, my ATAC BAM seems to only count the overlapping bases in the middle.

Create a mini json parser utility to extract certain fields from json qc,
such as the TSS plot etc. Also can be used to merge multiple reports.

Have a section in the README explaining the json format.

Make an R shiny script to visualize results?

Seems to work just fine with CRAM. It complains about not being able
to access references hosted on www.ebi.ac.uk/ena, but this does not appear
to actually prevent it from working properly. To enable curl again,
just delete the configure --disable-libcurl line in the Makefile.


lldb howto:

  lldb bin/quaqc
  breakpoint set -f quaqc.c -l ###
  run test/ex1.bam
  ..
  v
  bt
  frame var -T ...


