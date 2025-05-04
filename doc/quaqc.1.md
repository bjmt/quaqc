% QUAQC(1) Version 1.3a | quaqc User Manual

# NAME

quaqc - Quick ATAC-seq QC

# SYNOPSIS

quaqc \[options] file1.bam \[file2.bam \[...]]

# DESCRIPTION

quaqc allows for ATAC-seq-specific quality control and read filtering of NGS
data with minimal processing time and extremely low memory overhead. Any number
of samples can be processed, using multiple threads if desired. quaqc outputs
a comprehensive set of aligned read metrics, including alignment size, fragment
size, percent duplicates, mapq scores, read depth, GC content, and others.
Although designed for ATAC-seq data, quaqc can also be used for other unspliced
DNA-seq experiments (such as ChIP-seq) as many of the metrics are related to
genereal sequencing quality.  Non-nuclear aligned reads can also be separately
considered if so desired, helpful when determining the level of chloroplast
contamination in plant ATAC-seq.

When aligning ATAC-seq reads, it is recommended not to perform any filtering
before using quaqc, as well as marking duplicates and fixing mate coordinates
for paired end data. The BAM must also be coordinate sorted. This can all be
accomplished using `samtools` in the following manner during alignment:

    [bowtie2/etc ... ] \
        | samtools fixmate -m -u - - \
        | samtools sort -u - \
        | samtools markdup --write-index - file.bam

# OPTIONS

**-m** *STR*, **\--mitochondria**=*STR*
:   Comma-separated mitochondrial reference names (no spaces). This causes quaqc
    to skip reads mapping to these references for certain metrics. The remaining
    metrics will be merged for all matching references under a single set
    of metrics for "mitochondrial" reads. By default quaqc
    checks for references with names matching to chrM, ChrM, Mt and mitochondria.
    To skip this check entirely, simply provide an empty string (e.g. `-m ' '`).

**-p** *STR*, **\--plastids**=*STR*
:   This option behaves the exact same as the `--mitochondria` option. Metrics
    for all matching references will be merged under a single set of metrics
    for "plastidic" reads. By default quaqc checks for references with names
    matching to chrC, ChrC, Pt and chloroplast. To skip
    this check entirely, simply provide an empty string (e.g. `-p ' '`).

**-P** *FILE*, **\--peaks**=*FILE*
:   Filename of a BED file containing ranges that will be use to calculate
    the fraction of reads in peaks (FRIP). The impact on the runtime of
    quaqc when using this option is generally negligible. Can be gzipped.

**-T** *FILE*, **\--tss**=*FILE*
:   Filename of a BED file containing ranges that will be used to construct
    a vector of read densities around the center of the ranges.
    Afterwards a final TSS enrichment score (TES) is calculated by dividing
    the average read depth around the middle 10% of the ranges by the
    average read depth in the first 25% of the ranges. This function is
    strand-aware. If ranges are missing strand information, they are assumed
    to be on the forward/positive strand. The impact on the runtime of
    quaqc when using this option is generally negligible. Can be gzipped.

**-n** *STR*, **\--target-names**=*STR*
:   Comma-separated reference names that quaqc will be restricted to (no spaces).
    This can include reference names also provided to `--mitochondria`
    and `--plastids`. The effective genome size for metric calculations
    will also be adjusted. Using this option will also proportionally
    scale the runtime down to the fraction of the total genome being
    scanned, since quaqc will simply skip over the entire section of the
    BAM not containing matching reads.

**-t** *FILE*, **\--target-list**=*FILE*
:   Filename of a BED file containing ranges that quaqc will be restricted
    to. The effective genome size for metric calculations will also be
    adjusted. Cannot be used simultaneously with `--blacklist`. Can be gzipped.
    To be used for QC, reads must be fully contained within these ranges.

**-b** *FILE*, **\--blacklist**=*FILE*
:   Filename of a BED file containing ranges in which all read aligning to
    these regions are to be ignored by quaqc. The effective genome size for
    metric calculations will also be adjusted. Cannot be used simultaneously
    with `--target-list`. Can be gzipped. To be used for QC, reads must be
    fully outside these ranges.

**-r** *STR*, **\--rg-names**=*STR*
:   Comma-separated read group identifiers. This will cause quaqc to discard
    reads which do not have matching read group tags. Cannot be used
    simultaneously with `--rg-list`. If a BAM header does not contain any of
    the target read groups, the run will end with an error.

**-R** *STR*, **\--rg-list**=*STR*
:   Filename of a file containing read group identifiers, one per line. This
    will cause quaqc to discard reads which do not have a matching read group
    tag. Cannot be used simultaneously with `--rg-names`. If a BAM header does
    not contain any of the target read groups, the run will end with an error.
    Can be gzipped.

**\--rg-tag**=*STR*
:   Filter reads using a different SAM tag than the read group (RG) tag. This
    can be used, for example, to filter single-cell BAMs using the barcode (CB)
    tag. Please note that the SAM spec requires all possible RG groups to be
    encoded in the header, meaning quaqc can check that the provided groups
    from `--rg-list` exist in the BAM before starting. No such check is performed
    when using a different tag, meaning mistakes can lead to quaqc filtering
    out all reads.

**-2**, **\--use-secondary**
:   Do not filter out alignments marked with the secondary alignment flag
    for use in the calculation of metrics done after read filtering. Useful
    when you wish to include multi-mapped reads (such as when `bowtie2` is run
    with the `-k` flag).

**-N**, **\--use-nomate**
:   Do not filter out paired-end reads when it is properly mapped but
    its mate read is not (for use in the calculation of metrics done
    after read filtering). Note that quaqc gets fragment and mate
    information from the individual read metadata and SAM flags
    (i.e. it does not wait to see the corresponding read to know if it
    will be properly mapped or calculate the fragment size). Depending on the
    sequence aligner which created the BAM, this information may be
    absent without using `samtools fixmate` or equivalent. If this is
    the case, using this argument may be required to prevent quaqc from
    classifying all paired-ends as filtered.

**-d**, **\--use-dups**
:   Do not filter out duplicate reads for use in the calculation of metrics
    done after read filtering. Note that quaqc can only know if reads are
    duplicates based on their SAM flags, which may require use of
    `samtools markdup` or equivalent.

**\--use-chimeric**
:   Do not filter out alignments marked with the supplementary alignment flag
    for use in the calculation of metrics done after read filtering. Generally
    this argument should never be set outside of special cases.

**-D**, **\--use-dovetails**
:   Do not filter out dovetailing alignments from paired reads. Dovetailing
    is when a read end coordinate extends past the start coordinate of its mate.

**\--no-se**
:   Filter out single end (SE) reads, if a BAM happens to contain both SE and
    PE reads.

**-q** *INT*, **\--mapq**=*INT*
:   Filter out reads with a MAPQ below this value for use in the calculation
    of metrics done after read filtering. The default value is 30. Note that
    reads with MAPQ scores of 255 are classified as unscored and
    will ignore this setting.

**\--min-qlen**=*INT*
:   Filter out reads which cover fewer bases on the reference than this value
    for use in the calculation of metrics done after read filtering. The
    default value is 1.

**\--min-flen**=*INT*
:   Filter out paired-end reads with a fragment size below this value
    for use in the calculation of metrics done after read filtering. The
    default value is 1. Note that this argument is effectively ignored
    when `--use-nomate` is set.

**\--max-qlen**=*INT*
:   Filter out reads which cover more bases on the reference than this value
    for use in the calculation of metrics done after read filtering. The
    default value is 250.

**\--max-flen**=*INT*
:   Filter out paired-end reads with a fragment size higher than this value
    for use in the calculation of metrics done after read filtering. The
    default value is 2000. Note that this argument is effectively ignored
    when `--use-nomate` is set.

**-A**, **\--use-all**
:   Do not apply ANY filters to mapped reads. Note that unmapped reads
    are always ignored, no matter which options are set.

**\--max-depth**=*INT*
:   The maximum allowed value for the read depth histogram. Depths with
    a higher value than this will be counted as this value. The default
    is 100,000. Be aware that this value is used directly during memory
    allocation. To estimate how much memory will be used for the depths
    histogram, multiply the value by 4 bytes (e.g. the default uses 400 kB).
    If quaqc detects depth values equal to the set maximum during QC
    report generation, it will add a note warning that some depth values
    may have been truncated (which will result in incorrect stats). The
    average depth can be calculated without this data so is unaffected
    by this issue.

**\--max-qhist**=*INT*
:   The maximum allowed value for the covered bases per read histogram.
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

**\--max-fhist**=*INT*
:   The maximum allowed value for the fragment size histogram.
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

**\--tss-size**=*INT*
:   The size of the density vector range generated when `--tss` is set, in bases.
    Ranges are first centered at their midpoints, then resized in both
    directions to a final width of the set value. The default is 2000.
    Be aware that this value is used directly during memory
    allocation. To estimate how much memory will be used for the TSS density
    values, multiply the value by 4 bytes (e.g. the default uses 8 kB).

**\--tss-qlen**=*INT*
:   The final size of adjusted read coordinates when generating the read
    density values when `--tss` is set. Reads are first set to size 1 (anchored
    from their five-prime ends), then resized in both directions to a final
    width of the set value. The default is 100. To prevent read resizing and
    instead use the actual coordinates of the reads, set this value to 0.

**\--tss-tn5**
:   When resizing the reads as described in the `--tss-qlen` option, adjust
    the read five-prime coordinates forward 4 bases (to center the coordinate
    in the middle of the Tn5 transposase binding site). This option is ignored
    when `--tss-qlen` is set to 0.

**\--omit-gc**
:   Skip GC content metrics. This can shave off a small percentage of the
    runtime for regular short read experiments (<10%). The savings may be
    more substantial for long read experiments, as quaqc has to iterate
    over every base in the alignments to count GC bases.

**\--omit-depth**
:   Skip generation of the read depths histogram. This can shave off a
    small percentage of the runtime for regular short read experiments (<10%).
    The savings may be more substantial for long read experiments, as quaqc
    has to iterate over the entire alignment length to count per-base depths.

**-f**, **\--fast**
:   Set `--omit-gc` and `--omit-depth`, thus skipping the two metric which
    require iterating over the entire read lengths. Together this can shave off
    about 15% of the runtime for regular short read experiments. The savings
    may be more substantial for long read experiments.

**\--lenient**
:   Set `--use-nomate`, `--use-dups`, `--use-dovetails`, and `--mapq=10`. This
    relaxes the filtering parameters, allowing a greater number of reads to be
    counted for QC.

**\--strict**
:   Set `--min-flen=50`, `--max-flen=150`, and `--mapq=40`. This
    restricts the filtering parameters to keep only the highest quality
    reads.

**\--nfr**
:   Set `--no-se`, `--max-flen=120`, and `--tss-tn5`. These filters enrich for
    reads found within nucleosome free regions (NFR), as well as shifting the
    start sites to account for the Tn5 transposase insertion.

**\--nbr**
:   Set `--no-se`, `--min-flen=150`, `--max-flen=1000`, and `--tss-qlen=0`.
    These filters enrich for reads in nucleosome bound regions (NBR). In
    addition, the read coordinates are maintained as is when generating the
    TSS pileup.

**\--footprint**
:   Set `--tss-qlen=1`, `--tss-size=501`, and `--tss-tn5`. This generates a
    smaller TSS pileup with single base pair resolution  of Tn5 transposase
    insertion frequency.

**\--chip**
:   Set `--tss-qlen=0` and `--tss-size=5001`. Additionally, any BED file
    provided with the `--peaks` option is used for generating the pileup
    (which is normally generated from `--tss`).

**-o** *DIR*, **\--output-dir**=*DIR*
:   Directory where the QC reports will be saved. By default, these are
    saved in the same directory as the input BAMs.

**-O** *STR*, **\--output-ext**=*STR*
:   Filename extension of the QC report, replacing the previous ".bam" of
    the input BAMs. By default ".quaqc.txt" is used.

**-0**, **\--no-output**
:   Suppress the generation of QC reports.

**-J** *FILE*, **\--json**=*FILE*
:   Save all QC reports for all samples into a single JSON file for further
    processing. This format, while not intended to be human readable, contains
    additional data such as the full alignment size, fragment size, GC content,
    mapq, and read depth histograms, as well as the TSS pileup. To save to
    standard output, provide `-J-`. To compress the output JSON, add the ".gz"
    extension to the supplied filename.

**-S**, **\--keep**
:   Save the nuclear reads passing all filters in a new BAM. This will
    significantly slow down quaqc.

**-k** *DIR*, **\--keep-dir**=*DIR*
:   By default, when `--keep` is set a new filtered BAM is
    created in the same directory as the input BAM. Setting this will
    change the final directory where the new BAM will be written.

**-K** *STR*, **\--keep-ext**=*STR*
:   By default, when `--keep` is set a new filtered BAM is
    created with the text ".filt.bam" appended to the file name. Use this
    argument to change it. If an existing ".bam" or ".cram" extension
    exists, it will be stripped.

**-B**, **\--bedGraph**
:   Output a Gzipped bedGraph of the alignments passing all filters within
    target regions. Only limited memory is used to store bedGraph records
    which are output on the fly. This means that some consecutive positions
    with identical scores will sometimes be present as distinct ranges if
    quaqc ran out of memory to store such records simultaneously. 

**\--bedGraph-qlen**=*INT*
:   When outputting the alignments in bedGraph format, they are resized from
    the 5-prime position. This option controls the final size of the 5-prime
    centered alignment. To instead use the original start and end coordinates
    of the alignment, set this option to 0.

**\--bedGraph-tn5**
:   Shift the 5-prime alignment coordinates of each read when generating the
    bedGraph to account for the transposase offset (+4/-5), as per the
    `--tss-tn5` option.

**\--bedGraph-dir**=*DIR*
:   As per the `--keep-dir` option, change the output directory of the 
    bedGraph files.

**\--bedGraph-ext**=*STR*
:   As per the `--keep-ext` option, change the default bedGraph filename
    extension. The bedGraph will always be Gzipped, so not including ".gz"
    will still lead to a compressed file.

**\--tn5-fwd**=*INT*
:   Alter the default value used to shift the 5-prime ends of reads when either
    `--tss-tn5` or `--bedGraph-tn5` are used. This value is added to the start
    coordinate of forward strand reads.

**\--tn5-rev**=*INT*
:   Alter the default value used to shift the 5-prime ends of reads when either
    `--tss-tn5` or `--bedGraph-tn5` are used. This value is substracted from the
    end coordinate of reverse strand reads.

**-j** *INT*, **\--threads**=*INT*
:   Set the number of child threads used to process input BAMs. At minimum,
    one child thread is launched (meaning quaqc technically uses two threads,
    though not simultaneously), and at maximum, one child thread per sample
    is launched (in addition to the main parent thread). All of the data
    structures are duplicated for each new thread, meaning memory usage will
    increase linearly with increasing thread count. When using default settings,
    the `--max-depth` option has the biggest impact on memory growth. Set this
    to a lower value to mitigate this.

**-i** *STR*, **\--title**=*STR*
:   Assign a title to the run. All output reports will contain this title.

**-c**, **\--continue**
:   If set when processing more than one input file, quaqc will keep running if
    it encounters errors processing individual files (e.g. one file is unsorted).

**-v**, **\--verbose**
:   Print progress messages during runtime. This flag can be used a second time
    to further increase verbosity.

**\--version**
:   Print the version number of quaqc to `stdout` and exit.

**-h**, **\--help**
:   Print a help message with a brief description of all available commands.

# BUGS

Please report bugs on GitHub: <https://github.com/bjmt/quaqc/issues>

# AUTHOR

quaqc was created by Benjamin Jean-Marie Tremblay.

