# quaqc metrics

quaqc generates three main categories of metrics: genome stats,
total read stats, and nuclear read stats after filtering. This last one
is the most comprehensive.

Some terminology:

* Read: Should need no explanation. Can be mapped or unmapped.
* Alignment: A read aligned to the genome. Note that only part of a
  read can be aligned, for example if the outer edges are soft-clipped.
  The distinction between read and alignment is important, because
  quaqc calculates many of the metrics using the coordinates of the
  actual alignment itself, and not those of the entire read (for example
  the length of the alignment is more useful then that of the read
  when filtering for QC).
* Read mate: When reads are paired-end (PE), then each read has a matching
  mate from each end of the sequenced fragment.
* Fragment: The initial piece of sequenced DNA, which may be longer than
  the actual read (or reads if PE). For PE experiments this information
  can be determined by taking the start coordinates of both reads in a pair.

## Target genome stats

```
┌──────────────────────────────────────────────────────────┐
│                   Target genome stats                    │
└──────────────────────────────────────────────────────────┘

  Total sequences:                             7
  ├── Nuclear sequences:                       5
  ├── Mitochondrial sequences:                 1
  └── Plastid sequences:                       1

  Total genome size:                 119,667,750
  ├── Nuclear genome size:           119,146,348 (99.6%)
  ├── Mitochondrial genome size:         366,924 (0.3%)
  └── Plastid genome size:               154,478 (0.1%)

  Effective nuclear genome size:           5,000 (0.0%)
```

This genome information is taken from the BAM header. Mitochondrial
and plastid sequences are recognized using the provided names from
the `-m/--mitochondria` and `-p/--plastids` flags, respectively.
Finally, an effective genome
size is printed: this is the actual size of the nuclear genome that
is considered by quaqc. If this number is smaller than the actual
size, then it is because any of the `-n/--target-names`, `-t/--target-list`,
or `-b/--blacklist` flags have been set, restricting the available
ranges of the genome to scan.

## Total read stats

```
┌──────────────────────────────────────────────────────────┐
│                     Total read stats                     │
└──────────────────────────────────────────────────────────┘

  Total reads:                             3,116
  ├── Mapped:                              3,114 (99.9%)
  └── Unmapped:                                2 (0.1%)

  Effective reads:                         3,114 (99.9%)
  ├┬─ Nuclear:                             3,114 (100.0%)
  │└─ Duplicated:                            605 (19.4%)
  ├┬─ Mitochondrial:                           0 (0.0%)
  │└─ Duplicated:                              0 (0.0%)
  ├┬─ Plastid:                                 0 (0.0%)
  │└─ Duplicated:                              0 (0.0%)
  ├── SE reads:                                0 (0.0%)
  ├┬─ PE reads:                            3,114 (100.0%)
  │└─ Properly mated:                      3,112 (99.9%)
  ├┬─ Primary alignments:                  3,114 (100.0%)
  │└─ Non redundant fraction:              2,509 (80.6%)
  ├── Secondary alignments:                    0 (0.0%)
  └── Supplementary alignments:                0 (0.0%)
```

This section starts with the total number of mapped and unmapped
reads in the BAM, which is information taken from the BAM index.
Afterwards, an effective read count is printed: this is the number
of mapped reads which fall within the effective nuclear genome,
as well as any mitochondria or plastid sequences if specified.
Read duplication stats are provided individually for nuclear,
mitochondrial, and plastid reads, as well as additional general
read metrics. Typically when performing plant ATAC-seq it is important
to distinguish between nuclear and plastid read duplication, since
the latter can be much higher and is not very informative of the
quality of the nuclear reads.


The non redundant fraction (NRF) refers to the number of
de-duplicated primary alignments.

PE reads are considered properly
mated if they are aligned such that they are oriented towards each
other, for example:

```
Proper:     |--->  <---|
Improper:   |--->  |--->
```

Improperly mated reads can be preserved using the `-N/--use-nomate` flag.
An additional flag can also specifically preserve those which are dovetailing
with `-D/--use-dovetails`, which is when the end coordinates of alignments
extend past the start coordinates of the mate alignment:

```
Proper:       |--->  <---|

Dovetailing:      |--->
                 <---|
```

## Nuclear read stats after filtering

```
┌──────────────────────────────────────────────────────────┐
│            Nuclear read stats after filtering            │
└──────────────────────────────────────────────────────────┘

  Reads passing filters:                   2,487 (79.9%)
```

For ATAC-seq, only nuclear reads are relevant. As a result only
these are used for the majority of the metrics.

```
  Alignment size
  ├── Min:                                    26
  ├── 1st pctile:                             36
  ├── Average:                               104.363
  ├── SD:                                     46.141
  ├── 99th pctile:                           150
  └── Max:                                   150
```

quaqc generates a histogram of the size of all alignments from reads
passing filters. (The raw data can be obtained from the JSON report.)
Some summary statistics are printed, including minimum/maximum,
1st and 99th percentiles, and average and standard deviation. Note
that the percentiles will always be printed as integers, meaning
there is a chance of slight error.

```
  Fragments passing filters:               1,245

  Fragment size
  ├── Min:                                    26
  ├── 1st pctile:                             38
  ├── Average:                               165.670
  ├── SD:                                    145.897
  ├── 99th pctile:                           776
  └── Max:                                   915
```

The same stats are also printed for fragments in the case of PE
experiments.

```
  Genome coverage:                           100.000%

  Read depth
  ├── Min:                                     0
  ├── 1st pctile:                             16
  ├── Average:                                51.910
  ├── SD:                                     15.885
  ├── 99th pctile:                           116
  └── Max:                                   123
```

The genome coverage is calculated using the effective nuclear
genome size. A read depth histogram is also generated by quaqc.

```
  Alignment MAPQ
  ├── Min:                                    31
  ├── 1st pctile:                             38
  ├── Average:                                41.879
  ├── SD:                                      0.849
  ├── 99th pctile:                            42
  ├── Max:                                    42
  └── MAPQ missing:                            0 (0.0%)
```

Similar statistic are available for alignment MAPQ scores.
MAPQ scores range from 0 to 255, with 255 signifying the
absence of a score (i.e., missing).

```
  Alignment GC
  ├── Min:                                    18%
  ├── 1st pctile:                             22%
  ├── Average:                                40.174%
  ├── SD:                                      7.167%
  ├── 99th pctile:                            56%
  └── Max:                                    62%
```

Finally, an alignment GC content histogram is generated. This
data is approximative, since as alignments can be of different
length, quaqc calculates an integer percent value between 0%
and 100%.

```
  Number of peaks:                             1 (100.0%)
  Peak coverage:                               5.000%
  Fraction of reads in peaks:                 12.505%
```

If the `-P/--peaks` flag is used, the number of alignments overlapping
the peaks will be counted. The number of peaks refers to the count
of peaks within the effective nuclear genome, and the peak coverage
refers to the fraction of the effective nuclear genome which is
covered by those peaks. Finally the fraction of all passing reads which align
inside the peaks is calculated.

```
  Number of TSSs:                              1 (100.0%)
  TSS enrichment score:                        3.046
```

Lastly, if the `-T/--tss` flag is used, a TSS pileup is generated
(available in the JSON output) and a TSS enrichment score (TES)
is provided. This is calculated as the max read density within the middle
10% of the TSS pileup range divided by the first 25% of the range (i.e.,
background levels in the promoter region).

