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
read metrics. The non redundant fraction (NRF) refers to the number of
de-duplicated primary alignments.

## Nuclear read stats after filtering

```
┌──────────────────────────────────────────────────────────┐
│            Nuclear read stats after filtering            │
└──────────────────────────────────────────────────────────┘

  Reads passing filters:                   2,487 (79.9%)
```


