#!/bin/bash
# quaqc comprehensive test suite.
# Requires samtools and python3 (standard library only) on PATH.
# Run from the test/ directory: bash test_suite.sh

set -uo pipefail

QUAQC=$(cd "$(dirname "$0")/.." && pwd)/quaqc
TMPROOT=$(mktemp -d /tmp/quaqc_suite_XXXXX)
PASS=0; FAIL=0

trap 'rm -rf "$TMPROOT"' EXIT

# ---------------------------------------------------------------------------
# Test framework
# ---------------------------------------------------------------------------

_cur_test=""

run_test() {
  _cur_test="$1"
  local d="$TMPROOT/$1"
  mkdir -p "$d"
  printf 'Testing %-50s ' "$1 ..."
  if "test_$1" "$d" >"$d/stdout" 2>"$d/stderr"; then
    echo "PASS"
    PASS=$((PASS + 1))
  else
    echo "FAIL"
    echo "    stderr: $(tail -3 "$d/stderr" | tr '\n' ' ')"
    FAIL=$((FAIL + 1))
  fi
}

# Assert helpers used inside test functions (return 1 on failure).
assert_eq() {
  local label="$1" actual="$2" expected="$3"
  if [ "$actual" = "$expected" ]; then return 0; fi
  echo "  assert_eq $label: expected='$expected' got='$actual'" >&2
  return 1
}

assert_contains() {
  local file="$1" pattern="$2"
  if grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_contains: '$pattern' not found in $file" >&2
  return 1
}

assert_not_contains() {
  local file="$1" pattern="$2"
  if ! grep -qF "$pattern" "$file"; then return 0; fi
  echo "  assert_not_contains: '$pattern' unexpectedly found in $file" >&2
  return 1
}

assert_file_exists() {
  if [ -s "$1" ]; then return 0; fi
  echo "  assert_file_exists: '$1' missing or empty" >&2
  return 1
}

# ---------------------------------------------------------------------------
# BAM construction helpers
# ---------------------------------------------------------------------------

# seq N char — print N copies of char
_seq() { python3 -c "import sys; print(sys.argv[2]*int(sys.argv[1]), end='')" "$1" "$2"; }

# pe_pair name chr r1_pos r1_len r2_pos r2_len mapq [extra_flag_r1] [extra_flag_r2]
# Outputs two SAM lines for a proper PE pair (R1 fwd, R2 rev).
# Positions are 1-based (SAM convention).
pe_pair() {
  local name="$1" chr="$2" r1p="$3" r1l="$4" r2p="$5" r2l="$6" mapq="${7:-42}"
  local xf1="${8:-0}" xf2="${9:-0}"
  local tlen=$((r2p + r2l - r1p))
  local f1=$(( 99 + xf1 ))    # 1+2+32+64 = paired+proper+mate_rev+r1
  local f2=$(( 147 + xf2 ))   # 1+2+16+128 = paired+proper+rev+r2
  local s1 s2
  s1=$(_seq "$r1l" A)
  s2=$(_seq "$r2l" T)
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%s\t*\n' \
    "$name" "$f1" "$chr" "$r1p" "$mapq" "$r1l" "$r2p" "$tlen" "$s1"
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t-%d\t%s\t*\n' \
    "$name" "$f2" "$chr" "$r2p" "$mapq" "$r2l" "$r1p" "$tlen" "$s2"
}

# se_read name chr pos len mapq [flag]
se_read() {
  local name="$1" chr="$2" pos="$3" len="$4" mapq="${5:-42}" flag="${6:-0}"
  local s; s=$(_seq "$len" A)
  printf '%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t*\n' \
    "$name" "$flag" "$chr" "$pos" "$mapq" "$len" "$s"
}

# dup_pair — same as pe_pair but with duplicate bit (0x400) set on both reads
dup_pair() {
  pe_pair "$1" "$2" "$3" "$4" "$5" "$6" "${7:-42}" 1024 1024
}

# secondary_read name chr pos len mapq — secondary alignment (FLAG=256)
secondary_read() {
  se_read "$1" "$2" "$3" "$4" "${5:-42}" 256
}

# make_bam outfile — reads SAM from stdin, sorts, indexes
make_bam() {
  local out="$1"
  samtools sort -o "$out" - 2>/dev/null
  samtools index "$out" 2>/dev/null
}

# chr_header name len [name len ...]
chr_header() {
  printf '@HD\tVN:1.6\tSO:queryname\n'
  while [ $# -ge 2 ]; do
    printf '@SQ\tSN:%s\tLN:%d\n' "$1" "$2"
    shift 2
  done
}

# ---------------------------------------------------------------------------
# JSON extraction helper
# ---------------------------------------------------------------------------

# jq_get json_file python_expr — extract a scalar from JSON using only stdlib
jq_get() {
  python3 - "$1" "$2" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    d = json.load(f)
print(eval(sys.argv[2], {"d": d, "round": round}))
PY
}

# ---------------------------------------------------------------------------
# Extract a number from quaqc text output
# ---------------------------------------------------------------------------

# field_n "Fixed text pattern" file — return the first integer on the matching line
field_n() {
  grep -F "$1" "$2" | grep -oE '[0-9,]+' | head -1 | tr -d ','
}

# field_pct "Fixed text pattern" file — return the percent (without %) on matching line
field_pct() {
  grep -F "$1" "$2" | grep -oE '\([0-9.]+%\)' | head -1 | tr -d '()%'
}

# ---------------------------------------------------------------------------
# Create shared test BAMs once; each test uses its own copy or the shared one
# ---------------------------------------------------------------------------

SHARED="$TMPROOT/shared"
mkdir -p "$SHARED"

# basic_pe.bam — 10 PE pairs on chr1 (2000bp), MAPQ=42, 50+50bp reads, 150bp fragments
# Pairs are spread across chr1.
{
  chr_header chr1 2000
  for i in $(seq 1 10); do
    local_r1=$((100 + (i-1) * 150))
    local_r2=$((local_r1 + 100))
    pe_pair "rd${i}" chr1 "$local_r1" 50 "$local_r2" 50 42
  done
} | make_bam "$SHARED/basic_pe.bam"

# mapq_mixed.bam — 5 PE pairs MAPQ=42 + 3 PE pairs MAPQ=20 on chr1
{
  chr_header chr1 2000
  for i in $(seq 1 5); do
    r1=$((100 + (i-1)*150)); r2=$((r1+100))
    pe_pair "hi${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  for i in $(seq 1 3); do
    r1=$((900 + (i-1)*150)); r2=$((r1+100))
    pe_pair "lo${i}" chr1 "$r1" 50 "$r2" 50 20
  done
} | make_bam "$SHARED/mapq_mixed.bam"

# secondary.bam — 4 PE pairs MAPQ=42 + 2 secondary reads MAPQ=42, all on chr1
{
  chr_header chr1 2000
  for i in $(seq 1 4); do
    r1=$((100 + (i-1)*150)); r2=$((r1+100))
    pe_pair "pri${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  secondary_read sec1 chr1 1600 50 42
  secondary_read sec2 chr1 1700 50 42
} | make_bam "$SHARED/secondary.bam"

# dups.bam — 5 PE pairs MAPQ=42 (3 normal, 2 marked as PCR duplicates)
{
  chr_header chr1 2000
  for i in $(seq 1 3); do
    r1=$((100 + (i-1)*150)); r2=$((r1+100))
    pe_pair "nrm${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  dup_pair dup1 chr1 550 50 650 50 42
  dup_pair dup2 chr1 700 50 800 50 42
} | make_bam "$SHARED/dups.bam"

# multitype.bam — 4 PE pairs on chr1 (nuclear) + 2 PE pairs on chrM (mito)
{
  chr_header chr1 2000 chrM 1000
  for i in $(seq 1 4); do
    r1=$((100 + (i-1)*150)); r2=$((r1+100))
    pe_pair "nuc${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  pe_pair mito1 chrM 101 50 201 50 42
  pe_pair mito2 chrM 251 50 351 50 42
} | make_bam "$SHARED/multitype.bam"

# nested.bam — 2 PE pairs on chr1: pair1 large [100,300), pair2 nested [150,200)
# pair1 R1 at pos=101 len=200M (qend=300), pair2 R1 at pos=151 len=50M (qend=200).
# The inner read has a strictly shorter qend, making the coverage bug observable.
{
  chr_header chr1 1000
  pe_pair nest1 chr1 101 200 351 50 42
  pe_pair nest2 chr1 151  50 351 50 42
} | make_bam "$SHARED/nested.bam"

# plastid.bam — 3 PE pairs on chr1 (nuclear) + 2 PE pairs on chrC (plastid)
{
  chr_header chr1 2000 chrC 1000
  for i in $(seq 1 3); do
    r1=$((100 + (i-1)*150)); r2=$((r1+100))
    pe_pair "nuc${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  pe_pair plast1 chrC 101 50 201 50 42
  pe_pair plast2 chrC 251 50 351 50 42
} | make_bam "$SHARED/plastid.bam"

# mapq_bimodal.bam — 10 PE pairs MAPQ=30 + 10 PE pairs MAPQ=50 on chr1 (5000bp)
# Produces: avg=40.0, min=30, max=50, SD=10.0, 1st_pctile=30, 99th_pctile=50
{
  chr_header chr1 5000
  for i in $(seq 1 10); do
    r1=$((100 + (i-1)*200)); r2=$((r1+100))
    pe_pair "lo${i}" chr1 "$r1" 50 "$r2" 50 30
  done
  for i in $(seq 1 10); do
    r1=$((2200 + (i-1)*200)); r2=$((r1+100))
    pe_pair "hi${i}" chr1 "$r1" 50 "$r2" 50 50
  done
} | make_bam "$SHARED/mapq_bimodal.bam"

# flen_bimodal.bam — 5 PE pairs TLEN=100 + 5 PE pairs TLEN=300 on chr1 (2000bp)
# Produces: avg=200.0, min=100, max=300, SD=100.0, 1st_pctile=100, 99th_pctile=300
{
  chr_header chr1 2000
  for i in $(seq 1 5); do
    r1=$((100 + (i-1)*100)); r2=$((r1+50))
    pe_pair "sh${i}" chr1 "$r1" 50 "$r2" 50 42
  done
  for i in $(seq 1 5); do
    r1=$((700 + (i-1)*200)); r2=$((r1+250))
    pe_pair "lg${i}" chr1 "$r1" 50 "$r2" 50 42
  done
} | make_bam "$SHARED/flen_bimodal.bam"

# pileup3.bam — 3 PE pairs with R1 at identical position (depths_max=3)
{
  chr_header chr1 1000
  pe_pair pile1 chr1 101 50 301 50 42
  pe_pair pile2 chr1 101 50 301 50 42
  pe_pair pile3 chr1 101 50 301 50 42
} | make_bam "$SHARED/pileup3.bam"

# short_frags.bam — reads with small TLEN (fragment size 10) — below flen_min=15
# Using PE pairs where R1 and R2 overlap; TLEN = 10
{
  chr_header chr1 1000
  for i in $(seq 1 5); do
    r1=$((100 + (i-1)*20))
    r2=$((r1 + 5))
    printf 'short%d\t99\tchr1\t%d\t42\t10M\t=\t%d\t10\t%s\t*\n' "$i" "$r1" "$r2" "AAAAAAAAAA"
    printf 'short%d\t147\tchr1\t%d\t42\t5M\t=\t%d\t-10\t%s\t*\n' "$i" "$r2" "$r1" "TTTTT"
  done
} | make_bam "$SHARED/short_frags.bam"

# ---------------------------------------------------------------------------
# Individual test functions
# ---------------------------------------------------------------------------

# --- Regression: bug #1 (missing break in RG_TAG case) ---

test_rg_tag_no_use_secondary() {
  local d="$1"
  cp "$SHARED/secondary.bam" "$d/s.bam"
  cp "$SHARED/secondary.bam.bai" "$d/s.bam.bai"
  # Without --rg-tag: secondary excluded by default
  "$QUAQC" -o "$d" "$d/s.bam" >/dev/null
  # With --rg-tag XX: should NOT silently enable use-secondary (bug #1)
  "$QUAQC" -o "$d" --output-ext .rgtag.txt --rg-tag XX "$d/s.bam" >/dev/null
  assert_contains "$d/s.quaqc.txt"  "Include secondary alignments:              no" || return 1
  assert_contains "$d/s.rgtag.txt"  "Include secondary alignments:              no" || return 1
  # Reads passing filters must be identical
  local filt_base filt_rgtag
  filt_base=$(field_n   "Reads passing filters" "$d/s.quaqc.txt")
  filt_rgtag=$(field_n  "Reads passing filters" "$d/s.rgtag.txt")
  assert_eq "filt_n baseline vs rg-tag" "$filt_rgtag" "$filt_base"
}

# --- Regression: bug #2 (double-free destroy_depths) ---

test_no_crash_bed_bedgraph() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # Both --bedGraph and -e --bed-ins active: the bug would double-free bedGraph.
  # Should exit cleanly without memory corruption / crash.
  "$QUAQC" --no-out --bedGraph -e --bed-ins "$d/b.bam" >/dev/null
}

# --- Regression: bug #3 (sec_pct used se_total instead of sec_total) ---

test_secondary_pct_correct() {
  local d="$1"
  cp "$SHARED/secondary.bam" "$d/s.bam"
  cp "$SHARED/secondary.bam.bai" "$d/s.bam.bai"
  "$QUAQC" -o "$d" "$d/s.bam" >/dev/null
  local report="$d/s.quaqc.txt"
  # secondary.bam: 4 PE pairs = 8 primary reads + 2 secondary reads = 10 r_seen
  # sec_pct = 2/10 * 100 = 20.0 %
  local sec_pct; sec_pct=$(field_pct "Secondary alignments:" "$report")
  assert_eq "secondary pct" "$sec_pct" "20.0"
}

# --- Regression: bug #4 (JSON depth histogram range used fhist_max not depth_max) ---

test_json_depth_range() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  # Default depth_histogram_max = 100000; frag_histogram_max = 2000
  # Before fix, range[1] would be 2000 (wrong). After fix it is 100000.
  local range1
  range1=$(jq_get "$d/out.json" \
    "d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['depth']['depths_histogram']['range'][1]")
  assert_eq "depth histogram range[1]" "$range1" "100000"
}

# --- Regression: bug #5 (warn_flen !!! triple-negation) ---

test_warn_flen_populated() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # basic_pe.bam has fragments of 150bp; set --max-fhist 100 so bin[100] is non-zero
  "$QUAQC" --no-out --json "$d/out.json" --max-fhist 100 "$d/b.bam" >/dev/null
  local suspect
  suspect=$(jq_get "$d/out.json" \
    "d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['fragment']['addn_stats_are_suspect']")
  # Fragment sizes (150) exceed histogram max (100), so last bin is populated → true
  assert_eq "addn_stats_are_suspect" "$suspect" "True"
}

# --- Regression: bug #6 (coverage undercount for nested reads) ---
# Two nuclear PE pairs where the inner R1 has a strictly shorter qend than the outer R1.
# Outer: R1=[100,300), inner: R1=[150,200). Both R2s at [350,400).
# Fixed non-overlapping cov = 200+50 = 250; bug shrinks last_end to 200 → cov = 100+50 = 150.
# genome_cov = non-overlapping-bases / chrom-size = 250/1000 = 0.25 (fixed), 0.15 (bug).

test_nested_read_coverage() {
  local d="$1"
  {
    chr_header chr1 1000
    # Outer: R1 pos=101 200M (qend=300), R2 pos=351 50M (qend=400); TLEN=300
    pe_pair nest1 chr1 101 200 351 50 42
    # Inner: R1 pos=151 50M (qend=200, shorter than outer's qend=300!), same R2; TLEN=250
    pe_pair nest2 chr1 151  50 351 50 42
  } | make_bam "$d/nested.bam"
  "$QUAQC" --no-out --json "$d/out.json" "$d/nested.bam" >/dev/null
  local genome_cov
  genome_cov=$(jq_get "$d/out.json" \
    "round(d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['depth']['genome_cov'] * 100)")
  # Fixed: genome_cov = 250/1000 = 0.25 → round(0.25*100) = 25
  # Bug:   genome_cov = 150/1000 = 0.15 → round(0.15*100) = 15
  assert_eq "genome_cov (nested reads)" "$genome_cov" "25"
}

# --- Regression: bug #8 (gzerror called on NULL handle) ---

test_gzopen_fail_graceful() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # Write bedGraph to a nonexistent directory — should print an error, not crash/segfault
  "$QUAQC" --no-out --bedGraph --bedGraph-dir "$d/nonexistent_dir" "$d/b.bam" >/dev/null 2>&1 || true
  # The test passes as long as quaqc doesn't segfault (exit code != 139)
  local rc=$?
  [ "$rc" -ne 139 ]
}

# --- Regression: bugs #13/#14 (bed-ins coordinate and off-by-one) ---

test_bed_ins_coordinates() {
  local d="$1"
  # One PE pair, forward R1 at pos=101 (1-based, 0-based pos=100), 50M
  # With --bed-tn5 (tn5_fwd=4): insertion at pos=100+4=104 (0-based)
  # BED line should be: chr1 \t 104 \t 105
  {
    chr_header chr1 1000
    pe_pair ins1 chr1 101 50 201 50 42
  } | make_bam "$d/ins.bam"
  # -e enables BED output; --bed-ins switches to BED3 insertion format.
  # BED file goes next to BAM by default (in $d).
  "$QUAQC" --no-out -e --bed-ins --bed-tn5 "$d/ins.bam" >/dev/null
  local bedfile="$d/ins.bed.gz"
  assert_file_exists "$bedfile"
  local entries; entries=$(gzip -dc "$bedfile" | wc -l | tr -d ' ')
  [ "$entries" -ge 1 ] || { echo "  no BED entries written" >&2; return 1; }
  # Forward strand insertion should be at 104 (0-based start)
  if gzip -dc "$bedfile" | grep -q $'^chr1\t104\t105$'; then
    return 0
  fi
  echo "  expected chr1:104-105 in BED; got: $(gzip -dc "$bedfile")" >&2
  return 1
}

# ---------------------------------------------------------------------------
# Feature correctness tests
# ---------------------------------------------------------------------------

test_basic_read_counts() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" -o "$d" "$d/b.bam" >/dev/null
  local report="$d/b.quaqc.txt"
  # 10 PE pairs = 20 mapped reads
  assert_eq "total reads"   "$(field_n 'Total reads:'   "$report")" "20" || return 1
  assert_eq "mapped reads"  "$(field_n '├── Mapped:'     "$report")" "20" || return 1
  # 10 PE pairs all pass filters (MAPQ=42, qlen=50, flen=150)
  assert_eq "passing reads" "$(field_n 'Reads passing filters:' "$report")" "20" || return 1
  # 10 fragments (one per PE pair)
  assert_eq "fragments"     "$(field_n 'Fragments passing filters:' "$report")" "10"
}

test_mapq_filter() {
  local d="$1"
  cp "$SHARED/mapq_mixed.bam" "$d/m.bam"
  cp "$SHARED/mapq_mixed.bam.bai" "$d/m.bam.bai"
  "$QUAQC" -o "$d" "$d/m.bam" >/dev/null
  # 5 high-MAPQ pairs (10 reads) pass; 3 low-MAPQ pairs (6 reads) filtered
  assert_eq "reads passing MAPQ filter" \
    "$(field_n 'Reads passing filters:' "$d/m.quaqc.txt")" "10"
}

test_mapq_threshold_change() {
  local d="$1"
  cp "$SHARED/mapq_mixed.bam" "$d/m.bam"
  cp "$SHARED/mapq_mixed.bam.bai" "$d/m.bam.bai"
  # Lowering MAPQ threshold to 10 should include the MAPQ=20 reads too
  "$QUAQC" -o "$d" --mapq 10 "$d/m.bam" >/dev/null
  assert_eq "reads with mapq>=10" \
    "$(field_n 'Reads passing filters:' "$d/m.quaqc.txt")" "16"
}

test_flen_filter() {
  local d="$1"
  cp "$SHARED/short_frags.bam" "$d/sf.bam"
  cp "$SHARED/short_frags.bam.bai" "$d/sf.bam.bai"
  "$QUAQC" -o "$d" "$d/sf.bam" >/dev/null
  # Fragment size = 10 < min_flen = 15 → all filtered
  assert_eq "short frags filtered" \
    "$(field_n 'Reads passing filters:' "$d/sf.quaqc.txt")" "0"
}

test_flen_min_override() {
  local d="$1"
  cp "$SHARED/short_frags.bam" "$d/sf.bam"
  cp "$SHARED/short_frags.bam.bai" "$d/sf.bam.bai"
  # Short reads have qlen=10; need --min-qlen 5 in addition to --min-flen 5
  "$QUAQC" -o "$d" --min-flen 5 --min-qlen 5 "$d/sf.bam" >/dev/null
  assert_eq "short frags with min-flen 5" \
    "$(field_n 'Reads passing filters:' "$d/sf.quaqc.txt")" "10"
}

test_mito_classification() {
  local d="$1"
  cp "$SHARED/multitype.bam" "$d/mt.bam"
  cp "$SHARED/multitype.bam.bai" "$d/mt.bam.bai"
  "$QUAQC" -o "$d" "$d/mt.bam" >/dev/null
  local report="$d/mt.quaqc.txt"
  # 4 nuclear pairs = 8 nuclear reads; 2 mito pairs = 4 mito reads
  assert_eq "nuclear reads" \
    "$(grep '├┬─ Nuclear:'  "$report" | grep -oE '[0-9,]+' | head -1 | tr -d ',')" "8" || return 1
  assert_eq "mito reads" \
    "$(grep '├┬─ Mitochondrial:' "$report" | grep -oE '[0-9,]+' | head -1 | tr -d ',')" "4"
}

test_use_secondary() {
  local d="$1"
  cp "$SHARED/secondary.bam" "$d/s.bam"
  cp "$SHARED/secondary.bam.bai" "$d/s.bam.bai"
  # Without --use-secondary: only 8 primary reads pass
  "$QUAQC" -o "$d" --output-ext .nosec.txt "$d/s.bam" >/dev/null
  # With --use-secondary: secondary reads now pass MAPQ/qlen filters
  # but still fail flen filter (isize=0 < 15). Use --min-flen 0 to allow them.
  "$QUAQC" -o "$d" --output-ext .sec.txt --use-secondary --min-flen 0 "$d/s.bam" >/dev/null
  local without with
  without=$(field_n 'Reads passing filters:' "$d/s.nosec.txt")
  with=$(field_n    'Reads passing filters:' "$d/s.sec.txt")
  # With --use-secondary + --min-flen 0, the 2 secondary reads are also included
  [ "$with" -gt "$without" ] || {
    echo "  --use-secondary had no effect (without=$without, with=$with)" >&2
    return 1
  }
}

test_use_dups() {
  local d="$1"
  cp "$SHARED/dups.bam" "$d/d.bam"
  cp "$SHARED/dups.bam.bai" "$d/d.bam.bai"
  # Without --use-dups: only 3 non-dup pairs = 6 reads pass
  "$QUAQC" -o "$d" --output-ext .nodups.txt "$d/d.bam" >/dev/null
  # With --use-dups: 5 pairs = 10 reads pass
  "$QUAQC" -o "$d" --output-ext .dups.txt --use-dups "$d/d.bam" >/dev/null
  local without with
  without=$(field_n 'Reads passing filters:' "$d/d.nodups.txt")
  with=$(field_n    'Reads passing filters:' "$d/d.dups.txt")
  assert_eq "no-dup baseline"   "$without" "6" || return 1
  assert_eq "with --use-dups"   "$with"    "10"
}

test_pe_fragment_count() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" -o "$d" "$d/b.bam" >/dev/null
  # 10 PE pairs → 10 R1 reads → 10 fragments
  local frags; frags=$(field_n 'Fragments passing filters:' "$d/b.quaqc.txt")
  assert_eq "fragment count" "$frags" "10"
}

test_fragment_avg_size() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" -o "$d" "$d/b.bam" >/dev/null
  # All fragments are 150bp (r1=50M, gap=50, r2=50M, TLEN=150)
  local avg; avg=$(grep 'Average:' "$d/b.quaqc.txt" | head -2 | tail -1 | grep -oE '[0-9]+\.[0-9]+')
  # Fragment average should be 150.000
  [[ "$avg" == "150.000" ]] || {
    echo "  fragment avg: expected 150.000, got '$avg'" >&2
    return 1
  }
}

# ---------------------------------------------------------------------------
# Output mode tests
# ---------------------------------------------------------------------------

test_json_valid() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  python3 -c "import json; json.load(open('$d/out.json'))" || {
    echo "  JSON failed to parse" >&2; return 1
  }
}

test_no_output_flag() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" --no-out -o "$d" "$d/b.bam" >/dev/null
  # Report file must NOT be created
  [ ! -f "$d/b.quaqc.txt" ] || {
    echo "  report file unexpectedly created" >&2; return 1
  }
}

test_output_ext() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" -o "$d" --output-ext .myreport.txt "$d/b.bam" >/dev/null
  assert_file_exists "$d/b.myreport.txt"
}

test_save_bam_passthrough() {
  local d="$1"
  cp "$SHARED/mapq_mixed.bam" "$d/m.bam"
  cp "$SHARED/mapq_mixed.bam.bai" "$d/m.bam.bai"
  "$QUAQC" --no-out --keep -o "$d" "$d/m.bam" >/dev/null
  local saved="$d/m.filt.bam"
  assert_file_exists "$saved" || return 1
  # Saved BAM must contain exactly the passing reads (10, not 16)
  local saved_count; saved_count=$(samtools view -c "$saved" 2>/dev/null)
  assert_eq "saved BAM read count" "$saved_count" "10"
}

test_bedgraph_output() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  "$QUAQC" --no-out --bedGraph --bedGraph-dir "$d" "$d/b.bam" >/dev/null
  assert_file_exists "$d/b.bedGraph.gz"
  # Must be readable gzip
  local lines; lines=$(gzip -dc "$d/b.bedGraph.gz" | wc -l | tr -d ' ')
  [ "$lines" -gt 0 ] || { echo "  empty bedGraph" >&2; return 1; }
}

test_multiple_bams() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"   "$d/a.bam";  cp "$SHARED/basic_pe.bam.bai"   "$d/a.bam.bai"
  cp "$SHARED/mapq_mixed.bam" "$d/b.bam";  cp "$SHARED/mapq_mixed.bam.bai" "$d/b.bam.bai"
  "$QUAQC" -o "$d" "$d/a.bam" "$d/b.bam" >/dev/null
  assert_file_exists "$d/a.quaqc.txt" || return 1
  assert_file_exists "$d/b.quaqc.txt"
}

test_multithreaded_same_output() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/a.bam"; cp "$SHARED/basic_pe.bam.bai" "$d/a.bam.bai"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"; cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # Single-threaded
  "$QUAQC" -o "$d" --output-ext .s1.txt "$d/a.bam" "$d/b.bam" >/dev/null
  # Two-threaded (one thread per file)
  "$QUAQC" -o "$d" --output-ext .t2.txt -j 2 "$d/a.bam" "$d/b.bam" >/dev/null
  # Reads passing filters must be identical in both runs
  local s1 t2
  s1=$(field_n 'Reads passing filters:' "$d/a.s1.txt")
  t2=$(field_n 'Reads passing filters:' "$d/a.t2.txt")
  assert_eq "single vs multi-thread" "$s1" "$t2"
}

test_fast_mode() {
  local d="$1"
  cp "$SHARED/basic_pe.bam" "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # --fast limits to first chromosome; still produces valid output
  "$QUAQC" -o "$d" --fast --output-ext .fast.txt "$d/b.bam" >/dev/null
  assert_file_exists "$d/b.fast.txt"
}

# ---------------------------------------------------------------------------
# Additional feature tests
# ---------------------------------------------------------------------------

test_blacklist_filtering() {
  local d="$1"
  # 4 PE pairs on chr1. Pair 1 R1=[100,150) R2=[200,250) fall inside blacklist [0,300).
  # Pairs 2-4 start at positions 350+ and are unaffected.
  {
    chr_header chr1 2000
    pe_pair bl1 chr1 101 50 201 50 42
    pe_pair ok1 chr1 351 50 451 50 42
    pe_pair ok2 chr1 551 50 651 50 42
    pe_pair ok3 chr1 751 50 851 50 42
  } | make_bam "$d/bl.bam"
  printf 'chr1\t0\t300\n' > "$d/blacklist.bed"
  "$QUAQC" -o "$d" --output-ext .nbl.txt "$d/bl.bam" >/dev/null
  "$QUAQC" -o "$d" --output-ext .bl.txt --blacklist "$d/blacklist.bed" "$d/bl.bam" >/dev/null
  local without with
  without=$(field_n 'Reads passing filters:' "$d/bl.nbl.txt")
  with=$(field_n 'Reads passing filters:' "$d/bl.bl.txt")
  assert_eq "without blacklist" "$without" "8" || return 1
  assert_eq "with blacklist"    "$with"    "6"
}

test_target_list_filtering() {
  local d="$1"
  # 4 PE pairs: pairs 1-2 are fully within target [50,600); pairs 3-4 are outside.
  # bed_overlap_within requires the read to be entirely enclosed in a target region.
  {
    chr_header chr1 2000
    pe_pair in1 chr1 101 50 251 50 42   # R1=[100,150), R2=[250,300): within [50,600)
    pe_pair in2 chr1 351 50 451 50 42   # R1=[350,400), R2=[450,500): within [50,600)
    pe_pair out1 chr1 701 50 801 50 42  # R1=[700,750): outside [50,600)
    pe_pair out2 chr1 901 50 1001 50 42 # R1=[900,950): outside [50,600)
  } | make_bam "$d/tl.bam"
  printf 'chr1\t50\t600\n' > "$d/target.bed"
  "$QUAQC" -o "$d" --output-ext .notl.txt "$d/tl.bam" >/dev/null
  "$QUAQC" -o "$d" --output-ext .tl.txt --target-list "$d/target.bed" "$d/tl.bam" >/dev/null
  local without with
  without=$(field_n 'Reads passing filters:' "$d/tl.notl.txt")
  with=$(field_n 'Reads passing filters:' "$d/tl.tl.txt")
  assert_eq "without target" "$without" "8" || return 1
  assert_eq "with target"    "$with"    "4"
}

test_plastid_classification() {
  local d="$1"
  cp "$SHARED/plastid.bam"     "$d/p.bam"
  cp "$SHARED/plastid.bam.bai" "$d/p.bam.bai"
  "$QUAQC" -o "$d" "$d/p.bam" >/dev/null
  local report="$d/p.quaqc.txt"
  # 3 nuclear pairs = 6 reads; 2 plastid pairs = 4 reads
  assert_eq "nuclear reads" \
    "$(grep '├┬─ Nuclear:'  "$report" | grep -oE '[0-9,]+' | head -1 | tr -d ',')" "6" || return 1
  assert_eq "plastid reads" \
    "$(grep '├┬─ Plastid:' "$report" | grep -oE '[0-9,]+' | head -1 | tr -d ',')" "4"
}

test_max_qlen_filter() {
  local d="$1"
  # 5 SE reads with qlen=50 (pass max_qlen=250) + 3 SE reads with qlen=300 (fail).
  {
    chr_header chr1 2000
    for i in $(seq 1 5); do se_read "s${i}" chr1 $((100+(i-1)*60)) 50 42; done
    for i in $(seq 1 3); do se_read "l${i}" chr1 $((500+(i-1)*400)) 300 42; done
  } | make_bam "$d/q.bam"
  "$QUAQC" -o "$d" --output-ext .dflt.txt "$d/q.bam" >/dev/null
  "$QUAQC" -o "$d" --output-ext .mx.txt --max-qlen 350 "$d/q.bam" >/dev/null
  assert_eq "default max-qlen (5 pass)" \
    "$(field_n 'Reads passing filters:' "$d/q.dflt.txt")" "5" || return 1
  assert_eq "with --max-qlen 350 (8 pass)" \
    "$(field_n 'Reads passing filters:' "$d/q.mx.txt")" "8"
}

test_max_flen_filter() {
  local d="$1"
  # 3 pairs TLEN=150 (pass default max_flen=2000) + 2 pairs TLEN=2100 (fail).
  {
    chr_header chr1 3000
    pe_pair p1 chr1 101 50 201 50 42   # TLEN=150
    pe_pair p2 chr1 301 50 401 50 42   # TLEN=150
    pe_pair p3 chr1 501 50 601 50 42   # TLEN=150
    pe_pair p4 chr1 701 50 2751 50 42  # TLEN=701+50+2751-701=2100 → 2751+50-701=2100
    pe_pair p5 chr1 801 50 2851 50 42  # TLEN=2100
  } | make_bam "$d/f.bam"
  "$QUAQC" -o "$d" --output-ext .dflt.txt "$d/f.bam" >/dev/null
  "$QUAQC" -o "$d" --output-ext .mx.txt --max-flen 3000 "$d/f.bam" >/dev/null
  assert_eq "default max-flen (6 pass)" \
    "$(field_n 'Reads passing filters:' "$d/f.dflt.txt")" "6" || return 1
  assert_eq "with --max-flen 3000 (10 pass)" \
    "$(field_n 'Reads passing filters:' "$d/f.mx.txt")" "10"
}

test_frip_calculation() {
  local d="$1"
  # 3 PE pairs; pairs 1-2 have all reads within peak [0,500); pair 3 is outside.
  {
    chr_header chr1 2000
    pe_pair p1 chr1 101 50 251 50 42   # R1=[100,150), R2=[250,300): within peak
    pe_pair p2 chr1 351 50 401 50 42   # R1=[350,400), R2=[400,450): within peak
    pe_pair p3 chr1 601 50 701 50 42   # R1=[600,650), R2=[700,750): outside peak
  } | make_bam "$d/frip.bam"
  printf 'chr1\t0\t500\n' > "$d/peaks.bed"
  "$QUAQC" --no-out --json "$d/out.json" --peaks "$d/peaks.bed" "$d/frip.bam" >/dev/null
  # 4 reads overlap the peak; 2 do not → FRIP = 4/6
  local n_in_peaks
  n_in_peaks=$(jq_get "$d/out.json" \
    "round(d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['peaks']['fraction_of_reads_in_peaks'] * 6)")
  assert_eq "reads in peaks" "$n_in_peaks" "4"
}

test_quant_output() {
  local d="$1"
  # 2 PE pairs fully within peak; 1 pair outside. --quant produces a TSV count matrix.
  {
    chr_header chr1 2000
    pe_pair p1 chr1 101 50 201 50 42   # within peak [0,500)
    pe_pair p2 chr1 301 50 401 50 42   # within peak [0,500)
    pe_pair p3 chr1 601 50 701 50 42   # outside peak
  } | make_bam "$d/q.bam"
  printf 'chr1\t0\t500\n' > "$d/peaks.bed"
  "$QUAQC" -0 --peaks "$d/peaks.bed" --quant "$d/counts.tsv" "$d/q.bam" >/dev/null
  assert_file_exists "$d/counts.tsv" || return 1
  # 2 rows: header + 1 peak
  local lines; lines=$(wc -l < "$d/counts.tsv" | tr -d ' ')
  assert_eq "quant row count" "$lines" "2" || return 1
  # Count for the peak must be > 0
  local count; count=$(tail -1 "$d/counts.tsv" | cut -f2)
  [ "$count" -gt 0 ] || { echo "  quant count is 0" >&2; return 1; }
}

test_json_multiple_reports() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/a.bam"; cp "$SHARED/basic_pe.bam.bai"     "$d/a.bam.bai"
  cp "$SHARED/mapq_mixed.bam"   "$d/b.bam"; cp "$SHARED/mapq_mixed.bam.bai"   "$d/b.bam.bai"
  "$QUAQC" --no-out --json "$d/out.json" "$d/a.bam" "$d/b.bam" >/dev/null
  local nreports
  nreports=$(jq_get "$d/out.json" "len(d['quaqc_reports'])")
  assert_eq "report count" "$nreports" "2"
}

test_strict_preset() {
  local d="$1"
  # 3 pairs MAPQ=42 TLEN=150 — pass both default and strict filters.
  # 2 pairs MAPQ=35 TLEN=150 — pass default (mapq≥30) but fail strict (mapq≥40).
  # 2 pairs MAPQ=42 TLEN=200 — pass default (max_flen=2000) but fail strict (max_flen=150).
  {
    chr_header chr1 2000
    pe_pair p1 chr1  101 50  201 50 42   # TLEN=150, MAPQ=42
    pe_pair p2 chr1  301 50  401 50 42   # TLEN=150, MAPQ=42
    pe_pair p3 chr1  501 50  601 50 42   # TLEN=150, MAPQ=42
    pe_pair p4 chr1  701 50  801 50 35   # TLEN=150, MAPQ=35
    pe_pair p5 chr1  901 50 1001 50 35   # TLEN=150, MAPQ=35
    pe_pair p6 chr1 1101 50 1251 50 42   # TLEN=200, MAPQ=42
    pe_pair p7 chr1 1301 50 1451 50 42   # TLEN=200, MAPQ=42
  } | make_bam "$d/s.bam"
  "$QUAQC" -o "$d" --output-ext .dflt.txt "$d/s.bam" >/dev/null
  "$QUAQC" -o "$d" --output-ext .strict.txt --strict "$d/s.bam" >/dev/null
  assert_eq "default (14 reads)" \
    "$(field_n 'Reads passing filters:' "$d/s.dflt.txt")" "14" || return 1
  assert_eq "strict (6 reads)" \
    "$(field_n 'Reads passing filters:' "$d/s.strict.txt")" "6"
}

test_bed_output() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # -e writes one BED6 line per passing read; 10 PE pairs = 20 reads
  "$QUAQC" --no-out -e --bed-dir "$d" "$d/b.bam" >/dev/null
  assert_file_exists "$d/b.bed.gz" || return 1
  local lines; lines=$(gzip -dc "$d/b.bed.gz" | wc -l | tr -d ' ')
  assert_eq "BED line count" "$lines" "20"
}

# ---------------------------------------------------------------------------
# Output statistics correctness tests
# ---------------------------------------------------------------------------

test_qlen_stats_uniform() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # All 20 reads are 50bp → avg=50.0, min=50, max=50, SD=0.0
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  local nuc="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['alignment']"
  assert_eq "avg qlen" "$(jq_get "$d/out.json" "${nuc}['size_average']")" "50.0" || return 1
  assert_eq "min qlen" "$(jq_get "$d/out.json" "${nuc}['size_min']")"     "50"   || return 1
  assert_eq "max qlen" "$(jq_get "$d/out.json" "${nuc}['size_max']")"     "50"   || return 1
  assert_eq "sd qlen"  "$(jq_get "$d/out.json" "${nuc}['size_sd']")"      "0.0"
}

test_flen_stats_uniform() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # All 10 fragments are TLEN=150 → avg=150.0, min=150, max=150, SD=0.0
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  local frag="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['fragment']"
  assert_eq "avg flen" "$(jq_get "$d/out.json" "${frag}['size_average']")" "150.0" || return 1
  assert_eq "min flen" "$(jq_get "$d/out.json" "${frag}['size_min']")"     "150"   || return 1
  assert_eq "max flen" "$(jq_get "$d/out.json" "${frag}['size_max']")"     "150"   || return 1
  assert_eq "sd flen"  "$(jq_get "$d/out.json" "${frag}['size_sd']")"      "0.0"
}

test_mapq_avg_mixed() {
  local d="$1"
  cp "$SHARED/mapq_mixed.bam"     "$d/m.bam"
  cp "$SHARED/mapq_mixed.bam.bai" "$d/m.bam.bai"
  # --mapq 10 lets all 16 reads pass: 10×MAPQ42 + 6×MAPQ20 → avg = 540/16 = 33.75
  "$QUAQC" --no-out --json "$d/out.json" --mapq 10 "$d/m.bam" >/dev/null
  local avg
  avg=$(jq_get "$d/out.json" \
    "d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['mapq']['score_average']")
  assert_eq "avg mapq" "$avg" "33.75"
}

test_read_depth_calculation() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # 20 reads × 50bp = 1000 non-overlapping bases on 2000bp chr1.
  # read_depth_average = 1000/2000 = 0.5; genome_cov = 1000/2000 = 0.5.
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  local depth="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['depth']"
  assert_eq "read_depth_average" \
    "$(jq_get "$d/out.json" "${depth}['read_depth_average']")" "0.5" || return 1
  assert_eq "genome_cov" \
    "$(jq_get "$d/out.json" "${depth}['genome_cov']")"         "0.5"
}

test_gc_content_known() {
  local d="$1"
  # 4 SE reads, each 50bp = 25'A' + 25'C': exactly 50% GC per read.
  local seq; seq=$(python3 -c "print('A'*25+'C'*25, end='')")
  {
    chr_header chr1 1000
    for i in $(seq 1 4); do
      printf 'gc%d\t0\tchr1\t%d\t42\t50M\t*\t0\t0\t%s\t*\n' \
        "$i" $((100 + (i-1)*100)) "$seq"
    done
  } | make_bam "$d/gc.bam"
  "$QUAQC" --no-out --json "$d/out.json" "$d/gc.bam" >/dev/null
  local gc="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['gc']"
  assert_eq "gc pct_average" "$(jq_get "$d/out.json" "${gc}['pct_average']")" "50.0" || return 1
  assert_eq "gc pct_min"     "$(jq_get "$d/out.json" "${gc}['pct_min']")"     "50"   || return 1
  assert_eq "gc pct_max"     "$(jq_get "$d/out.json" "${gc}['pct_max']")"     "50"
}

test_nrf_calculation() {
  local d="$1"
  cp "$SHARED/dups.bam"     "$d/d.bam"
  cp "$SHARED/dups.bam.bai" "$d/d.bam.bai"
  # 3 non-dup pairs (6 primary reads) + 2 dup pairs (4 dup-flagged reads) = 10 total.
  # NRF = (10 - 4) / 10 = 0.6
  "$QUAQC" --no-out --json "$d/out.json" "$d/d.bam" >/dev/null
  local nrf
  nrf=$(jq_get "$d/out.json" \
    "d['quaqc_reports'][0]['report']['unfiltered_read_stats']['effective']['nuclear']['non_redundant_fraction']")
  assert_eq "NRF" "$nrf" "0.6"
}

test_qlen_mixed_stats() {
  local d="$1"
  # 5 PE pairs with 50bp reads + 5 PE pairs with 100bp reads.
  # Total: 10 reads of 50bp + 10 reads of 100bp = 20 reads.
  # avg = (10×50 + 10×100)/20 = 75.0; min=50; max=100.
  {
    chr_header chr1 2000
    for i in $(seq 1 5); do
      pe_pair "s${i}" chr1 $((100+(i-1)*100)) 50 $((150+(i-1)*100)) 50 42
    done
    for i in $(seq 1 5); do
      pe_pair "l${i}" chr1 $((700+(i-1)*200)) 100 $((850+(i-1)*200)) 100 42
    done
  } | make_bam "$d/mix.bam"
  "$QUAQC" --no-out --json "$d/out.json" "$d/mix.bam" >/dev/null
  local nuc="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['alignment']"
  assert_eq "avg qlen (mixed)" "$(jq_get "$d/out.json" "${nuc}['size_average']")" "75.0" || return 1
  assert_eq "min qlen (mixed)" "$(jq_get "$d/out.json" "${nuc}['size_min']")"     "50"   || return 1
  assert_eq "max qlen (mixed)" "$(jq_get "$d/out.json" "${nuc}['size_max']")"     "100"
}

test_mapq_stats_bimodal() {
  local d="$1"
  cp "$SHARED/mapq_bimodal.bam"     "$d/m.bam"
  cp "$SHARED/mapq_bimodal.bam.bai" "$d/m.bam.bai"
  # 20×MAPQ30 + 20×MAPQ50 reads: avg=40.0, min=30, max=50, SD=10.0
  "$QUAQC" --no-out --json "$d/out.json" "$d/m.bam" >/dev/null
  local mapq="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['mapq']"
  assert_eq "mapq avg" "$(jq_get "$d/out.json" "${mapq}['score_average']")" "40.0"  || return 1
  assert_eq "mapq min" "$(jq_get "$d/out.json" "${mapq}['score_min']")"     "30"    || return 1
  assert_eq "mapq max" "$(jq_get "$d/out.json" "${mapq}['score_max']")"     "50"    || return 1
  assert_eq "mapq sd"  "$(jq_get "$d/out.json" "${mapq}['score_sd']")"      "10.0"
}

test_mapq_percentiles() {
  local d="$1"
  cp "$SHARED/mapq_bimodal.bam"     "$d/m.bam"
  cp "$SHARED/mapq_bimodal.bam.bai" "$d/m.bam.bai"
  # 20×MAPQ30 + 20×MAPQ50: mass-weighted 1st pctile=30, 99th pctile=50
  "$QUAQC" --no-out --json "$d/out.json" "$d/m.bam" >/dev/null
  local mapq="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['mapq']"
  assert_eq "mapq 1st pctile"  "$(jq_get "$d/out.json" "${mapq}['score_1st_pctile']")"  "30" || return 1
  assert_eq "mapq 99th pctile" "$(jq_get "$d/out.json" "${mapq}['score_99th_pctile']")" "50"
}

test_flen_mixed_stats() {
  local d="$1"
  cp "$SHARED/flen_bimodal.bam"     "$d/f.bam"
  cp "$SHARED/flen_bimodal.bam.bai" "$d/f.bam.bai"
  # 5×TLEN100 + 5×TLEN300: avg=200.0, min=100, max=300, SD=100.0
  "$QUAQC" --no-out --json "$d/out.json" "$d/f.bam" >/dev/null
  local frag="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['fragment']"
  assert_eq "flen avg" "$(jq_get "$d/out.json" "${frag}['size_average']")" "200.0" || return 1
  assert_eq "flen min" "$(jq_get "$d/out.json" "${frag}['size_min']")"     "100"   || return 1
  assert_eq "flen max" "$(jq_get "$d/out.json" "${frag}['size_max']")"     "300"   || return 1
  assert_eq "flen sd"  "$(jq_get "$d/out.json" "${frag}['size_sd']")"      "100.0"
}

test_flen_percentiles() {
  local d="$1"
  cp "$SHARED/flen_bimodal.bam"     "$d/f.bam"
  cp "$SHARED/flen_bimodal.bam.bai" "$d/f.bam.bai"
  # 5×TLEN100 + 5×TLEN300: mass-weighted 1st pctile=100, 99th pctile=300
  "$QUAQC" --no-out --json "$d/out.json" "$d/f.bam" >/dev/null
  local frag="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['fragment']"
  assert_eq "flen 1st pctile"  "$(jq_get "$d/out.json" "${frag}['size_1st_pctile']")"  "100" || return 1
  assert_eq "flen 99th pctile" "$(jq_get "$d/out.json" "${frag}['size_99th_pctile']")" "300"
}

test_depths_stats_basic() {
  local d="$1"
  cp "$SHARED/basic_pe.bam"     "$d/b.bam"
  cp "$SHARED/basic_pe.bam.bai" "$d/b.bam.bai"
  # 1000 covered bases + 1000 uncovered out of 2000bp: min=0, max=1, SD=0.5
  "$QUAQC" --no-out --json "$d/out.json" "$d/b.bam" >/dev/null
  local dep="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['depth']"
  assert_eq "depths_min" "$(jq_get "$d/out.json" "${dep}['depths_min']")" "0"   || return 1
  assert_eq "depths_max" "$(jq_get "$d/out.json" "${dep}['depths_max']")" "1"   || return 1
  assert_eq "depths_sd"  "$(jq_get "$d/out.json" "${dep}['depths_sd']")"  "0.5"
}

test_depths_max_pileup() {
  local d="$1"
  cp "$SHARED/pileup3.bam"     "$d/p.bam"
  cp "$SHARED/pileup3.bam.bai" "$d/p.bam.bai"
  # 3 PE pairs with R1 at identical position → max read depth of 3
  "$QUAQC" --no-out --json "$d/out.json" "$d/p.bam" >/dev/null
  local dep="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['depth']"
  assert_eq "depths_max" "$(jq_get "$d/out.json" "${dep}['depths_max']")" "3"
}

test_gc_mixed_extremes() {
  local d="$1"
  local seq_a; seq_a=$(python3 -c "print('A'*50, end='')")
  local seq_c; seq_c=$(python3 -c "print('C'*50, end='')")
  {
    chr_header chr1 1000
    for i in $(seq 1 4); do
      printf 'ga%d\t0\tchr1\t%d\t42\t50M\t*\t0\t0\t%s\t*\n' "$i" $((100+(i-1)*100)) "$seq_a"
    done
    for i in $(seq 1 4); do
      printf 'gc%d\t0\tchr1\t%d\t42\t50M\t*\t0\t0\t%s\t*\n' "$i" $((550+(i-1)*100)) "$seq_c"
    done
  } | make_bam "$d/gc.bam"
  # 4 all-A (0% GC) + 4 all-C (100% GC): avg=50.0, min=0, max=100, SD=50.0
  "$QUAQC" --no-out --json "$d/out.json" "$d/gc.bam" >/dev/null
  local gc="d['quaqc_reports'][0]['report']['filtered_read_stats']['nuclear']['gc']"
  assert_eq "gc avg" "$(jq_get "$d/out.json" "${gc}['pct_average']")" "50.0" || return 1
  assert_eq "gc min" "$(jq_get "$d/out.json" "${gc}['pct_min']")"     "0"    || return 1
  assert_eq "gc max" "$(jq_get "$d/out.json" "${gc}['pct_max']")"     "100"  || return 1
  assert_eq "gc sd"  "$(jq_get "$d/out.json" "${gc}['pct_sd']")"      "50.0"
}

test_unfiltered_secondary_counts() {
  local d="$1"
  cp "$SHARED/secondary.bam"     "$d/s.bam"
  cp "$SHARED/secondary.bam.bai" "$d/s.bam.bai"
  # 4 PE pairs (8 primary reads) + 2 secondary SE reads = 10 total nuclear reads
  "$QUAQC" --no-out --json "$d/out.json" "$d/s.bam" >/dev/null
  local nuc="d['quaqc_reports'][0]['report']['unfiltered_read_stats']['effective']['nuclear']"
  assert_eq "nuclear n"           "$(jq_get "$d/out.json" "${nuc}['n']")"                    "10" || return 1
  assert_eq "primary_alignments"  "$(jq_get "$d/out.json" "${nuc}['primary_alignments']")"   "8"  || return 1
  assert_eq "secondary_alignments" "$(jq_get "$d/out.json" "${nuc}['secondary_alignments']")" "2"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

echo
echo "quaqc test suite"
echo "================"
echo

for t in \
  rg_tag_no_use_secondary \
  no_crash_bed_bedgraph \
  secondary_pct_correct \
  json_depth_range \
  warn_flen_populated \
  nested_read_coverage \
  gzopen_fail_graceful \
  bed_ins_coordinates \
  basic_read_counts \
  mapq_filter \
  mapq_threshold_change \
  flen_filter \
  flen_min_override \
  mito_classification \
  use_secondary \
  use_dups \
  pe_fragment_count \
  fragment_avg_size \
  json_valid \
  no_output_flag \
  output_ext \
  save_bam_passthrough \
  bedgraph_output \
  multiple_bams \
  multithreaded_same_output \
  fast_mode \
  blacklist_filtering \
  target_list_filtering \
  plastid_classification \
  max_qlen_filter \
  max_flen_filter \
  frip_calculation \
  quant_output \
  json_multiple_reports \
  strict_preset \
  bed_output \
  qlen_stats_uniform \
  flen_stats_uniform \
  mapq_avg_mixed \
  read_depth_calculation \
  gc_content_known \
  nrf_calculation \
  qlen_mixed_stats \
  mapq_stats_bimodal \
  mapq_percentiles \
  flen_mixed_stats \
  flen_percentiles \
  depths_stats_basic \
  depths_max_pileup \
  gc_mixed_extremes \
  unfiltered_secondary_counts \
; do
  run_test "$t"
done

echo
echo "================================"
echo "Results: $PASS passed, $FAIL failed out of $((PASS + FAIL)) tests"
echo

[ "$FAIL" -eq 0 ]
