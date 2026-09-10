# Read sampling and target sequence size

The default `getfastq --max_bp` is 999,999,999,999,999 bp. Ordinary runs therefore
extract all available spots before the configured filters. The considerations
below matter when the total input exceeds an explicitly smaller target.

`max_bp` is an approximate total target, divided between runs. Average spot
length estimates the initial number of spots. Filtering and a compensatory
second round can leave output below or above the target. It is not a strict
output cap, nor a promise of equal bases, fragments, independent molecules, or
statistical power. For example, 100 million bp corresponds to one million
100-base single-end records, but about 333,333 pairs of 2 x 150 bases.

## Methods

The default `--sampling_method contiguous` retains the previous behavior. A
limited run usually starts at spot 10,000 and takes a continuous interval;
smaller inputs use a tail interval or all spots. A second round continues beyond
the first interval. Input ordering may be associated with lanes, tiles, time,
quality, or expression composition. A continuous interval is not a random
sample. The magnitude of this effect on biological results has not been
established here.

Opt in to reproducible, non-replacement random spot sampling:

```bash
amalgkit getfastq --out_dir ./work --max_bp 1000000000 --sampling_method random --sampling_seed 42
```

A spot is one single-end record or an aligned pair. Both mates are selected
together. Sampling occurs **before** minimum-length filtering, fastp, and the
optional rRNA/contamination filters. If either mate is too short, the entire
selected pair is rejected. `--treat_identical_paired_as_single yes` is currently
incompatible with random sampling; preserve both mates for this mode.

The versioned SHA-256 random stream and partial Fisher-Yates permutation select
each fixed-size sample with equal spot inclusion probabilities. Further rounds
continue the same permutation without replacement. Adaptive compensation based
on filtering is distinct from fixed-size sampling: its final inclusion
probabilities must not be interpreted as a fixed-size uniform sample.

Within a round, output remains in input order; second-round output is appended.
The same validated input order, run ID, seed, algorithm version, budget, and
filter configuration reproduce the same selected spots. Changing the order of
metadata rows or the number of workers does not change the selection. gzip
headers, timing fields, and filesystem identities need not be byte-identical.
Switching between SRA and original FASTQ sources may change order or headers and
is not guaranteed to reproduce the same sample.

## Inputs and cost

SRA (including the ENA/Trace original-FASTQ path) and public GSA use the random
mode when requested. Private files remain exempt unless explicitly enabled:

```bash
amalgkit getfastq --out_dir ./work --max_bp 1000000000 --sampling_method random --sampling_seed 42 --sampling_private yes
```

In random mode only participating runs share the budget and contribute to its
attainment rate. Private source files are read, not overwritten. Keep them outside
their managed getfastq run directory; this is checked before resume cleanup. Input record or
pair counts must agree with the declared spot count. Malformed FASTQ, mate ID
or order mismatch, and unsupported SRA layouts that do not yield one record or
pair per declared spot fail validation instead of silently changing the sampling
population. Update incorrect private metadata through `integrate` first.

Random sampling validates the full population while avoiding redundant intermediates:

- Every extraction round scans the complete input, including unselected reads.
- SRA is fully expanded, or original FASTQ is downloaded and validated. SRA
  expansion is sampled before compression; only selected output is compressed.
  GSA streams directly from its validated cache while holding the input lock,
  without generating a full temporary FASTQ. Each second round scans again.
- Ordinary gzip does not provide general indexed random access. A small sample
  can still require nearly complete decompression.
- Selection state uses memory proportional to the number of permutation ranks
  visited, except for the all-spots first-round shortcut. High sampling fractions
  on very large inputs can require substantial memory.
- Random private-input resume fingerprints hash the source files, so skipping
  extraction still involves reading them.

A small random sample can still cost substantially more than contiguous
extraction because all candidates are read and checked. The default remains
contiguous pending representative biological and resource evaluation.

## Provenance and resume

Each random run writes `getfastq_sampling.json`. Its sampling schema version 1
records the algorithm, run ID, seed, candidate counts and measured bases,
normalized input FASTQ SHA-256, budget provenance, round rank ranges, selected
position digests, and pre/post-minimum-length counts. Rank ranges plus the
versioned algorithm regenerate the selected input positions. The digest uses
canonical LF-terminated FASTQ records, so a missing final newline does not
corrupt concatenated rounds.

`getfastq_stats.tsv` and merged metadata carry `getfastq_sampling_*` fields for
method, seed, algorithm, input digest, candidate spots/bases, and round ranges.
This namespace preserves the existing `sampling_seed` used for selecting runs.
`spot_start_*` and `spot_end_*` are zero (not applicable) for random sampling;
consumers must use the permutation-rank fields and manifest instead. `num_dumped`
and `bp_dumped` count selected candidates, not full materialization I/O;
`num_written` and `bp_written` are after minimum-length filtering. Later filter
counters retain their existing meanings. `bp_specified_for_extraction` and
`bp_still_available` use measured selected/remaining bases in random mode.

The existing resume schema remains compatible for contiguous results. Random
sampling adds a versioned fingerprint component and protects the sampling
manifest with SHA-256. Sampling fingerprint version 2 also rejects results
created before whole-input validation required nonempty IDs and sequences.
Those random runs are recomputed once; contiguous resume is unchanged.
It cannot adopt unmarked legacy output. A changed seed,
algorithm, budget, participating run set, private source contents, or manifest
invalidates affected processing. An interrupted second round restarts the run
under the existing state machine. A zero-yield second round retains the first
round's output; a zero-yield first round retains the existing stop behavior.

Changing sampled reads requires rerunning quantification and any downstream
aggregates that use those reads. Previously complete inputs do not require
scientific recalculation solely because this option was added.

## Validation before changing the default

Regression tests cover selection prefixes, non-replacement, all-spots and tiny
inputs, paired IDs/order, variable lengths, whole-input corruption detection,
measured counts, source preservation, manifest invalidation, round interruption,
and worker/metadata ordering. A fixed-seed synthetic check compares inclusion
frequencies and block composition with binomial/hypergeometric expectations.
These checks establish selection behavior, not biological equivalence to full
data or equal detection power.

Before changing the default, compare full data, contiguous extraction, random
sampling, and seeded stratified/distributed alternatives over multiple seeds.
Include gradients, blocks, periodic ordering, short/long and single/paired inputs.
Measure transcript abundance error, mapping rate, low-expression dropout and
numbers of detected genes against an independent fixed-size random reference.
Distinguish unavoidable information loss from selection bias. Predeclare
scientific tolerances using reference-sampling variation.

Separately measure cold/warm cache runs at 1%, 10%, 50%, and 100% fractions for
SRA, original FASTQ, GSA, and private inputs: wall/CPU time, maximum RSS, bytes
read/written/downloaded, and peak temporary storage, including compensation.
Compare implementations on identical selected sets and equivalent decompressed
outputs. Stratified/block access may reduce seeks on supported sources, but
periodic structure can defeat simple evenly spaced sampling; do not enable it
based on I/O measurements alone. Broader biological comparisons and cold-cache,
network-inclusive measurements remain prerequisites to any default change. Local synthetic timing
measurements alone do not establish biological equivalence.
