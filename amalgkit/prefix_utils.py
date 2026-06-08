from bisect import bisect_left


def find_prefixed_entries(entries, prefix, entries_sorted=False):
    if entries_sorted:
        left = bisect_left(entries, prefix)
        right = bisect_left(entries, prefix + '\uffff')
        matched = []
        for i in range(left, right):
            entry = entries[i]
            if entry.startswith(prefix):
                matched.append(entry)
        return matched
    return sorted([
        entry for entry in entries
        if entry.startswith(prefix)
    ])


def find_species_prefixed_entries(entries, species_prefix, entries_sorted=False):
    # Match exact species prefix and only delimiter-based variants (e.g., ".idx", "_k31.idx", "-v1.fa").
    matched = find_prefixed_entries(entries, species_prefix, entries_sorted=entries_sorted)
    allowed = []
    for entry in matched:
        if entry == species_prefix:
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '.'):
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '_'):
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '-'):
            allowed.append(entry)
            continue
    return allowed


def find_run_prefixed_entries(entries, run_id, entries_sorted=False):
    # Match exact run ID and only known output delimiters (e.g., ".fastq.gz", "_1.fastq.gz").
    matched = find_prefixed_entries(entries, run_id, entries_sorted=entries_sorted)
    allowed = []
    for entry in matched:
        if entry == run_id:
            allowed.append(entry)
            continue
        if entry.startswith(run_id + '.'):
            allowed.append(entry)
            continue
        if entry.startswith(run_id + '_'):
            allowed.append(entry)
            continue
    return allowed
