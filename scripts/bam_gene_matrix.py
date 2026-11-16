# Ver 0.0.6
# Giang le

#!/usr/bin/env python3

import sys
import os
from collections import Counter, defaultdict
import difflib

import pysam


def usage():
    print(
        "\nUsage:\n"
        "  python bam_gene_matrix.py <input.bam> <reference.fa> <chrom> "
        "<start_1based> <end_1based> <output_prefix>\n\n"
        "Example:\n"
        "  python bam_gene_matrix.py sample_hap1.bam sample_hap1_TRA_contig.fa "
        "hap1_sample 94610322 94610597 sample_TRA\n"
    )
    sys.exit(1)


def print_contig_diagnostics(
    chrom,
    bam,
    ref,
    bam_path,
    ref_fasta_path,
    start0,
    end0,
):
    bam_abs_path = os.path.realpath(bam_path)
    fasta_abs_path = os.path.realpath(ref_fasta_path)

    bam_contigs = list(bam.references)
    fasta_contigs = list(ref.references)

    print()
    print("ERROR: Requested contig was not found or could not be fetched.")
    print()

    print("Input files:")
    print(f"  BAM:   {bam_abs_path}")
    print(f"  FASTA: {fasta_abs_path}")
    print()

    print("Requested region:")
    print(f"  Contig:        {chrom}")
    print(f"  Start 1-based: {start0 + 1}")
    print(f"  End 1-based:   {end0}")
    print()

    print("Contig present in BAM?")
    print(f"  {'YES' if chrom in bam_contigs else 'NO'}")
    print()

    print("Contig present in FASTA?")
    print(f"  {'YES' if chrom in fasta_contigs else 'NO'}")
    print()

    bam_matches = difflib.get_close_matches(
        chrom, bam_contigs, n=20, cutoff=0.2
    )

    fasta_matches = difflib.get_close_matches(
        chrom, fasta_contigs, n=20, cutoff=0.2
    )

    print("Closest matching BAM contigs:")
    if bam_matches:
        for contig in bam_matches:
            print(f"  {contig}")
    else:
        print("  No close BAM contig matches found.")
    print()

    print("Closest matching FASTA contigs:")
    if fasta_matches:
        for contig in fasta_matches:
            print(f"  {contig}")
    else:
        print("  No close FASTA contig matches found.")
    print()

    print("All BAM contigs:")
    for contig in bam_contigs:
        print(f"  {contig}")
    print()

    print("All FASTA contigs:")
    for contig in fasta_contigs:
        print(f"  {contig}")
    print()

    sys.exit(2)


def validate_contig(
    chrom,
    bam,
    ref,
    bam_path,
    ref_fasta_path,
    start0,
    end0,
):
    bam_contigs = set(bam.references)
    fasta_contigs = set(ref.references)

    if chrom not in bam_contigs or chrom not in fasta_contigs:
        print_contig_diagnostics(
            chrom,
            bam,
            ref,
            bam_path,
            ref_fasta_path,
            start0,
            end0,
        )


def get_read_name(read):
    read_name = read.query_name

    if read.is_read1:
        read_name += "/1"
    elif read.is_read2:
        read_name += "/2"

    return read_name


def build_read_calls(read, seq, start0, end0, gene_len):
    """
    Build reference-position calls and collect insertions.

    Reference symbols:
        A/C/G/T/N = nucleotide
        -         = deletion
        ~         = reference skip
        .         = uncovered

    Insertions are anchored after a reference position.

    Example:
        Reference: A T G C A
        Read:      A T G G G C A

        Insertion GG is anchored after reference position 3.

    Only insertions strictly between positions inside the
    requested region are included in the expanded matrix.
    """

    calls = ["."] * gene_len
    insertions = {}

    ref_pos = read.reference_start
    query_pos = 0

    if read.cigartuples is None:
        return calls, insertions

    for operation, length in read.cigartuples:

        # M, =, X
        if operation in (0, 7, 8):

            for j in range(length):
                current_ref = ref_pos + j
                current_query = query_pos + j

                if start0 <= current_ref < end0:
                    gene_index = current_ref - start0
                    calls[gene_index] = seq[current_query].upper()

            ref_pos += length
            query_pos += length

        # Insertion
        elif operation == 1:

            insertion_seq = seq[
                query_pos:query_pos + length
            ].upper()

            # Anchor insertion after preceding reference base.
            anchor_ref = ref_pos - 1

            # Include internal insertion sites.
            if start0 <= anchor_ref < end0 - 1:
                anchor_index = anchor_ref - start0

                if anchor_index in insertions:
                    insertions[anchor_index] += insertion_seq
                else:
                    insertions[anchor_index] = insertion_seq

            query_pos += length

        # Deletion
        elif operation == 2:

            for j in range(length):
                current_ref = ref_pos + j

                if start0 <= current_ref < end0:
                    gene_index = current_ref - start0
                    calls[gene_index] = "-"

            ref_pos += length

        # Reference skip
        elif operation == 3:

            for j in range(length):
                current_ref = ref_pos + j

                if start0 <= current_ref < end0:
                    gene_index = current_ref - start0
                    calls[gene_index] = "~"

            ref_pos += length

        # Soft clipping
        elif operation == 4:
            query_pos += length

        # Hard clipping or padding
        elif operation in (5, 6):
            continue

        else:
            raise ValueError(
                f"Unsupported CIGAR operation: {operation}"
            )

    return calls, insertions


def summarize_read(calls, ref_seq, insertions):
    covered_positions = 0
    matches = 0
    mismatches = 0
    deletions = 0
    reference_skips = 0

    for i, base in enumerate(calls):

        if base == ".":
            continue

        if base == "~":
            reference_skips += 1
            continue

        covered_positions += 1

        ref_base = ref_seq[i]

        if base == "-":
            deletions += 1

        elif base == ref_base:
            matches += 1

        else:
            mismatches += 1

    inserted_bases = sum(
        len(insertion_seq)
        for insertion_seq in insertions.values()
    )

    insertion_events = len(insertions)

    return (
        covered_positions,
        matches,
        mismatches,
        deletions,
        reference_skips,
        inserted_bases,
        insertion_events,
    )


def collect_insertion_columns(all_read_data):
    """
    Determine the maximum insertion length at each
    reference position.

    Example:
        Read 1: GG
        Read 2: GGA
        Read 3: G

    Maximum insertion length = 3.
    """

    max_insertion_lengths = defaultdict(int)

    for read_data in all_read_data:

        insertions = read_data["insertions"]

        for anchor_index, insertion_seq in insertions.items():

            insertion_length = len(insertion_seq)

            if insertion_length > max_insertion_lengths[anchor_index]:
                max_insertion_lengths[anchor_index] = insertion_length

    return dict(max_insertion_lengths)


def build_expanded_columns(gene_len, max_insertion_lengths):
    """
    Build ordered matrix column definitions.

    Each column is:
        ("REF", reference_index, None)

    or:
        ("INS", anchor_index, insertion_base_index)
    """

    columns = []

    for i in range(gene_len):

        columns.append(("REF", i, None))

        max_length = max_insertion_lengths.get(i, 0)

        for insertion_index in range(max_length):
            columns.append(
                ("INS", i, insertion_index)
            )

    return columns


def insertion_site_is_covered(calls, anchor_index):
    """
    A read must have reference alignment information
    on both sides of an insertion site.

    A deletion counts as reference alignment information.
    A reference skip does not.
    """

    if anchor_index + 1 >= len(calls):
        return False

    left_base = calls[anchor_index]
    right_base = calls[anchor_index + 1]

    if left_base in (".", "~"):
        return False

    if right_base in (".", "~"):
        return False

    return True


def build_expanded_read(
    calls,
    insertions,
    columns,
):
    """
    Construct a read row for the expanded matrix.

    In insertion columns:
        A/C/G/T = inserted nucleotide
        -       = no inserted nucleotide
        .       = insertion site not covered
    """

    expanded = []

    for column_type, ref_index, insertion_index in columns:

        if column_type == "REF":

            expanded.append(calls[ref_index])

        else:

            insertion_seq = insertions.get(ref_index, "")

            if not insertion_site_is_covered(calls, ref_index):
                expanded.append(".")

            elif insertion_index < len(insertion_seq):
                expanded.append(insertion_seq[insertion_index])

            else:
                expanded.append("-")

    return expanded


def write_gene_matrix(
    output_path,
    all_read_data,
    ref_seq,
    columns,
):
    with open(output_path, "w") as matrix_fh:

        header = ["read_id"]
        reference_row = ["REFERENCE"]

        for column_type, ref_index, insertion_index in columns:

            if column_type == "REF":

                header.append(
                    f"gene_pos_{ref_index + 1}"
                )

                reference_row.append(ref_seq[ref_index])

            else:

                header.append(
                    f"INS_{ref_index + 1}_{insertion_index + 1}"
                )

                reference_row.append("-")

        matrix_fh.write(
            "\t".join(header) + "\n"
        )

        matrix_fh.write(
            "\t".join(reference_row) + "\n"
        )

        for read_data in all_read_data:

            expanded = build_expanded_read(
                read_data["calls"],
                read_data["insertions"],
                columns,
            )

            matrix_fh.write(
                read_data["read_name"]
                + "\t"
                + "\t".join(expanded)
                + "\n"
            )


def write_read_summary(
    output_path,
    all_read_data,
    ref_seq,
    gene_len,
):
    with open(output_path, "w") as read_fh:

        read_fh.write(
            "read_id\t"
            "covered_positions\t"
            "matches_to_reference\t"
            "mismatches_to_reference\t"
            "deletions\t"
            "reference_skips\t"
            "inserted_bases\t"
            "insertion_events\t"
            "percent_identity_vs_reference\t"
            "percent_identity_including_insertions\t"
            "percent_gene_covered\n"
        )

        for read_data in all_read_data:

            (
                covered,
                matches,
                mismatches,
                deletions,
                reference_skips,
                inserted_bases,
                insertion_events,
            ) = summarize_read(
                read_data["calls"],
                ref_seq,
                read_data["insertions"],
            )

            # Original identity calculation.
            percent_identity = (
                matches / covered * 100
                if covered > 0
                else 0.0
            )

            # Insertion-aware identity.
            identity_denominator = (
                matches
                + mismatches
                + deletions
                + inserted_bases
            )

            percent_identity_with_insertions = (
                matches / identity_denominator * 100
                if identity_denominator > 0
                else 0.0
            )

            percent_gene_covered = (
                covered / gene_len * 100
            )

            read_fh.write(
                f"{read_data['read_name']}\t"
                f"{covered}\t"
                f"{matches}\t"
                f"{mismatches}\t"
                f"{deletions}\t"
                f"{reference_skips}\t"
                f"{inserted_bases}\t"
                f"{insertion_events}\t"
                f"{percent_identity:.2f}\t"
                f"{percent_identity_with_insertions:.2f}\t"
                f"{percent_gene_covered:.2f}\n"
            )


def write_position_summary(
    output_path,
    all_read_calls,
    ref_seq,
    chrom,
    gene_start_1based,
    gene_len,
):
    with open(output_path, "w") as pos_fh:

        pos_fh.write(
            "gene_pos\t"
            "chrom\t"
            "ref_pos_1based\t"
            "ref_base\t"
            "depth_covered\t"
            "A_count\tC_count\tG_count\tT_count\tN_count\tDEL_count\t"
            "A_percent\tC_percent\tG_percent\tT_percent\tN_percent\tDEL_percent\t"
            "major_base\tmajor_base_percent\n"
        )

        for i in range(gene_len):

            ref_pos_1based = gene_start_1based + i
            ref_base = ref_seq[i]

            # Exclude uncovered and reference-skipped positions.
            bases_at_position = [
                calls[i]
                for calls in all_read_calls
                if calls[i] not in (".", "~")
            ]

            counts = Counter(bases_at_position)
            depth = len(bases_at_position)

            if depth > 0:

                a_pct = counts["A"] / depth * 100
                c_pct = counts["C"] / depth * 100
                g_pct = counts["G"] / depth * 100
                t_pct = counts["T"] / depth * 100
                n_pct = counts["N"] / depth * 100
                del_pct = counts["-"] / depth * 100

                major_base, major_count = counts.most_common(1)[0]
                major_pct = major_count / depth * 100

            else:

                a_pct = c_pct = g_pct = t_pct = n_pct = del_pct = 0.0
                major_base = "."
                major_pct = 0.0

            pos_fh.write(
                f"{i + 1}\t"
                f"{chrom}\t"
                f"{ref_pos_1based}\t"
                f"{ref_base}\t"
                f"{depth}\t"
                f"{counts['A']}\t"
                f"{counts['C']}\t"
                f"{counts['G']}\t"
                f"{counts['T']}\t"
                f"{counts['N']}\t"
                f"{counts['-']}\t"
                f"{a_pct:.2f}\t"
                f"{c_pct:.2f}\t"
                f"{g_pct:.2f}\t"
                f"{t_pct:.2f}\t"
                f"{n_pct:.2f}\t"
                f"{del_pct:.2f}\t"
                f"{major_base}\t"
                f"{major_pct:.2f}\n"
            )


def write_insertion_summary(
    output_path,
    all_read_data,
    chrom,
    gene_start_1based,
):
    """
    Summarize insertion sequences at each reference anchor.

    Frequency denominator:
        Reads spanning both sides of the insertion site.
    """

    insertion_counts = defaultdict(Counter)
    insertion_depth = Counter()

    all_anchors = set()

    for read_data in all_read_data:
        all_anchors.update(read_data["insertions"].keys())

    for anchor_index in sorted(all_anchors):

        for read_data in all_read_data:

            calls = read_data["calls"]
            insertions = read_data["insertions"]

            if not insertion_site_is_covered(calls, anchor_index):
                continue

            insertion_depth[anchor_index] += 1

            insertion_seq = insertions.get(anchor_index)

            if insertion_seq:
                insertion_counts[anchor_index][insertion_seq] += 1

    with open(output_path, "w") as ins_fh:

        ins_fh.write(
            "chrom\t"
            "gene_pos_after\t"
            "ref_pos_after_1based\t"
            "insertion_sequence\t"
            "insertion_length\t"
            "supporting_reads\t"
            "depth_spanning_site\t"
            "insertion_percent\n"
        )

        for anchor_index in sorted(insertion_counts):

            depth = insertion_depth[anchor_index]

            for insertion_seq, count in insertion_counts[
                anchor_index
            ].most_common():

                insertion_percent = (
                    count / depth * 100
                    if depth > 0
                    else 0.0
                )

                ins_fh.write(
                    f"{chrom}\t"
                    f"{anchor_index + 1}\t"
                    f"{gene_start_1based + anchor_index}\t"
                    f"{insertion_seq}\t"
                    f"{len(insertion_seq)}\t"
                    f"{count}\t"
                    f"{depth}\t"
                    f"{insertion_percent:.2f}\n"
                )


def main():

    if len(sys.argv) != 7:
        usage()

    bam_path = sys.argv[1]
    ref_fasta_path = sys.argv[2]
    chrom = sys.argv[3]

    gene_start_1based = int(sys.argv[4])
    gene_end_1based = int(sys.argv[5])

    output_prefix = sys.argv[6]

    # Preserve original coordinate handling.
    # Input coordinates are already adjusted by the user.
    start0 = gene_start_1based
    end0 = gene_end_1based

    gene_len = end0 - start0

    if gene_len <= 0:

        print()
        print("ERROR: End coordinate must be greater than start coordinate.")
        print()
        print(f"Start coordinate: {gene_start_1based}")
        print(f"End coordinate:   {gene_end_1based}")
        print()

        sys.exit(1)

    matrix_out = f"{output_prefix}.gene_read_matrix.tsv"

    read_summary_out = (
        f"{output_prefix}.read_vs_reference.tsv"
    )

    position_summary_out = (
        f"{output_prefix}.position_base_percent.tsv"
    )

    insertion_summary_out = (
        f"{output_prefix}.insertion_summary.tsv"
    )

    all_read_data = []
    all_read_calls = []

    total_reads_seen = 0
    reads_removed_full_deletion = 0
    reads_written = 0

    total_insertion_events = 0
    total_inserted_bases = 0

    try:

        with pysam.AlignmentFile(
            bam_path, "rb"
        ) as bam, pysam.FastaFile(
            ref_fasta_path
        ) as ref:

            validate_contig(
                chrom,
                bam,
                ref,
                bam_path,
                ref_fasta_path,
                start0,
                end0,
            )

            try:

                ref_seq = ref.fetch(
                    chrom,
                    start0,
                    end0,
                ).upper()

            except ValueError:

                print_contig_diagnostics(
                    chrom,
                    bam,
                    ref,
                    bam_path,
                    ref_fasta_path,
                    start0,
                    end0,
                )

            if len(ref_seq) != gene_len:
                raise ValueError(
                    "Requested reference interval exceeds the "
                    "available FASTA sequence or has an "
                    "unexpected length. "
                    f"Expected {gene_len} bases, "
                    f"received {len(ref_seq)}."
                )

            try:

                reads = bam.fetch(
                    chrom,
                    start0,
                    end0,
                )

            except ValueError:

                print_contig_diagnostics(
                    chrom,
                    bam,
                    ref,
                    bam_path,
                    ref_fasta_path,
                    start0,
                    end0,
                )

            # First pass:
            # Collect reference calls and insertion events.

            for read in reads:

                if read.is_unmapped or read.query_sequence is None:
                    continue

                total_reads_seen += 1

                read_name = get_read_name(read)

                calls, insertions = build_read_calls(
                    read,
                    read.query_sequence,
                    start0,
                    end0,
                    gene_len,
                )

                # Preserve original full-deletion filtering.
                if all(base == "-" for base in calls):

                    reads_removed_full_deletion += 1
                    continue

                read_data = {
                    "read_name": read_name,
                    "calls": calls,
                    "insertions": insertions,
                }

                all_read_data.append(read_data)
                all_read_calls.append(calls)

                reads_written += 1

                total_insertion_events += len(insertions)

                total_inserted_bases += sum(
                    len(seq)
                    for seq in insertions.values()
                )

        # Determine insertion columns across all reads.

        max_insertion_lengths = collect_insertion_columns(
            all_read_data
        )

        columns = build_expanded_columns(
            gene_len,
            max_insertion_lengths,
        )

        # Write expanded matrix.

        write_gene_matrix(
            matrix_out,
            all_read_data,
            ref_seq,
            columns,
        )

        # Write read-level statistics.

        write_read_summary(
            read_summary_out,
            all_read_data,
            ref_seq,
            gene_len,
        )

        # Write reference-position statistics.

        write_position_summary(
            position_summary_out,
            all_read_calls,
            ref_seq,
            chrom,
            gene_start_1based,
            gene_len,
        )

        # Write insertion statistics.

        write_insertion_summary(
            insertion_summary_out,
            all_read_data,
            chrom,
            gene_start_1based,
        )

        print()
        print("SUCCESS")
        print()

        print(f"Wrote: {matrix_out}")
        print(f"Wrote: {read_summary_out}")
        print(f"Wrote: {position_summary_out}")
        print(f"Wrote: {insertion_summary_out}")

        print()

        print("Read filtering summary:")
        print(f"  Reads examined:              {total_reads_seen}")
        print(f"  Reads written:               {reads_written}")
        print(f"  100% deletion reads removed: {reads_removed_full_deletion}")

        print()

        print("Insertion summary:")
        print(f"  Total insertion events:      {total_insertion_events}")
        print(f"  Total inserted bases:        {total_inserted_bases}")
        print(f"  Insertion sites:             {len(max_insertion_lengths)}")

        print()

        print("Matrix dimensions:")
        print(f"  Reference positions:         {gene_len}")
        print(f"  Insertion columns:           {len(columns) - gene_len}")
        print(f"  Total matrix columns:        {len(columns)}")

        print()

    except FileNotFoundError as e:

        print()
        print("ERROR: Input or output path was not found.")
        print()
        print(f"Details: {e}")
        print()

        print("Input files:")
        print(f"  BAM:   {os.path.realpath(bam_path)}")
        print(f"  FASTA: {os.path.realpath(ref_fasta_path)}")
        print()

        sys.exit(1)

    except OSError as e:

        print()
        print("ERROR: File or index problem.")
        print()
        print(f"Details: {e}")
        print()

        print("Possible causes:")
        print("  1. BAM file is missing.")
        print("  2. BAM index is missing.")
        print("  3. FASTA file is missing.")
        print("  4. FASTA index .fai is missing.")
        print("  5. Output directory does not exist.")
        print()

        print("Input files:")
        print(f"  BAM:   {os.path.realpath(bam_path)}")
        print(f"  FASTA: {os.path.realpath(ref_fasta_path)}")
        print()

        sys.exit(1)

    except ValueError as e:

        print()
        print("ERROR: Invalid alignment or reference interval.")
        print()
        print(f"Details: {e}")
        print()

        sys.exit(1)


if __name__ == "__main__":
    main()