blastmetag pseudo-code:

Configuration and setup
    Define input file paths:
        TSV: core OTUs table.
        ZOTUS_FA: fasta with zotu sequences.
        REF_FA: cleaned UNITE reference fasta.
        OLD_DB: old taxonomy database fasta.
        OUTPUT: path for new updated database fasta.

    Define thresholds:
        G_CUT: identity cutoff for genus-level.
        S_CUT: identity cutoff for species-level.
        MIN_COV: minimum alignment coverage fraction.
        MAX_SAMPLES: maximum number of samples (default 899 if not set).

    Create a temporary working directory and register cleanup on exit.
        Inside TMPDIR, define paths for all intermediate files
        (filtered TSVs, selected FASTA, BLAST outputs, taxa lists, final DB, etc.).


Step 1: filter OTU table
    Read TSV header and write it to SELECTED_TSV.
    For each data row in TSV:
        Let f = relative frequency column.
        Let s = number of samples column.
        Keep the row if:
            f ≥ 0.001, or
            0.0005 ≤ f < 0.001 and s ≥ MAX_SAMPLES / 2, or
            0.00025 ≤ f < 0.0005 and s ≥ MAX_SAMPLES / 4.
    Save these kept rows (including header) to SELECTED_TSV.
    From SELECTED_TSV, create NONSPEC_TSV:
        Copy header.
        For each subsequent row, inspect taxonomy column.
        Keep only rows where no species label of the form s:... appears.
        Result: OTUs with no species‑level annotation.


Step 2: extract zotu sequences
    From NONSPEC_TSV, collect unique IDs from the first column into ZOTU_IDS_FILE.
    Read ZOTUS_FA as FASTA:
        For each record:
            Extract sequence identifier from header (first whitespace‑separated token).
            If the ID is present in ZOTU_IDS_FILE, copy this record into ZOTUS_SELECTED_FA.
    Result: FASTA file of only those zotus that lack species annotation but pass abundance filters.


Step 3: run BLAST with coverage information
    Build a nucleotide BLAST database from REF_FA with makeblastdb, storing it under BLAST_DB_PREFIX.
    Run blastn:
        Query = ZOTUS_SELECTED_FA.
        Database = the newly built UNITE DB.
        Output format = tabular with fields:
        qseqid, stitle, pident, length, bitscore, evalue, qlen.
        Use up to CPU_THREADS threads and allow up to 500 targets per query.
    Save BLAST tabular output to BLAST_OUT.
​

Step 4: keep best hits under cutoffs
    For each line in BLAST_OUT:
        Extract:
            q: query ID.
            st: subject title (full UNITE taxonomy string).
            pid: percent identity.
            alen: alignment length.
            bits: bitscore.
            qlen: query length.

        If qlen == 0, skip.
        Compute cov = alen / qlen.
        If cov < MIN_COV, skip hit.
        If pid ≥ S_CUT:
            Maintain, for each query, the best species‑level hit (≥S_CUT) by:
                Higher identity first.
                If ties in identity, higher bitscore.
        If pid ≥ G_CUT:
            Maintain, for each query, the best genus‑level hit (≥G_CUT) similarly.

    After processing all hits:
        Output, for each query, the best species‑level hit if it exists.
        If no species hit exists but a genus‑level hit exists, output that genus‑level hit.
    Save these best hits to BEST_HITS.
    ​

Step 5: derive genera and species lists
    For each line in BEST_HITS:
        Take subject title string.
        Parse taxonomy to find:
            g:<genus> tag if present.
            s:<species> tag if present.

        If pid ≥ S_CUT and a species name exists:
            Append species name (no prefix) to SPECIES_FILE.
        Else if pid ≥ G_CUT and a genus name exists:
            Append genus name to GENERA_FILE.

    Sort both GENERA_FILE and SPECIES_FILE in place to keep only unique names.


Step 6: extract matching reference sequences and build new DB
    Concatenate unique genera and species names into TAXA_SEARCH.
    For each taxon name in TAXA_SEARCH:
        Generate two search patterns:
            "g:<name>;"
            "s:<name>;"
        Write all patterns to patterns.txt.

    Scan the reference fasta REF_FA:
        For any header line matching any pattern (case‑insensitive, fixed string):
            Output the header and its following sequence line into FINAL_DB.

    In FINAL_DB, normalize headers:
        Replace any header like >seqNNNN with a shorter > plus rest to strip numeric prefix.

    Write the result to OUTPUT as the candidate updated database.
    Deduplicate sequences in OUTPUT:
        Treat FASTA as records.
        For each record:
        Extract header and sequence.
        If sequence not seen before:
        Record it as seen and keep this record.
        If sequence already present, discard duplicate.

    Overwrite OUTPUT with the deduplicated FASTA.
    Print total number of sequences now in OUTPUT.


Step 7: compare against previous database and report gains
    Extract raw sequences (no headers) from OLD_DB into old_seqs.txt:
        For each FASTA record, output only the sequence part, removing carriage returns.
    Do the same for OUTPUT into new_seqs.txt.
    Compute new unique sequences:
        Set new_unique_seqs.txt to sequences that appear in new_seqs.txt but not in old_seqs.txt.
        Count lines in new_unique_seqs.txt as new_total.

    If new_total == 0:
        Report that no new unique sequences were added.

    Else:
        Report how many unique sequences were added.
        Build a mapping table seq2header.tsv from sequence string to header in OUTPUT.
        Using new_unique_seqs.txt and this mapping, obtain headers of newly added sequences into new_headers.txt.

        From these headers:
            Count how many distinct genera appear (by tags g:<genus>).
            Count how many distinct species appear (by tags s:<species>).

        Print:
            Number of new genera added.
            Number of new species added.

    Finally, print the path to the new database file and a completion message.