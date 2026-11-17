These are empty dummy files used in cases where a process expects a file input _optionally_ (in this case, the `BOLTZ` multiple sequence alignments).

`null` or non-existent `file()` or `path()` objects are not allowed in Nextflow, so we use these as 'flag file' placeholders. Internally, the logic uses the file name to determine what to do.

None of this seems ideal, but the alternative would seem to be duplicating the process definitions with different combinations of path inputs.

`combine_placeholder_*` files are distinct zero-byte placeholders for `COMBINE_RFD3_SCORES`: Nextflow rejects staging the same path (`empty`) more than once per task even when `stageAs` differs, so each optional input uses a unique basename.

`no_refold_target.fasta` means "no `--refold_target_fasta`": `REFOLD_PREPARE_INPUTS` then reads the target back out of the design structures. It has its own name rather than reusing `empty.fa` so that a user's own file called `empty.fa` cannot be mistaken for the placeholder, which is decided on the staged basename.

`empty_pae.npz` means "no PAE matrix for this Boltz sample". `bin/parse_boltz_confidence.py --pae-npz` treats a zero-byte file as absent and leaves the `pae` column empty, rather than failing the task.
