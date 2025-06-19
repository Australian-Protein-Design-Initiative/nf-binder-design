These are empty dummy files used in cases where a process expects a file input _optionally_ (in this case, the `BOLTZ` multiple sequence alignments).

`null` or non-existent `file()` or `path()` objects are not allowed in Nextflow, so we use these as 'flag file' placeholders. Internally, the logic uses the file name to determine what to do.

None of this seems ideal, but the alternative would seem to be duplicating the process definitions with different combinations of path inputs.
