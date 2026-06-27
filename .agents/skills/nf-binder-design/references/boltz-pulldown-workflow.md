# Boltz Pulldown Workflow — Full Reference

## Overview

An AlphaPulldown-like protocol using Boltz-2 for multimer predictions from sequence. Assess binding interactions between designed binders and targets.

## Parameters

Use `--method boltz_pulldown --help` for the full list.

### MSA Options

By default, **no MSAs are used**.

- `--create_target_msa=true` — Generate MSAs for targets (recommended)
- `--create_binder_msa=false` — Skip binder MSAs (default; no homologues for de novo designs)

### MSA Generation

**Remote server**: `--use_msa_server` (public ColabFold MMSeqs2 server)

**Local databases**:
```bash
--uniref30 /path/to/uniref30
--colabfold_envdb /path/to/colabfold_envdb
```

### Custom Boltz Arguments

Via `ext.args` in `nextflow.config`:

```groovy
process {
    withName: BOLTZ {
        ext.args = "--accelerator cpu"
    }
}
```

## Output

Default: `results/boltz_pulldown/`

- `boltz_pulldown.tsv` — summary table
- `boltz_pulldown_report.html` — report with ipTM score plots
- Boltz predictions with scores and structures
