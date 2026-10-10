Empty placeholders with the nine names AlphaFold3 validates against `--db_dir`
before it reads its input. They exist only so `--af3_run_data_pipeline true` can
be compile-tested (`-preview`, no process runs) without the real ~630 GB set.

Do NOT point a real run at this directory.
