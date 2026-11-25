# No Hydra overrides

Staged in place of `assets/rfd4_proteina/hydra_overrides` when `--rfd4_metrics_config_dir false`
is set. It holds no Hydra configs, so adding it to `hydra.searchpath` changes nothing and
`rfd4-generate` runs with the container's own configs.
