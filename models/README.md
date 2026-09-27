Put model weights here (or symlinks to them).

You can download weights with the helper scripts:
```bash
./download_rfd_weights.sh
./download_af2_weights.sh
```

## AlphaFold3

AlphaFold3 weights are not included in any container - they are subject to the
[AlphaFold3 Model Parameters Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)
(non-commercial use only, no redistribution outside your organisation). After reading
the terms, download them into `models/alphafold3/` (the default `--af3_model_dir`):

```bash
./download_af3_weights.sh                 # interactive terms prompt
./download_af3_weights.sh -o /shared/af3  # elsewhere; then pass --af3_model_dir /shared/af3
```

The directory must contain exactly one model file (`af3.bin.zst` or `af3.bin`).
