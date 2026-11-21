# BindCraft2 test run (PD-L1, small)

This minimal `--method bc2` example runs four trajectories against PD-L1 using the default stage budget.

See `examples/pdl1-bc2/` for a larger example covering config mode and VHH modality.

```bash
./run.sh

# or for a dual-GPU workstation:
# ./run-dual-gpu.sh
```

Accepted designs are written as mmCIF files in `results/bc2/ranked/`, with `results/bc2/summary.csv` and `results/bc2/bc2_report.html` summarising all four batches. With only four trajectories, getting zero accepted designs is normal even on PD-L1 and does not indicate an error.

## Test settings (not for production)

| Lever | Typical Production | Here |
|---|---|---|
| Trajectories | 100s to 1000s | 4 |
| Binder length | 60-150 | 55-65 |

Short binders keep the complex near 180 residues, which uses less VRAM, so the example should run in less than ~6 GB VRAM.

To make it faster at the cost of the accept rate, cut the gradient stage budget, which
defaults to screen 50, refine 25, anneal 45, harden 5, mutate 15 = 140 steps per trajectory:

```bash
--bc2_set 'screen_steps=10;refine_steps=5;anneal_steps=10;harden_steps=2;mutate_steps=3'
```

Note that this is only a pipeline-plumbing test, not setting you'd expect to work for a successful design run.
