# BoltzGen protein binder design example.

To run locally on a single GPU:

```bash
./run.sh
```

To run _locally_ on a dual-GPU machine:

```bash
./run-dual-gpu.sh
```

After a run, you may like to tweak the filtering of the final designs output to `results/boltzgen/filtered/`:

```bash
./run-filter.sh
```