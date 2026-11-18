/*
Which engines take a list of seeds in one invocation, and how to fold the
per-seed job fan-out into one job for them.

--n_predictions is split into jobs of <engine>_batch_size samples each, and job i
runs seed base+i (see each subworkflow's foldPredictionBatches). With a batch size
of 1 - what refold scoring asks for, so that samples differ by more than diffusion
noise - that is N jobs of one sample, and N model loads. On a small complex the
load is most of the wall time.

Some engines will take every seed in one process, keeping each seed's outputs
apart themselves, so N replicates cost one load. Where that holds, the N jobs
collapse into one carrying the whole seed list.

Membership is a claim about the engine, checked against its source rather than its
documentation:

  esmfold2, esmfold2_fast  bin/fold/run_esmfold2.py loops fold() over --seed with
                           the model resident; names files <name>_seed_<S>_sample_<i>_*
  openfold3               runner YAML experiment_settings.seeds is a list; writes
                          output/<name>/seed_<S>/<name>_seed_<S>_sample_<i>_*
  protenix                runner/inference.py loops `for seed in seeds:` over
                          --seeds with the model resident; writes
                          <name>/seed_<N>/predictions/<name>_sample_<r>.cif, which
                          is why PROTENIX_FOLD has to put the seed back into the
                          names it publishes and scores

af3 also takes a seed list (modelSeeds), but by its own parameter: a comma list in
--af3_seeds already means "these seeds, one job", and ALPHAFOLD3_FOLD suppresses
the --n_predictions fan-out when it sees one. It is not collapsed here, because
that would be two mechanisms for the same thing.

boltz and rf3 are absent because they are single-seed per process: boltz calls
seed_everything(seed) once (boltz/main.py), and rf3's inference engine config has
`seed:` as a scalar that it stamps into the output directory name.
*/

def foldSeedsListEngines() {
    return [
        'esmfold2',
        'esmfold2_fast',
        'openfold3',
        'protenix',
    ] as Set
}

// Whether one invocation of this engine takes several seeds.
boolean foldSeedsTakesSeedList(String tool) {
    return foldSeedsListEngines().contains(tool)
}

/*
The seeds one job should run, or null to keep the per-batch fan-out.

Collapsing needs a base seed to count from - without one the jobs are
indistinguishable and running them together would give the same structure N
times - and every batch the same size, so that one sample count describes the
single job that replaces them.
*/
List<Integer> foldSeedsCollapse(String tool, List batches, Integer base_seed, boolean enabled = true) {
    if (!enabled || !foldSeedsTakesSeedList(tool) || base_seed == null) {
        return null
    }
    if (batches.size() < 2 || batches.toSet().size() != 1) {
        return null
    }
    return (0..<batches.size()).collect { i -> base_seed + i }
}
