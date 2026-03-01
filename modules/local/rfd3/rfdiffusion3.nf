process RFDIFFUSION3 {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/rfdiffusion3", pattern: 'output/*.cif.gz', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3/rfdiffusion3", pattern: 'output/*.json', mode: 'copy'

    input:
    path config_json
    path input_pdb
    val design_name
    val batch_size
    val n_batches
    val design_startnum
    val step_scale
    val gamma_0
    val extra_args

    output:
    path 'output/*.cif.gz', emit: cifs
    path 'output/*.json', emit: json_metrics
    path 'rfd3_*_batch*_backbone.tsv', emit: scores

    script:
    // Always set global_prefix so all batches share the {design_name}_ prefix.
    // rfd3 appends the config key after global_prefix, so we include a trailing
    // separator to keep the batch index visually distinct.
    def global_prefix = "global_prefix=${design_name}_b${design_startnum}_"
    def extra_args_str = extra_args ?: ''
    def uid = design_name.replaceFirst(/^rfd3_/, '')
    """
    set -euo pipefail

    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

    # Find least-used GPU and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    # Rewrite input paths in config to use staged basenames (copy to avoid modifying staged symlink)
    ${projectDir}/bin/rfd3/stage_rfd3_config.py stage ${config_json} -o ${design_name}.json

    # Run RFDiffusion3
    rfd3 design \
        out_dir=output \
        inputs=${design_name}.json \
        diffusion_batch_size=${batch_size} \
        n_batches=${n_batches} \
        inference_sampler.step_scale=${step_scale} \
        inference_sampler.gamma_0=${gamma_0} \
        ${global_prefix} \
        ${extra_args_str} \
        ${task.ext.args ?: ''}

    for j in output/*.json; do
      [[ -f "\$j" ]] || continue
      batch=\$(basename "\$j" .json | sed -n 's/.*batch\\([0-9]*\\).*/\1/p; s/.*_b\\([0-9]*\\)_.*/\1/p' | head -1)
      [[ -z "\$batch" ]] && batch=0
      python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rfdiffusion3 "\$j" -o rfd3_${uid}_batch\${batch}_backbone.tsv
    done
    """
}

/*
Usage:

https://rosettacommons.github.io/foundry/models/rfd3/index.html
https://github.com/RosettaCommons/foundry/tree/production/models/rfd3

$ python models/rfd3/src/rfd3/run_inference.py --help

run_inference is powered by Hydra.

== Configuration groups ==
Compose your configuration from those groups (group=option)

callbacks: design_callbacks, metrics_logging, train_logging
dataloader: default, fast
datasets: design_base
datasets/conditions: dna_condition, island, ppi, sequence_design, tipatom, unconditional
datasets/train: rfd3_monomer_distillation
datasets/train/pdb: af3_train_interface, af3_train_pn_unit, base, base_no_weights, base_transform_args, na_complex_distillation, pdb_base, rfd3_train_interface, rfd3_train_pn_unit
datasets/val: bcov_ppi_easy_medium, design_validation_base, dna_binder_design5, dna_binder_long, dna_binder_short, indexed, mcsa_41, mcsa_41_short_rigid, ppi_inference, sm_binder_hbonds, sm_binder_hbonds_short, unconditional, unconditional_deep, unindexed
datasets/val/val_examples: bcov_ppi_easy_medium_with_ori, bcov_ppi_easy_medium_with_ori_spoof_helical_bundle, bcov_ppi_easy_medium_with_ori_varying_lengths, bpem_ori_hb
debug: default, train_specific_examples
experiment: debug, pretrain, test-uncond, test-unindexed
inference_engine: base, dev, rfdiffusion3
logger: csv, default, wandb
model: rfd3_base
model/components: ema, rfd3_net
model/optimizers: adam
model/samplers: edm, symmetry
model/schedulers: af3
paths: default
paths/data: default
trainer: cpu, ddp, rfd3_base, xpu
trainer/loss/losses: diffusion_loss, sequence_loss
trainer/metrics: design_metrics


== Config ==
Override anything in the config (foo.bar=value)

ckpt_path: rfd3
num_nodes: 1
devices_per_node: 1
verbose: false
seed: null
inputs: ???
out_dir: ???
_target_: rfd3.engine.RFD3InferenceEngine
json_keys_subset: null
skip_existing: true
specification: {}
diffusion_batch_size: 8
n_batches: 1
inference_sampler:
  kind: default
  cfg_features:
  - active_donor
  - active_acceptor
  - ref_atomwise_rasa
  use_classifier_free_guidance: false
  cfg_t_max: null
  cfg_scale: 1.5
  center_option: all
  s_trans: 1.0
  inference_noise_scaling_factor: 1.0
  allow_realignment: false
  num_timesteps: 200
  step_scale: 1.5
  noise_scale: 1.003
  p: 7
  gamma_0: 0.6
  gamma_min: 1.0
  s_jitter_origin: 0.0
cleanup_guideposts: true
cleanup_virtual_atoms: true
read_sequence_from_sequence_head: true
output_full_json: true
global_prefix: null
dump_prediction_metadata_json: true
dump_trajectories: false
align_trajectory_structures: false
prevalidate_inputs: false
low_memory_mode: false


Powered by Hydra (https://hydra.cc)
Use --hydra-help to view Hydra specific help
*/