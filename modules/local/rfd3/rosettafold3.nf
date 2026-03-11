process ROSETTAFOLD3 {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/rosettafold3/output", pattern: 'output/*', mode: 'copy', saveAs: { it.replaceFirst('output/', '') }

    input:
    tuple val(meta), path(structure_cif), path(target_msa), path(target_templates), path(rf3_input_json)
    val(uid)

    output:
    path 'output/*', emit: results
    path 'output/*/*_summary_confidences.json', emit: confidence_json
    tuple val(meta), path('output/*/*_model.cif'), emit: refolded_cif
    path 'rfd3_*_rf3_*.tsv', emit: scores
    tuple val(meta), path('rfd3_*_rf3_*.tsv'), emit: scores_with_meta

    script:
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

    mkdir -p output

    rf3 fold \\
        inputs=${rf3_input_json} \\
        out_dir=output \\
        ckpt_path=${params.rf3_ckpt_path} \\
        num_steps=${params.rf3_num_steps} \\
        n_recycles=${params.rf3_n_recycles} \\
        diffusion_batch_size=${params.rf3_diffusion_batch_size} \\
        early_stopping_plddt_threshold=${params.rf3_early_stopping_plddt_threshold} \\
        annotate_b_factor_with_plddt=true \\
        one_model_per_file=true \\
        ${task.ext.args ?: ''}

    full_conf=\$(find output -maxdepth 2 -name '*_confidences.json' ! -name '*summary*' -print | head -1)
    summary_conf=\$(find output -maxdepth 2 -name '*_summary_confidences.json' -print | head -1)
    cif=\$(find output -maxdepth 2 -name '*_model.cif' -print | head -1)
    if [[ -n "\$full_conf" && -n "\$summary_conf" && -n "\$cif" ]]; then
        outdir=\$(dirname "\$cif")
        (cd "\$outdir" && python ${projectDir}/bin/ipsae.py \\
            --format rf3 \\
            --update-summary "\$(basename "\$summary_conf")" \\
            --binder-chain B --target-chain A \\
            "\$(basename "\$full_conf")" "\$(basename "\$cif")" 10 10)
    fi

    if [[ -n "\$summary_conf" ]]; then
        suffix=\$(basename "\$summary_conf" _summary_confidences.json | sed -n 's/.*\\(cif_b[0-9]*_d[0-9]*\\)/\\1/p')
        [[ -z "\$suffix" ]] && suffix=rf3
        python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rosettafold3 "\$summary_conf" -o rfd3_${uid}_batch0_rf3_\${suffix}.tsv
    fi
    """
}

/*
Usage:

https://rosettacommons.github.io/foundry/models/rf3/index.html
https://github.com/RosettaCommons/foundry/tree/production/models/rf3

$ python models/rf3/src/rf3/inference.py --help

inference is powered by Hydra.

== Configuration groups ==
Compose your configuration from those groups (group=option)

callbacks: default, dump_validation_structures, metrics_logging, train_logging
dataloader: default
datasets: base, pdb_and_distillation, pdb_only
datasets/train: disorder_distillation, domain_distillation, monomer_distillation, na_complex_distillation, rna_monomer_distillation
datasets/train/pdb: af3_weighted_sampling, base, plinder, train_interface, train_pn_unit
datasets/val: af3_ab_set, af3_validation, base, runs_and_poses
debug: default, train_specific_examples
experiment: quick-rf3, quick-rf3-with-confidence
experiment/pretrained: rf3, rf3_with_confidence
inference_engine: base, rf3
logger: csv, default, wandb
model: rf3, rf3_with_confidence
model/components: ema, rf3_net, rf3_net_with_confidence_head
model/optimizers: adam
model/schedulers: af3
paths: default
paths/data: default
trainer: cpu, ddp, rf3, rf3_with_confidence, xpu
trainer/loss: structure_prediction, structure_prediction_with_confidence
trainer/loss/losses: confidence_loss, diffusion_loss, distogram_loss
trainer/metrics: structure_prediction


== Config ==
Override anything in the config (foo.bar=value)

ckpt_path: rf3_foundry_01_24_latest.ckpt
num_nodes: 1
devices_per_node: 1
compress_outputs: false
inputs: ???
out_dir: ???
dump_predictions: true
dump_trajectories: false
one_model_per_file: false
annotate_b_factor_with_plddt: true
sharding_pattern: null
skip_existing: false
template_selection: null
ground_truth_conformer_selection: null
cyclic_chains: []
_target_: rf3.inference_engines.rf3.RF3InferenceEngine
n_recycles: 10
diffusion_batch_size: 5
num_steps: 50
template_noise_scale: 1.0e-05
early_stopping_plddt_threshold: 0.5
seed: null
verbose: false
raise_if_missing_msa_for_protein_of_length_n: null
metrics_cfg:
  ptm:
    _target_: rf3.metrics.predicted_error.ComputePTM
  iptm:
    _target_: rf3.metrics.predicted_error.ComputeIPTM
  count_clashing_chains:
    _target_: rf3.metrics.clashing_chains.CountClashingChains


Powered by Hydra (https://hydra.cc)
Use --hydra-help to view Hydra specific help
*/