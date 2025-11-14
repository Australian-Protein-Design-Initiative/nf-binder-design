process BOLTZGEN_MERGE {
    tag "merge"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    publishDir path: "${params.outdir}/boltzgen", pattern: '**', mode: 'copy'

    input:
    path batch_dirs

    output:
    path 'merged', type: 'dir', emit: merged_dir

    script:
    """
    mkdir -p merged/intermediate_designs_inverse_folded
    
    # Find all batch directories (batch_*) - follow symlinks with -L
    # When collected, batch directories are staged as symlinks in the work directory
    find -L . -maxdepth 1 -type d -name "batch_*" | while read batch_dir; do
        batch_path="\${batch_dir}"
        
        # Merge intermediate_designs directory (from design step)
        if [ -d "\${batch_path}/intermediate_designs" ]; then
            mkdir -p merged/intermediate_designs
            # Copy all .cif and .npz files from intermediate_designs
            find "\${batch_path}/intermediate_designs" -maxdepth 1 \\( -name "*.cif" -o -name "*.npz" \\) -exec cp {} merged/intermediate_designs/ \\; 2>/dev/null || true
            
            # Merge molecules_out_dir (required for protein-small_molecule protocol)
            if [ -d "\${batch_path}/intermediate_designs/molecules_out_dir" ]; then
                mkdir -p merged/intermediate_designs/molecules_out_dir
                cp \${batch_path}/intermediate_designs/molecules_out_dir/*.pkl merged/intermediate_designs/molecules_out_dir/ 2>/dev/null || true
            fi
        fi
        
        # Merge inverse-folded design .cif and .npz files (required for analysis)
        # These are the main output files from inverse_folding + folding steps
        if [ -d "\${batch_path}/intermediate_designs_inverse_folded" ]; then
            # Copy all .cif files from the intermediate_designs_inverse_folded directory
            find "\${batch_path}/intermediate_designs_inverse_folded" -maxdepth 1 -name "*.cif" -exec cp {} merged/intermediate_designs_inverse_folded/ \\; 2>/dev/null || true
            # Copy all .npz files from the intermediate_designs_inverse_folded directory (required for analysis)
            find "\${batch_path}/intermediate_designs_inverse_folded" -maxdepth 1 -name "*.npz" -exec cp {} merged/intermediate_designs_inverse_folded/ \\; 2>/dev/null || true
            
            # Merge metrics_tmp directory (contains data_* and metrics_* files)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/metrics_tmp" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/metrics_tmp
                cp \${batch_path}/intermediate_designs_inverse_folded/metrics_tmp/*.npz merged/intermediate_designs_inverse_folded/metrics_tmp/ 2>/dev/null || true
            fi
            
            # Merge refold_cif directories (required for analysis)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/refold_cif" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/refold_cif
                cp \${batch_path}/intermediate_designs_inverse_folded/refold_cif/*.cif merged/intermediate_designs_inverse_folded/refold_cif/ 2>/dev/null || true
            fi
            
            # Merge refold_design_cif (if design_folding was run)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/refold_design_cif" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/refold_design_cif
                cp \${batch_path}/intermediate_designs_inverse_folded/refold_design_cif/*.cif merged/intermediate_designs_inverse_folded/refold_design_cif/ 2>/dev/null || true
            fi
            
            # Merge fold_out_npz (optional, for metrics)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/fold_out_npz" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/fold_out_npz
                cp \${batch_path}/intermediate_designs_inverse_folded/fold_out_npz/*.npz merged/intermediate_designs_inverse_folded/fold_out_npz/ 2>/dev/null || true
            fi
            
            # Merge fold_out_design_npz (if design_folding was run)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/fold_out_design_npz" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/fold_out_design_npz
                cp \${batch_path}/intermediate_designs_inverse_folded/fold_out_design_npz/*.npz merged/intermediate_designs_inverse_folded/fold_out_design_npz/ 2>/dev/null || true
            fi
            
            # Merge affinity_out_npz (if affinity step was run for protein-small_molecule protocol)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/affinity_out_npz" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/affinity_out_npz
                cp \${batch_path}/intermediate_designs_inverse_folded/affinity_out_npz/*.npz merged/intermediate_designs_inverse_folded/affinity_out_npz/ 2>/dev/null || true
            fi
            
            # Merge molecules_out_dir (required for protein-small_molecule protocol)
            if [ -d "\${batch_path}/intermediate_designs_inverse_folded/molecules_out_dir" ]; then
                mkdir -p merged/intermediate_designs_inverse_folded/molecules_out_dir
                cp \${batch_path}/intermediate_designs_inverse_folded/molecules_out_dir/*.pkl merged/intermediate_designs_inverse_folded/molecules_out_dir/ 2>/dev/null || true
            fi
        fi
        
        # Copy config directory and steps.yaml from first batch (they should be identical across batches)
        if [ ! -d "merged/config" ] && [ -d "\${batch_path}/config" ]; then
            cp -r \${batch_path}/config merged/ 2>/dev/null || true
        fi
        if [ ! -f "merged/steps.yaml" ] && [ -f "\${batch_path}/steps.yaml" ]; then
            cp \${batch_path}/steps.yaml merged/ 2>/dev/null || true
        fi
    done
    """
}
