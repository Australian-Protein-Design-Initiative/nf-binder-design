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
        
        # Merge inverse-folded design .cif and .npz files (required for analysis)
        # These are the main output files from inverse_folding + folding steps
        if [ -d "\${batch_path}/intermediate_designs_inverse_folded" ]; then
            # Copy all .cif files from the intermediate_designs_inverse_folded directory
            find "\${batch_path}/intermediate_designs_inverse_folded" -maxdepth 1 -name "*.cif" -exec cp {} merged/intermediate_designs_inverse_folded/ \\; 2>/dev/null || true
            # Copy all .npz files from the intermediate_designs_inverse_folded directory (required for analysis)
            find "\${batch_path}/intermediate_designs_inverse_folded" -maxdepth 1 -name "*.npz" -exec cp {} merged/intermediate_designs_inverse_folded/ \\; 2>/dev/null || true
            
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
        fi
    done
    """
}
