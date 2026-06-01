# nf-binder-design

[![DOI](https://zenodo.org/badge/921968505.svg)](https://doi.org/10.5281/zenodo.16809704)

Nextflow pipelines for de novo protein binder design.

## Overview

This project provides Nextflow workflows for _de novo_ design of protein binders:

- **RFdiffusion → ProteinMPNN → AlphaFold2(initial guess) → Boltz-2 refolding** pipeline
- **RFdiffusion Partial Diffusion → Boltz-2 refolding** for diversification and optimization
- **RFdiffusion3** — [RosettaCommons/foundry](https://rosettacommons.github.io/foundry/models/rfd3/protein_binder_design.html) binder workflow: `rfd3` (RFdiffusion3) → `mpnn` → `rf3` (RosettaFold3) (→ Boltz-2 refolding)
- **BindCraft** - parallel execution across multiple GPUs
- **Germinal** - antibody/nanobody design in parallel across multiple GPUs
- **BoltzGen** - design proteins and complexes using BoltzGen
- **Boltz Pulldown** - an AlphaPulldown-like protocol using Boltz-2

![RFdiffusion workflow](images/rfd-workflow.png)

![BindCraft workflow](images/bindcraft-workflow.png)

## Features

- **Flexible workflow options** for different binder design strategies
- **HPC-ready** with support for SLURM and other job schedulers
- **Multi-GPU parallelization** for BindCraft trajectories
- **Plugin system** for custom design filters
- **Comprehensive reporting** with HTML outputs and summary statistics

## Quick Links

- [Setup Instructions](setup.md)
- [RFdiffusion Workflows](workflows/rfdiffusion.md)
- [RFdiffusion3 Workflow](workflows/rfd3.md)
- [BindCraft Workflow](workflows/bindcraft.md)
- [Germinal Workflow](workflows/germinal.md)
- [BoltzGen Workflow](workflows/boltzgen.md)
- [Boltz Pulldown](workflows/boltz-pulldown.md)
- [FoldSeek](subworkflows/foldseek.md)
- [Testing](testing.md)
- [Development](extra/development.md)
- [M3 HPC Examples](extra/m3-hpc-examples.md)
- [Utility Scripts](extra/utility-scripts.md)
- [Examples](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/examples/README.md)
- [GitHub Repository](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design)
- [Release history](changelog.md)

## License

The pipeline code the comprises `nf-binder-design` is licensed under the MIT License.

**Note** that some dependencies, packaged externally as containers, are under less permissive licenses:

> ⚠️ "Commercial Use Restrictions"
> >
> Components of these workflows use RFdiffusion and BindCraft, which depend on PyRosetta/Rosetta. These are  **free for non-commercial use only**. Commercial use requires a paid license agreement with University of Washington. See:
>
> - [Rosetta License](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md)
> - [Rosetta Licensing FAQ](https://rosettacommons.org/software/licensing-faq/)

## Citing

If you use `nf-binder-design` in your research, please cite:

- Perry, A., Taveneau, C., & Knott, G. J. (2025). nf-binder-design: a Nextflow pipeline for protein binder design (`<include version used, eg 0.2.0>`). Zenodo. [https://doi.org/10.5281/zenodo.16809705](https://doi.org/10.5281/zenodo.16809705)

and include citations for the underlying tools used in the workflow as appropriate:
  
  - RFdiffusion
    - Watson, J.L., Juergens, D., Bennett, N.R. et al. "De novo design of protein structure and function with RFdiffusion.", _Nature_, **620**, 1089–1100 (2023). [https://doi.org/10.1038/s41586-023-06415-8](https://doi.org/10.1038/s41586-023-06415-8)
  
    - Bennett, N.R., Coventry, B., Goreshnik, I. et al. Improving de novo protein binder design with deep learning. _Nat Commun_, **14**, 2625 (2023). [https://doi.org/10.1038/s41467-023-38328-5](https://doi.org/10.1038/s41467-023-38328-5)
  
  - RFdiffusion3
    
    - Butcher, J., Krishna, R., Mitra, R. et al. De novo Design of All-atom Biomolecular Interactions with RFdiffusion3. _bioRxiv_ (2025). [https://doi.org/10.1101/2025.09.18.676967](https://doi.org/10.1101/2025.09.18.676967)

  - ProteinMPNN / LigandMPNN / SolubleMPNN
  
    - Dauparas, J. et al. Robust deep learning–based protein sequence design using ProteinMPNN. _Science_, **378**,49-56(2022). [https://doi.org/10.1126/science.add2187](https://doi.org/10.1126/science.add2187)

    - Dauparas, J. et al. Atomic context-conditioned protein sequence design using LigandMPNN. _Nature Methods_, **22**, 717–723 (2025). [https://doi.org/10.1038/s41592-025-02626-1](https://doi.org/10.1038/s41592-025-02626-1)

    - Goverde, C. A. et al. Computational design of soluble and functional membrane protein analogues. _Nature_, **631**, 449–458 (2024). [https://doi.org/10.1038/s41586-024-07601-y](https://doi.org/10.1038/s41586-024-07601-y)

  - RosettaFold3
    
    - Corley, N. et al. Accelerating Biomolecular Modeling with AtomWorks and RF3. _bioRxiv_ (2025). [https://doi.org/10.1101/2025.08.14.670328](https://doi.org/10.1101/2025.08.14.670328)

  - Alphafold2
     - Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). [https://doi.org/10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
   
   - BindCraft
     - Pacesa, M., Nickel, L., Schellhaas, C. et al. One-shot design of functional protein binders with BindCraft. _Nature_, **645**, 1005-1010 (2025). [https://doi.org/10.1038/s41586-025-09429-6](https://doi.org/10.1038/s41586-025-09429-6)

  - Boltz
    - Passaro, S., Corso, G., Wohlwend, J. et al. Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. _bioRxiv_ (2025). [https://doi.org/10.1101/2025.06.14.659707](https://doi.org/10.1101/2025.06.14.659707)
    - Wohlwend, J., Corso, G., Passaro, S. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. _bioRxiv_ (2024). [https://doi.org/10.1101/2024.11.19.624167](https://doi.org/10.1101/2024.11.19.624167)
  
  - ColabFold
    - Mirdita, M., Schütze, K., Moriwaki, Y. et al. ColabFold: making protein folding accessible to all. _Nature Methods_, **19**, 679-682 (2022). [https://doi.org/10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)

  - BoltzGen
    - Stark, H., et al. "BoltzGen: Toward Universal Binder Design." Preprint (2025). [https://hannes-stark.com/assets/boltzgen.pdf](https://hannes-stark.com/assets/boltzgen.pdf) (accessed November 10, 2025).
