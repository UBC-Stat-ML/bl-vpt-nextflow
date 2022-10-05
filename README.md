This repository contains the nextflow scripts accompanying the paper
"Parallel Tempering With a Variational Reference" by Surjanovic, Syed, Bouchard-Côté, and Campbell. NeurIPS (2022).
[https://arxiv.org/abs/2206.00080](https://arxiv.org/abs/2206.00080)

### Required software

- bash
- the following should be available in the PATH variable:
    - Java Open SDK, version **11** (use [sdkman](https://sdkman.io/))
    - Rscript
- the following R packages should be installed:
    - ggplot2
    - dplyr
    - ggridges

### Replicating the optimization (section 4.2) and PT topology (4.4) results

Results will be written in a subdirectory of ``ptgrad-nextflow/deliverables``.

For the PT topology comparison, use ``./nextflow run pt_topologies.nf -resume``

Due to the various factorial designs, grid searches and replicates, these experiments are best performed on a cluster. For example, if your cluster uses a PBSPRO scheduler and supports singularity and conda, on the login node, use 

```
./nextflow-pbspro.sh run optimization.nf -resume
```

which will submit for you to the queue parallel runs via ``qsub``, gather the result and create the plots.

For other types of cluster schedulers, edit ``pbspro.config`` based on [the nextflow documentation](https://www.nextflow.io/docs/latest/executor.html). 


