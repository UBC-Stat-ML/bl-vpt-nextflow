# ptgrad-nextflow
This repository contains the Blang PPL scripts accompanying the paper
"Parallel Tempering With a Variational Reference" by Surjanovic, Syed, Bouchard-Côté, and Campbell (2022).
The repository https://github.com/UBC-Stat-ML/ptanalysis contains Blang PPL code.

### Required software

- bash
- the following should be available in the PATH variable:
    - Java Open SDK, version **11** (use [sdkman](https://sdkman.io/))
    - Rscript
- the following R packages should be installed:
    - ggplot2
    - dplyr
    - ggridges


### Compiling

```
cd ptanalysis
./gradlew installDist
```


### Running the stabilized moment matching on one model

Still from the ``ptanalysis`` folder:

```
chmod 755 pt-matching.sh
./pt-matching.sh --model ptbm.models.ToyMix
```

To see available options, use ``./pt-matching.sh --model ptbm.models.ToyMix --help``

To find other existing available models, use ``ls src/main/java/ptbm/models/``


### Replicating the optimization (section 4.2) and PT topology (4.4) results

Still from the ``ptanalysis`` folder:

```
./gradlew clean
cd ../ptgrad-nextflow/
./nextflow run optimization.nf -resume
```

Results will be written in a subdirectory of ``ptgrad-nextflow/deliverables``.

For the PT topology comparison, use ``./nextflow run pt_topologies.nf -resume``

Due to the various factorial designs, grid searches and replicates, these experiments are best performed on a cluster. For example, if your cluster uses a PBSPRO scheduler and supports singularity and conda, on the login node, use 

```
./nextflow-pbspro.sh run optimization.nf -resume
```

which will submit for you to the queue parallel runs via ``qsub``, gather the result and create the plots.

For other types of cluster schedulers, edit ``pbspro.config`` based on [the nextflow documentation](https://www.nextflow.io/docs/latest/executor.html). 


### Adding a new model

- You can base your new model on one of the models in ``ptanalysis/src/main/java/ptbm/models/``
- Implement the model based off the [Blang documentation](https://arxiv.org/pdf/1912.10396.pdf) with the following amendments:
    - For each latent random variable that you would like to approach variationally, enclose the distribution declaration with ``Opt``, and pass in the target distribution as argument as in the following example: ``mu ~ Normal(0.0, 1.0)`` becomes ``mu ~ Opt(Normal::distribution(0.0, 1.0))``
    - Declare the type of such random variable as ``VariationalReal`` instead of ``RealVar``.
- Then follow the compilation instructions above.
