#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '00243a456796bc1264b43bcbecff057110f00a5e' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}



models = [:]
models['coll-rockets']   = ' --model ptbm.models.CollapsedHierarchicalRockets\\\$Builder --model.data data/failure_counts.csv       --engine.nChains 10 --engine.nScans 5_000 '
models['8-schools']      = ' --model ptbm.models.EightSchools\\\$Builder --model.data data/eight-schools.csv --model.useTPrior true --engine.nChains 10 --engine.nScans 5_000 '
models['titanic']      = ' --model ptbm.models.LogisticRegression\\\$Builder --model.data data/titanic/titanic-covariates-original.csv --model.instances.name Name --model.instances.maxSize 20 --model.labels.dataSource data/titanic/titanic.csv --model.labels.name Survived--model.useTPrior false --engine.nChains 10 --engine.nScans 5_000 '
models['mining']         = ' --model ptbm.models.Mining\\\$Builder --model.counts file data/mining-disasters.csv --engine.nChains 10 --engine.nScans 5_000 '
models['mrna-no-transf'] = ' --model ptbm.models.MRNATransfectionNoTransform\\\$Builder --model.data data/m_rna_transfection/processed.csv --engine.nChains 10 --engine.nScans 5_000 '
models['phylo']          = ' --model ptbm.models.PhylogeneticTree\\\$Builder --model.observations.file data/FES_8.g.fasta --model.observations.encoding DNA --engine.nChains 20 --engine.nScans 10_000 '
models['single-cell']    = ' --model ptbm.models.SingleCell\\\$Builder --model.data.source data/chromo/gc.csv --model.data.gcContents.name value --model.data.readCounts.name value --model.data.readCounts.dataSource data/chromo/7.csv --model.configs.annealingStrategy Exponentiation --model.configs.annealingStrategy.thinning 1'
models['sparse-car']     = ' --model ptbm.models.SparseCAR\\\$Builder --model.data data/scotland_lip_cancer/data.csv --model.spatialData.adjacency data/scotland_lip_cancer/adj.csv --engine.nChains 20 --engine.nScans 10_000'
models['toy-mix']        = ' --model ptbm.models.ToyMix\\\$Builder --engine.nChains 10 --engine.nScans 5_000 '
models['vaccines']       = ' --model ptbm.models.Vaccines\\\$Builder --model.data data/vaccines/data.csv --engine.nChains 20 --engine.nScans 10_000 '


algos = [:]  
algos['F--T  F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational INF '
algos['V--T* F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational 100 '
algos['V--T* F--T*'] = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational false '
algos['V--T*--F']    = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational true '
algos['F--T--F']     = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational INF --engine.doSwapFixedRefAndVariational true '


params.dryRun = false

if (params.dryRun) {
}

process runMatching {

  time '1h'
  errorStrategy 'ignore'
  
  input:
    each model from models.entrySet()
    each algo from algos.entrySet()
    file code
    file data
    
  output:
    file 'output' into results
    file 'fixedRefOutput' into fixedRefResults
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    --engine ptbm.OptPT \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    --engine.useFixedRefPT true \
    ${algo.value} \
    ${model.value} \
    --engine.nThreads single
  mkdir output
  mkdir fixedRefOutput
  
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/fixedReferencePT/monitoring/*.csv fixedRefOutput
  
  echo "\nmodelDescription\t${model.key}" >> results/latest/arguments.tsv
  echo "algo\t${algo.key}" >> results/latest/arguments.tsv
  cp results/latest/*.tsv output
  cp results/latest/*.tsv fixedRefOutput
  echo "fixedRefOutput\tfalse" >> output/arguments.tsv
  echo "fixedRefOutput\ttrue"  >> fixedRefOutput/arguments.tsv
  
  """ 

}

results_all = results.concat(fixedRefResults)




process aggregate {
  time '1h'
  echo false
  scratch false
  input:
    file 'exec_*' from results_all.toList()
  output:
    file 'results/aggregated/' into aggregated
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder lambdaInstantaneous.csv actualTemperedRestarts.csv globalLambda.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      modelDescription as model \
      algo \
      fixedRefOutput as fixedRefChain \
           from arguments.tsv
  mv results/latest results/aggregated
  """
}


process plot {
  scratch false
  input:
    file aggregated
  output:
    file '*.*'
    file 'aggregated' 
  afterScript 'rm Rplots.pdf; cp .command.sh rerun.sh'
  container 'cgrlab/tidyverse'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr") 
     
  lambda <- read.csv("${aggregated}/lambdaInstantaneous.csv.gz")
  lambda %>% 
    filter(isAdapt == "false") %>%
    mutate(beta2 = ifelse(fixedRefChain == "true", -beta, beta)) %>%
    ggplot(aes(x = beta2, y = value, color = fixedRefChain)) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      xlab("Beta") + 
      ylab("Intensity") + 
      theme_bw()
  ggsave(paste0("lambda.pdf"), width = 17, height = 8)
  
  restarts <- read.csv("${aggregated}/actualTemperedRestarts.csv.gz")
  restarts %>%
    ggplot(aes(x = round, y = count, color = fixedRefChain)) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      theme_bw()
  ggsave(paste0("actualTemperedRestarts.pdf"), width = 17, height = 8)
  
  global <- read.csv("${aggregated}/globalLambda.csv.gz")
  global %>%
    ggplot(aes(x = round, y = value, color = fixedRefChain)) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      theme_bw()
  ggsave(paste0("globalLambda.pdf"), width = 17, height = 8)
  
  """
  
}