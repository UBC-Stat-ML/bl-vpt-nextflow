#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '3daaeda130d8afec4b456f4b49d172022c6c1ee6' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}

seeds = (1..5)

models = []

class Model {
  String name
  String stat
  String args
  def int hashCode() { return [name, stat, args].hashCode() }
}

def addModel(String n, String s, String a) {
  m = new Model(name: n, stat: s, args: a)
  models.add(m)
}

nScans = 10000
nScans_ref = 100000

addModel('mrna-no-transf', 'beta',  ' --model ptbm.models.MRNATransfectionNoTransform --model.data data/m_rna_transfection/processed.csv --engine.nChains 30 ')
addModel('toy-mix',        'x',     ' --model ptbm.models.ToyMix --engine.nChains 10 ')
addModel('sparse-car',     'alpha', ' --model ptbm.models.SparseCAR --model.data data/scotland_lip_cancer/data.csv --model.spatialData.adjacency data/scotland_lip_cancer/adj.csv --engine.nChains 20 ')
addModel('coll-rockets',   'p0',    ' --model ptbm.models.CollapsedHierarchicalRockets --model.data data/failure_counts.csv       --engine.nChains 10 ')
addModel('8-schools',      'tau',   ' --model ptbm.models.EightSchools --model.data data/eight-schools.csv --model.useTPrior true --engine.nChains 10 ')
addModel('titanic',        'sigma', ' --model ptbm.models.LogisticRegression --model.data data/titanic/titanic-covariates-original.csv --model.instances.name Name --model.instances.maxSize 20 --model.labels.dataSource data/titanic/titanic.csv --model.labels.name Survived--model.useTPrior false --engine.nChains 10 ')
addModel('mining',         's',     ' --model ptbm.models.Mining --model.counts file data/mining-disasters.csv --engine.nChains 10 ')
addModel('phylo',          'rate',  ' --model ptbm.models.PhylogeneticTree --model.observations.file data/FES_8.g.fasta --model.observations.encoding DNA --engine.nChains 20 ')
addModel('vaccines',       'eff_m', ' --model ptbm.models.Vaccines --model.data data/vaccines/data.csv --engine.nChains 20 ')


algos = [:]  
reference = 'Reference'
algos['F--T  F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational INF '
algos[reference]     = algos['F--T  F--T'] 
algos['V--T* F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational 100 '
algos['V--T* F--T*'] = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational false '
algos['V--T*--F']    = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational true '
algos['F--T--F']     = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational INF --engine.doSwapFixedRefAndVariational true '

postprocessor = '' //' --postProcessor ptgrad.VariationalPostprocessor '

params.dryRun = false

if (params.dryRun) {
  nScans = 200
  nScans_ref = 400
  seeds = (1..2)
}

process runMatching {

  time '5h'
  errorStrategy 'ignore'
  
  input:
    each model from models
    each seed from seeds
    each algo from algos.entrySet()
    file code
    file data
    
  output:
    file 'output' into results
    file 'fixedRefOutput' into fixedRefResults
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    $postprocessor \
    --engine ptbm.OptPT \
    --engine.random $seed \
    --engine.nScans ${if (algo.key == reference) nScans_ref else nScans} \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    --engine.useFixedRefPT true \
    ${algo.value} \
    ${model.args} \
    --engine.nThreads single
  mkdir output
  mkdir fixedRefOutput
  
  mv results/latest/samples/${model.stat}.csv output/statistic.csv
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/fixedReferencePT/monitoring/*.csv fixedRefOutput
  mv results/latest/fixedReferencePT/samples/${model.stat}.csv fixedRefOutput/statistic.csv
  
  echo "\nmodelDescription\t${model.name}" >> results/latest/arguments.tsv
  echo "algo\t${algo.key}" >> results/latest/arguments.tsv
  echo "path\t\$(pwd)" >> results/latest/arguments.tsv
  echo "statisticName\t${model.stat}" >> results/latest/arguments.tsv
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
    --dataPathInEachExecFolder lambdaInstantaneous.csv actualTemperedRestarts.csv globalLambda.csv logNormalizationConstantProgress.csv statistic.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      modelDescription as model \
      engine.random as seed \
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
  require("ggridges")
  
  stat <- read.csv("${aggregated}/statistic.csv.gz")
  
  ks_distances <- data.frame(
    algo = factor(),
    model = factor(),
    seed = integer(),
    ks_p = double(),
    ks_stat = double()
  )

  for (a in unique(stat\$algo)) {
      for (m in unique(stat\$model)) {
        for (s in unique(stat\$seed)) {
          ref <- stat %>%
            filter(sample > ${nScans_ref/2}) %>%
            filter(algo == "Reference") %>%
            filter(fixedRefChain == "false") %>%
            filter(model == m) %>%
            filter(seed == s) %>%
            .\$value
          cur <- stat %>%
            filter(sample > ${nScans/2}) %>%
            filter(algo == a) %>%
            filter(fixedRefChain == "false") %>%
            filter(model == m) %>%
            filter(seed == s) %>%
            .\$value
          test <- ks.test(ref, cur)
          ks_distances <- ks_distances %>% add_row(algo = a, model = m, seed = s, ks_p = test\$p.value, ks_stat = as.numeric(test\$statistic))
        }
      }
  }
  
  ks_distances %>%
    ggplot(aes(x = seed, y = ks_stat)) +
      facet_grid(algo ~ model) +
      geom_point() +
      theme_bw()
  ggsave(paste0("ks.pdf"), width = 17, height = 8)
  
  stat <- stat %>% inner_join(ks_distances, by = c("algo", "model", "seed"))
  stat\$large_disc <- ifelse(stat\$ks_stat > 0.2, 1, 0)
  
  
  stat %>%
    filter(sample > ${nScans/2}) %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = value, y = factor(seed), fill = factor(large_disc))) +
      facet_grid(algo ~ model, scales="free") +
      geom_density_ridges(panel_scaling = TRUE) +
      theme_bw()
  ggsave(paste0("statistic.pdf"), width = 17, height = 8) 
 
  restarts <- read.csv("${aggregated}/actualTemperedRestarts.csv.gz")
  restarts %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = round, y = count, color = factor(seed), group = factor(seed))) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      theme_bw()
  ggsave(paste0("actualTemperedRestarts.pdf"), width = 17, height = 8)
  
  global <- read.csv("${aggregated}/globalLambda.csv.gz")
  global %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = round, y = value, color = factor(seed), group = factor(seed))) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      theme_bw()
  ggsave(paste0("globalLambda.pdf"), width = 17, height = 8)
  

  
  """
  
}

/*
  lambda <- read.csv("${aggregated}/lambdaInstantaneous.csv.gz")
  lambda %>% 
    filter(fixedRefChain == "false") %>%
    filter(isAdapt == "false") %>%
    mutate(beta2 = ifelse(fixedRefChain == "true", -beta, beta)) %>%
    ggplot(aes(x = beta2, y = value, color = factor(seed), group = factor(seed)) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      xlab("Beta") + 
      ylab("Intensity") + 
      theme_bw()
  ggsave(paste0("lambda.pdf"), width = 17, height = 8)
  
  cnst <- read.csv("${aggregated}/logNormalizationConstantProgress.csv.gz")
  cnst %>% 
    filter(fixedRefChain == "false") %>%
    filter(round > 8) %>%
    ggplot(aes(x = round, y = value, color = factor(seed), group = factor(seed))) +
      facet_grid(model ~ algo, scales="free_y") +
      geom_line() +
      theme_bw()
  ggsave(paste0("logNormalizationConstantProgress.pdf"), width = 17, height = 8)

*/

