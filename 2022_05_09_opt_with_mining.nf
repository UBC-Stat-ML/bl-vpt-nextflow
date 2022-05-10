#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '134d9821b708441e270bfcb4df1b3185b0ea5234' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)
nChains = 8
nOptIters = 100
nScansPerGradient = 20

nCPUs = 5

params.dryRun = false

if (params.dryRun) {
  nOptIters = 5
  seeds = seeds.subList(0, 1)
  nCPUs = 2
}

nScans = nOptIters * nScansPerGradient
maxBudget = nScans * nChains

model_name = 'ptbm.models.Mining'

model_match =  ' --model ' + model_name + '\\\$Builder '
model_match += ' --model.counts file data/mining-disasters.csv  '



model_opt = model_match.replace('--model', '--model.interpolation.target')

fixed_scedule_match  = ' --engine.adaptFraction 0.0 '
fixed_scedule_match += ' --engine.ladder FromAnotherExec '
fixed_scedule_match += ' --engine.ladder.annealingParameters data/schedules/' + model_name + '/annealingParameters.csv '

fixed_scedule_opt = fixed_scedule_match.replace('--engine.', '--engine.pt.')


process run {

  time '5m'
  cpus nCPUs
  errorStrategy 'ignore'
  
  input:
    each seed from seeds
    each obj from 'SKL', 'FKL'
    each opt from 'Adam'
    each stepScale from 0.1
    file code
    file data
    
  output:
    file 'output' into results
    
  """
  java -Xmx5g -cp code/lib/\\*  ptgrad.Variational \
    --treatNaNAsNegativeInfinity true \
    --experimentConfigs.resultsHTMLPage false \
    --engine.maxOptimizationRestarts 1 \
    --model.interpolation Automatic \
    --engine ptgrad.VariationalPT \
    --engine.detailedGradientInfo false \
    --engine.pt.nScans 10 \
    $fixed_scedule_opt \
    --engine.nScansPerGradient $nScansPerGradient \
    --engine.optimizer.progressCheckLag 0.0 \
    --engine.pt.scmInit.nParticles 10 \
    --engine.pt.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.pt.nPassesPerScan 1 \
    $model_opt \
    --engine.antithetics IS \
    --engine.pt.nChains $nChains \
    --engine.pt.nThreads fixed \
    --engine.pt.nThreads.number $nCPUs \
    --engine.pt.random $seed \
    --engine.objective $obj \
    --engine.optimizer $opt \
    --engine.optimizer.maxIters $nOptIters \
    --engine.optimizer.stepScale $stepScale 
  mkdir output
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/*.tsv output
  """
}


process runMatching {

  time '5m'
  cpus nCPUs
  errorStrategy 'ignore'
  
  input:
    each seed from seeds
    file code
    file data
    
  output:
    file 'output' into results_matching
    
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --treatNaNAsNegativeInfinity true \
    --experimentConfigs.resultsHTMLPage false \
    --engine ptbm.OptPT \
    --engine.nScans $nScans \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    --engine.useFixedRefPT false \
    $model_match \
    --engine.nChains $nChains \
    --engine.nThreads fixed \
    --engine.nThreads.number $nCPUs \
    --engine.random $seed \
    --engine.minSamplesForVariational 10
  mkdir output
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/*.tsv output
  echo "\nengine.pt.random\t$seed" >> output/arguments.tsv
  echo "engine.objective\tFKL" >> output/arguments.tsv
  echo "engine.optimizer\tMomentMatch" >> output/arguments.tsv
  """
}

results_all = results.concat(results_matching)

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
    --dataPathInEachExecFolder optimizationMonitoring.csv optimizationPath.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      model.interpolation.target as model \
      engine.optimizer as optimizer \
      engine.antithetics as antithetics \
      engine.objective as objective \
      engine.optimizer.stepScale as stepScale \
      engine.pt.random as random \
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
     
  paths <- read.csv("${aggregated}/optimizationPath.csv.gz")
  paths <- filter(paths, budget <= $maxBudget)
  ggplot(paths, aes(x = budget, y = value, color = factor(random))) +
    facet_grid(objective + optimizer + name ~ factor(stepScale), labeller = label_both) +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("Parameter") + 
    ylim(-100, 350) +
    geom_line(alpha = 0.5)  + 
    theme_bw()
  ggsave(paste0("optimizationPaths.pdf"), width = 17, height = 30)
  
  optmonitor <- read.csv("${aggregated}/optimizationMonitoring.csv.gz")
  optmonitor <- filter(optmonitor, name == "Rejection")
  optmonitor <- filter(optmonitor, budget <= $maxBudget) # when hitting NaN, budget can be 2x larger
  ggplot(optmonitor, aes(x = budget, y = value, color = factor(random))) +
    facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("Global Communication Barrier (GCB)") + 
    geom_line(alpha = 0.5)  + 
    theme_bw()
  ggsave(paste0("optimizationMonitoring.pdf"), width = 17, height = 7)
  
  optmonitor %>% 
    filter(is.finite(value)) %>% 
    group_by(budget, objective, optimizer, stepScale) %>%
    summarise(mean_GCB = mean(value)) %>%
    ggplot(aes(x = budget, y = mean_GCB, colour = optimizer)) +
      facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
      scale_x_log10() +
      xlab("Budget (number of exploration steps)") + 
      ylab("GCB (averaged over 10 restarts, ignoring failures)") + 
      geom_line(alpha = 1) + 
      theme_bw()
  ggsave(paste0("optimizationMonitoring-mean.pdf"), width = 17, height = 7)
  
  optmonitor\$isFinite <- is.finite(optmonitor\$value)
  optmonitor %>% 
    group_by(budget, objective, optimizer, stepScale) %>%
    summarise(mean_is_finite = sum(isFinite)/${seeds.size()}) %>%
    ggplot(aes(x = budget, y = mean_is_finite, colour = optimizer)) +
      facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
      scale_x_log10() +
      ylim(0.0, 1.0) + 
      xlab("Budget (number of exploration steps)") + 
      ylab("Fraction of runs with finite objective") + 
      geom_line(alpha = 1) + 
      theme_bw()
  ggsave(paste0("optimizationMonitoring-isFinite.pdf"), width = 17, height = 7)
  

  """
  
}
