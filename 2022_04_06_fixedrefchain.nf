#!/usr/bin/env nextflow

// modified from optimization.nf

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')
data = file("data")

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '2dff8cda5cb0bd73b5f8f2fd842f1ac994165225'
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)
Ns = (2..2).collect{Math.pow(2, it)}
nOptIters = 1000

model_match =  ' --model ptbm.models.CollapsedHierarchicalRockets\\\$Builder '
model_match += ' --model.data data/failure_counts.csv '

model_opt = model_match.replace('--model', '--model.interpolation.target')


process runMatching {

  time '20m'
  //errorStrategy 'ignore'
  cpus 4
  
  input:
    each seed from seeds
    each nChain from Ns
    each useRef from 'true', 'false'
    file code
    file data
    
  output:
    file 'output' into results
    
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    --engine ptbm.OptPT \
    --engine.nScans 10000 \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    $model_match \
    --engine.nChains $nChain \
    --engine.useFixedRefPT $useRef \
    --engine.nThreads fixed \
    --engine.nThreads.number 4 \
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

process aggregate {
  time '5m'
  echo false
  scratch false
  input:
    file 'exec_*' from results.toList()
  output:
    file 'results/aggregated/' into aggregated
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --experimentConfigs.tabularWriter.compressed true \
    --dataPathInEachExecFolder optimizationMonitoring.csv optimizationPath.csv \
    --keys \
      model.interpolation.target as model \
      engine.optimizer as optimizer \
      engine.antithetics as antithetics \
      engine.objective as objective \
      engine.useFixedRefPT as useFixedRef \
      engine.optimizer.stepScale as stepScale \
      engine.nChains as nChains \
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
  ggplot(paths, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
    facet_grid(name ~ useFixedRef, scales="free_y") +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("Parameter") + 
    geom_line(alpha = 0.5)  + 
    theme_bw()
  ggsave(paste0("optimizationPaths.pdf"), width = 9, height = 10)
  
  optmonitor <- read.csv("${aggregated}/optimizationMonitoring.csv.gz")
  optmonitor <- filter(optmonitor, name == "Rejection")
  ggplot(optmonitor, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
    scale_x_log10() +
    facet_grid(. ~ useFixedRef) +
    xlab("Budget (number of exploration steps)") + 
    ylab("Global Communication Barrier (GCB)") + 
    geom_line(alpha = 0.5)  + 
    theme_bw()
  ggsave(paste0("optimizationMonitoring.pdf"), width = 9, height = 7)
  
  optmonitor %>% 
    filter(is.finite(value)) %>% 
    group_by(budget, objective, optimizer, stepScale, useFixedRef, nChains) %>%
    summarise(mean_GCB = mean(value)) %>%
    ggplot(aes(x = budget, y = mean_GCB, colour = optimizer, linetype = useFixedRef)) +
      scale_x_log10() +
      xlab("Budget (number of exploration steps)") + 
      ylab("GCB (averaged over 10 restarts, ignoring failures)") + 
      geom_line(alpha = 1) + 
      theme_bw()
  ggsave(paste0("optimizationMonitoring-mean.pdf"), width = 9, height = 7)
  """
}