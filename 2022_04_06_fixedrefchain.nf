#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')
data = file("data")

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'bfb0467a8d53863e189ad8031b8a542d449d5763'
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)
Ns = (2..3).collect{Math.pow(2, it)}
nOptIters = 1000

model_match =  ' --model ptbm.models.CollapsedHierarchicalRockets\\\$Builder '
model_match += ' --model.data data/failure_counts.csv '

model_opt = model_match.replace('--model', '--model.interpolation.target')


process runMatching {

  time '2h'
  //errorStrategy 'ignore'
  
  cpu 4
  
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
    --engine.nScans 25000 \
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


process analysisCode {
  executor 'local'
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'a9abcc40abcfb285588cc4c312d8ecc0bbdad06e'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}


process aggregate {
  time '1h'
  echo false
  scratch false
  input:
    file analysisCode
    file 'exec_*' from results.toList()
  output:
    file 'results/latest/' into aggregated
  """
  java -Xmx5g -cp code/lib/\\*  flows.Aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder optimizationMonitoring.csv optimizationPath.csv \
    --keys \
      model.interpolation.target as model \
      engine.optimizer as optimizer \
      engine.antithetics as antithetics \
      engine.objective as objective \
      engine.useFixedRefPT as useFixedRef \
      engine.optimizer.stepScale as stepScale \
      engine.pt.nChains as nChains \
      engine.pt.random as random \
           from arguments.tsv
  """
}

process plot {
  scratch false
  input:
    file aggregated
  output:
    file '*.pdf'
  //  file '*.csv'
  afterScript 'rm Rplots.pdf'
  container 'cgrlab/tidyverse'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  
  
  paths <- read.csv("${aggregated}/optimizationPath.csv")
  ggplot(paths, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
    facet_grid(name ~ nChains) +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("Parameter") + 
    geom_line(alpha = 0.5)  + 
    theme_bw()
  ggsave(paste0("optimizationPaths.pdf"), width = 9, height = 10)
  
  optmonitor <- read.csv("${aggregated}/optimizationMonitoring.csv")
  optmonitor <- filter(optmonitor, name == "Rejection")
  ggplot(optmonitor, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
    scale_x_log10() +
    facet_grid(. ~ nChains) +
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
      facet_grid(. ~ nChains) +
      xlab("Budget (number of exploration steps)") + 
      ylab("GCB (averaged over 10 restarts, ignoring failures)") + 
      geom_line(alpha = 1) + 
      theme_bw()
  ggsave(paste0("optimizationMonitoring-mean.pdf"), width = 9, height = 7)
  


  """
  
}