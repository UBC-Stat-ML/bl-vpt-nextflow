#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')
data = file("data")

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'c47993a95342163c345147b5ffbd0ea91a933e07'
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)
Ns = (4..4)  //.collect{Math.pow(2, it)}
nOptIters = 1000

model_match =  ' --model ptbm.models.CollapsedHierarchicalRockets\\\$Builder '
model_match += ' --model.data data/failure_counts.csv '
//model_match += ' --model.filter Ariane '

model_opt = model_match.replace('--model', '--model.interpolation.target')


process run {

  time '1h'
  errorStrategy 'ignore'
  
  input:
    each seed from seeds
    each obj from 'SKL', 'Rejection'
    each opt from 'Adam', 'SGD --engine.optimizer.schedule.exponent -0.5 '
    each stepScale from 0.01, 0.1, 1.0, 10.0
    each nChain from Ns
    file code
    file data
    
  output:
    file 'output' into results
    
  """
  java -Xmx5g -cp code/lib/\\*  ptgrad.Variational \
    --experimentConfigs.resultsHTMLPage false \
    --engine.maxOptimizationRestarts 1 \
    --model.interpolation Automatic \
    --engine ptgrad.VariationalPT \
    --engine.detailedGradientInfo false \
    --engine.pt.nScans 100 \
    --engine.nScansPerGradient 20 \
    --engine.optimizer.progressCheckLag 0.0 \
    --engine.pt.scmInit.nParticles 10 \
    --engine.pt.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.pt.nPassesPerScan 1 \
    $model_opt \
    --engine.antithetics IS \
    --engine.pt.nChains $nChain \
    --engine.pt.nThreads single \
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

  time '1h'
  //errorStrategy 'ignore'
  
  input:
    each seed from seeds
    each nChain from Ns
    file code
    file data
    
  output:
    file 'output' into results_matching
    
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
    --engine.nThreads single \
    --engine.random $seed \
    --engine.minSamplesForVariational 10
  mkdir output
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/*.tsv output
  echo "\nengine.pt.random\t$seed" >> output/arguments.tsv
  echo "engine.objective\tForwardKL" >> output/arguments.tsv
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

results_all = results.concat(results_matching)

process aggregate {
  time '1h'
  echo false
  scratch false
  input:
    file analysisCode
    file 'exec_*' from results_all.toList()
  output:
    file 'results/latest/' into aggregated
  """
  java -Xmx5g -cp code/lib/\\*  flows.Aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder optimizationMonitoring.csv \
    --keys \
      model.interpolation.target as model \
      engine.optimizer as optimizer \
      engine.antithetics as antithetics \
      engine.objective as objective \
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
  
  optmonitor <- read.csv("${aggregated}/optimizationMonitoring.csv")
  optmonitor <- filter(optmonitor, name == "Rejection")
  optmonitor <- filter(optmonitor, budget <= 75000) # when hitting NaN, budget can be 2x larger
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