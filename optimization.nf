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

seeds = (1..1)
Ns = (4..4).collect{Math.pow(2, it)}
nOptIters = 2

model_match =  ' --model ptbm.models.CollapsedHierarchicalRockets\\\$Builder '
model_match += ' --model.data data/failure_counts.csv '
model_match += ' --model.filter Ariane '

model_opt = model_match.replace('--model', '--model.interpolation.target')




process run {

  time '30m'
  //errorStrategy 'ignore'
  
  input:
    each seed from seeds
    each obj from 'SKL' //, 'Rejection'
    each opt from 'Adam' //, 'SGD --engine.optimizer.schedule.exponent -0.5 '
    each nChain from Ns
    file code
    file data
    
  output:
    file 'output' into results
    
  """
  java -Xmx5g -cp code/lib/\\*  ptgrad.Variational \
    --experimentConfigs.resultsHTMLPage false \
    --model.interpolation Automatic \
    --engine ptgrad.VariationalPT \
    --engine.detailedGradientInfo false \
    --engine.pt.nScans 100 \
    --engine.nScansPerGradient 100 \
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
    --engine.optimizer.stepScale 1.0 
  mkdir output
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/*.tsv output
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
  echo true
  scratch false
  input:
    file analysisCode
    file 'exec_*' from results.toList()
  output:
    file 'results/latest/' into aggregated
  """
  java -Xmx5g -cp code/lib/\\*  flows.Aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder optimizationMonitoring.csv \
    --keys \
      model.interpolation as model \
      engine.antithetics as antithetics \
      engine.objective \
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
    file '*.csv'
  container 'cgrlab/tidyverse'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  
  timings <- read.csv("${aggregated}/roundTimings.csv")
  ggplot(timings, aes(x = nChains, y = value/1000)) +
    ylab("seconds") + 
    facet_grid(model ~ .) +
    geom_point() + 
    theme_bw()
  ggsave(paste0("timings.pdf"), width = 5, height = 5)
  
  data <- read.csv("${aggregated}/stochastic-gradient-evaluations.csv") %>% filter(evaluation == 0)
   
  byChain <- data %>% 
    group_by(model,chain,nChains,objectiveType,antithetics,dim,essn) %>%
    summarize(
      mean = mean(value),
      sd = sd(value),
      snr = abs(mean)/sd)
  write.csv(byChain, "chain-specific-summary.csv")
  
  for (stat in c("snr", "sd", "mean")) {  
    ggplot(byChain, aes(x = chain, y = get(stat), linetype = objectiveType, color = interaction(antithetics,essn))) +
      facet_grid(model + dim ~ nChains, scales = "free_y") +
      geom_line() + 
      theme_bw()
    ggsave(paste0("chain-specific-", stat, ".pdf"), width = 15, height = 15)
  }
  
  overall <- data %>% 
    group_by(model,nChains,objectiveType,antithetics,dim,essn,random) %>%
    mutate(sum = sum(value)) %>%
    group_by(model,nChains,objectiveType,antithetics,dim,essn) %>%
    summarize(
      mean = mean(sum),
      sd = sd(sum),
      snr = abs(mean)/sd)
  write.csv(overall, "overall-summary.csv")
  
  for (stat in c("snr", "sd", "mean")) {  
    ggplot(overall, aes(x = nChains, y = get(stat), color = interaction(antithetics,essn))) +
      facet_grid(model + dim ~ objectiveType, scales = "free_y") +
      geom_line() + 
      theme_bw()
    ggsave(paste0("overall-", stat, ".pdf"), width = 15, height = 15)
  }  
  """
  
}