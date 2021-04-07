#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')
data = file("data")

process buildCode {
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '52b3b10324bc5bfa83e49849d26477f07f6619f0'
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh'
}

seeds = (1..2)
Ns = (2..2).collect{Math.pow(2, it)}

normalModelArgs  = ' --engine.initialParameters 0.5'
normalModelArgs += ' --engine.pt.adaptFraction 0.0'
normalModelArgs += ' --engine.pt.nScans 1'
normalModelArgs += ' --engine.pt.initialization FORWARD'

process run {

  cpus 1
  time '1h'
  memory '2 GB'
  
  input:
    each seed from seeds
    each obj from 'Rejection' //, 'SKL', 'ApproxRejection'
    each antit from 'OFF' //, 'IS', 'MCMC'
    each essn from '0.5' //, '1.0'
    each model from 'ConjugateNormal' + normalModelArgs /*, 
                    'ToyNormal' + normalModelArgs, 
                    'CHRVariational' */
    each nChain from Ns
    file code
    file data
    
  output:
    file 'results/latest' into results
    
  """
  java -Xmx5g -cp code/lib/\\*  ptgrad.Variational \
    --experimentConfigs.resultsHTMLPage false \
    --model.interpolation $model \
    --engine.antithetics $antit \
    --engine ptgrad.VariationalPT \
    --engine.pt.nChains $nChain \
    --engine.pt.nThreads single \
    --engine.relativeESSNeighbourhoodThreshold $essn \
    --engine.pt.random $seed \
    --engine.objective $obj \
    --engine.nScansPerGradient 20 \
    --engine.optimizer.maxIters 0
  """
}

process analysisCode {
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'b53e3a302ceba4427c55b838afb3b8f3fdc23ec5'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}

process aggregate {
  echo true
  input:
    file analysisCode
    file 'exec_*' from results.toList()
  output:
    file 'results/latest/' into aggregated
  """
  code/bin/aggregate \
    --dataPathInEachExecFolder stochastic-gradient-evaluations.csv monitoring/roundTimings.csv \
    --keys \
      model.interpolation as model \
      engine.antithetics as antithetics \
      engine.objective as objectiveType \
      engine.pt.nChains as nChains \
      engine.relativeESSNeighbourhoodThreshold as essn \
      engine.pt.random as random \
           from arguments.tsv
  """
}

process plot {
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
    ggplot(overall, aes(x = nChains, y = get(stat), linetype = objectiveType, color = interaction(antithetics,essn))) +
      facet_grid(model + dim ~ ., scales = "free_y") +
      geom_line() + 
      theme_bw()
    ggsave(paste0("overall-", stat, ".pdf"), width = 5, height = 15)
  }  
  """
  
}