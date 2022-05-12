#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '9bbcb268dde6d746ea0ae1b044e0cca34aa588c0' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}


nScans = 1000
nScans_ref = 5000


algos = [:]  
reference = 'Reference'
algos[reference]  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational INF '
algos['V--T*--F']    = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational true '

postprocessor = ' --postProcessor ptgrad.VariationalPostprocessor '

params.dryRun = false

if (params.dryRun) {
  nScans = 500
  nScans_ref = 500
}

process runMatching {

  input:
    each algo from algos.entrySet()
    file code
    file data
        
  output:
    file 'output' into results
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    $postprocessor \
    --engine ptbm.OptPT \
    --engine.storeSamplesForAllChains true \
    --engine.random 1 \
    --engine.nScans ${if (algo.key == reference) nScans_ref else nScans} \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    --engine.useFixedRefPT true \
    ${algo.value} \
    --model ptbm.models.MRNATransfectionNoTransformForPlotting  \
    --model.data data/m_rna_transfection/processed.csv \
    --engine.nChains 30 \
    --engine.nThreads single
  mkdir output
  
  cp results/latest/executionInfo/stdout.txt output
  cp results/latest/executionInfo/stderr.txt output
  mv results/latest/samplesForAllChains/beta.csv output/statistic.csv
  
  echo "\nalgorithm\t${algo.key}" >> results/latest/arguments.tsv
  echo "path\t\$(pwd)" >> results/latest/arguments.tsv
  cp results/latest/*.tsv output
  
  """ 

}


process aggregate {
  echo false
  scratch false
  input:
    file 'exec_*' from results.toList()
  output:
    file 'statistic.csv.gz' into statistic
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder \
        statistic.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      algorithm \
           from arguments.tsv
  mv results/latest results/aggregated
  mv results/aggregated/statistic.csv.gz .
  """
}

process plots {
  scratch false
  input:
    file statistic
  output:
    file '*.*'
  afterScript 'rm Rplots.pdf; cp .command.sh rerun.sh'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  
  betas <- read.csv("${statistic}")
  
  betas %>% 
    filter(algorithm == 'V--T*--F') %>%
    filter(sample > ${nScans/2}) %>% 
    ggplot(aes(x = value, colour = chain, group = chain)) + 
      geom_density() +
      theme_minimal() + 
      xlab("parameter") + 
      ylab("density") + 
      theme(legend.position="none") +
      xlim(-0.4, 1.2) + 
      ylim(0, 7) +
      scale_colour_gradient(low = "grey30", high = "skyblue1")   
  ggsave("variational-paths.pdf", width = 2, height = 2)
  
  betas %>% 
    filter(algorithm == 'Reference') %>%
    filter(sample > ${nScans_ref/2}) %>% 
    ggplot(aes(x = value, colour = chain, group = chain)) + 
      geom_density() +
      theme_minimal() +
      xlab("parameter") + 
      ylab("density") + 
      theme(legend.position="none") +
      xlim(-0.4, 1.2) +
      ylim(0, 7) +  
      scale_colour_gradient(low = "grey30", high = "orange") 
  ggsave("fixed-ref-paths.pdf", width = 2, height = 2)
    
  """  
}



