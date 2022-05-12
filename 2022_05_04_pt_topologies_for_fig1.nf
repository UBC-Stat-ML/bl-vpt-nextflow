#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '51907cf40c409ffa352f404bf538ee0a6b3da4ac' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)

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

nScans = 40000
nScans_ref = 80000
ks_threshold = 0.1

addModel('mrna-no-transf', 'beta',  ' --model ptbm.models.MRNATransfectionNoTransform --model.data data/m_rna_transfection/processed.csv --engine.nChains 30 ')


algos = [:]  
reference = 'Reference'
algos['F--T  F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational INF '
algos[reference]     = algos['F--T  F--T'] 
algos['V--T* F--T']  = ' --engine.fullyIndepFixedRef true  --engine.minSamplesForVariational 100 '
algos['V--T*--F']    = ' --engine.fullyIndepFixedRef false --engine.minSamplesForVariational 100 --engine.doSwapFixedRefAndVariational true '

postprocessor = ' --postProcessor ptgrad.VariationalPostprocessor '

params.dryRun = false

if (params.dryRun) {
  nScans = 500
  nScans_ref = 500
  seeds = (1..1)
  models = models.subList(0, 1)
  //algos.retainAll{k, v -> ['V--T* F--T', reference].contains(k)}
}

process runMatching {

  input:
    each model from models
    each seed from seeds
    each algo from algos.entrySet()
    file code
    file data
    
  time '8h'
  errorStrategy 'ignore'
    
  output:
    file 'output' into results
    file 'fixedRefOutput' into fixedRefResults
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    $postprocessor \
    --engine ptbm.OptPT \
    --engine.storeSamplesForAllChains true \
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
  
  cp results/latest/executionInfo/stdout.txt output
  cp results/latest/executionInfo/stdout.txt fixedRefOutput
  cp results/latest/executionInfo/stderr.txt output
  cp results/latest/executionInfo/stderr.txt fixedRefOutput
  mv results/latest/samples/${model.stat}.csv output/statistic.csv
  mv results/latest/*.csv output
  mv results/latest/samplesForAllChains/*.csv output
  mv results/latest/fixedReferencePT/samplesForAllChains/*.csv fixedRefOutput
  mv results/latest/monitoring/*.csv output
  mv results/latest/fixedReferencePT/monitoring/*.csv fixedRefOutput
  mv results/latest/fixedReferencePT/samples/${model.stat}.csv fixedRefOutput/statistic.csv
  mv results/latest/ess/allEss.csv output
  mv results/latest/fixedReferencePT/ess/allEss.csv fixedRefOutput
  
  echo "\nmodelDescription\t${model.name}" >> results/latest/arguments.tsv
  echo "algorithm\t${algo.key}" >> results/latest/arguments.tsv
  echo "path\t\$(pwd)" >> results/latest/arguments.tsv
  echo "statisticName\t${model.stat}" >> results/latest/arguments.tsv
  cp results/latest/*.tsv output
  cp results/latest/*.tsv fixedRefOutput
  echo "fixedRefOutput\tfalse" >> output/arguments.tsv
  echo "fixedRefOutput\ttrue"  >> fixedRefOutput/arguments.tsv
  
  """ 

}


process aggregate {
  time '1h'
  echo false
  scratch false
  input:
    file 'exec_*' from results.toList()
  output:
    file 'results/aggregated/' into aggregated
    file 'statistic.csv.gz' into statistic
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder \
        lambdaInstantaneous.csv \
        actualTemperedRestarts.csv \
        beta.csv \
        globalLambda.csv \
        logNormalizationConstantProgress.csv \
        statistic.csv \
        allEss.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      modelDescription as model \
      engine.random as seed \
      algorithm \
      fixedRefOutput as fixedRefChain \
           from arguments.tsv
  mv results/latest results/aggregated
  mv results/aggregated/statistic.csv.gz .
  """
}

big_w = 3
big_h = 3
big_dims = "width = ${big_w}, height = ${big_h}"
colours = 'values = c(  "good" = "blue", "poor" = "red")'
custom_colours = "scale_fill_manual($colours) + scale_colour_manual($colours)"


process computeKS {
  time '5h'
  input:
    file statistic
  output:
    file 'ks_distances.csv' into ks_distances
    file 'statistic.pdf'
    file 'trace-plots.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  require("ggridges")
  
  stat <- read.csv("${statistic}")
  
  ks_distances <- data.frame(
    algorithm = factor(),
    model = factor(),
    seed = integer(),
    ks_p = double(),
    ks_stat = double()
  )
  
  for (a in unique(stat\$algorithm)) {
    for (m in unique(stat\$model)) {
      for (s in unique(stat\$seed)) {
        ref <- stat %>%
          filter(sample > ${nScans_ref/2}) %>%
          filter(algorithm == "Reference") %>%
          filter(fixedRefChain == "false") %>%
          filter(model == m) %>%
          .\$value
        cur <- stat %>%
          filter(sample > ${nScans/2}) %>%
          filter(algorithm == a) %>%
          filter(fixedRefChain == "false") %>%
          filter(model == m) %>%
          filter(seed == s) %>%
          .\$value
        if (length(ref) < 10 | length(cur) < 10) {
          ks_distances <- ks_distances %>% add_row(algorithm = a, model = m, seed = s, ks_p = 0, ks_stat = 1)
        } else {
          test <- ks.test(ref, cur)
          ks_distances <- ks_distances %>% add_row(algorithm = a, model = m, seed = s, ks_p = test\$p.value, ks_stat = as.numeric(test\$statistic))
        }
      }
    }
  }
  write.csv(ks_distances, 'ks_distances.csv')
  
  ks_distances\$quality <- ifelse(ks_distances\$ks_stat > ${ks_threshold}, "poor", "good")
  stat <- stat %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  
  remove_outliers <- function(x) {
    qnt <- quantile(x, probs=c(.01, .99))
    H <- 1.5 * IQR(x)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  stat %>%
    filter(sample > ${nScans/2}) %>%
    filter(fixedRefChain == "false") %>%
    group_by(seed, algorithm, model) %>%
    mutate(value = remove_outliers(value)) %>%
    ggplot(aes(x = value, y = factor(seed), fill = quality)) +
      facet_grid(algorithm ~ ., scales="free") +
      geom_density_ridges(panel_scaling = TRUE) +
      xlab("statistic") + 
      ylab("density") + 
      coord_flip() +
      $custom_colours +
      theme_classic()
  ggsave("statistic.pdf", $big_dims) 
  
  stat %>%
    filter(algorithm != 'Reference') %>%
    filter(seed == 1) %>%
    filter(sample < 1000) %>%
    ggplot(aes(x = sample, y = value, colour = quality)) +
      geom_point(size = 0.1) + geom_line(alpha = 0.5) + 
      theme_classic() + 
      ylim(-0.4, 1.2) + 
      facet_grid(algorithm ~ .) +
      $custom_colours +
      xlab("MCMC iteration") +
      ylab("sample") +
      theme(legend.position="none") +
  ggsave("trace-plots.pdf", width = 1.7, height = 2)  
  """

}


process subsampleSeeds {
  input:
    file aggregated
  output:
    file 'subsampled.csv.gz' into subsampled
  """
  #!/usr/bin/env Rscript
  require("dplyr")
  
  data <- read.csv("${aggregated}/beta.csv.gz")
  
  subsampled <- data %>% 
    filter(seed == 1)
    
  write.csv(subsampled, "subsampled.csv.gz")
  """
}

process costlyPlots {
  scratch false
  input:
    file aggregated
    file ks_distances
    file subsampled
  output:
    file '*.*'
    file 'aggregated'
  afterScript 'rm Rplots.pdf; cp .command.sh rerun.sh'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  require("ggridges")
  
 
  betas <- read.csv("${subsampled}")
  
  betas %>% 
    filter(algorithm == 'V--T*--F') %>%
    filter(sample > ${nScans/2}) %>% 
    ggplot(aes(x = value, colour = chain, group = chain)) + 
      geom_density() +
      theme_classic() + 
      xlab("parameter") + 
      ylab("density") + 
      theme(legend.position="none") +
      xlim(-0.4, 1.2) + 
      scale_colour_gradient(low = "grey30", high = "skyblue1")   
  ggsave("variational-paths.pdf", width = 2, height = 2)
  
  betas %>% 
    filter(algorithm == 'Reference') %>%
    filter(sample > ${nScans_ref/2}) %>% 
    ggplot(aes(x = value, colour = chain, group = chain)) + 
      geom_density() +
      theme_classic() +
      xlab("parameter") + 
      ylab("density") + 
      theme(legend.position="none") +
      xlim(-0.4, 1.2) +  
      scale_colour_gradient(low = "grey30", high = "orange") 
  ggsave("fixed-ref-paths.pdf", width = 2, height = 13, limitsize = FALSE)
    
  """  
}

process cheapPlots {
  scratch false
  input:
    file aggregated
    file ks_distances
  output:
    file '*.*'
    file 'aggregated'
  afterScript 'rm Rplots.pdf; cp .command.sh rerun.sh'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  require("ggridges")
  require("forcats")
  
  ks_distances <- read.csv("${ks_distances}")
  ks_distances\$quality <- ifelse(ks_distances\$ks_stat > ${ks_threshold}, "poor", "good")
  
  
  restarts <- read.csv("${aggregated}/actualTemperedRestarts.csv.gz")
  restarts <- restarts %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  
  restarts\$order <- 
    ifelse(
      restarts\$algorithm == "F--T F--T", 
      3,
      ifelse(
        restarts\$algorithm == "V--T* F--T",
        2,
        1
      )
    )
      
  restarts %>%
    filter(algorithm != "Reference") %>%
    group_by(quality, model, algorithm, seed, order) %>%
    summarize(total_count = sum(count) ) %>%
    ggplot(aes(x = reorder(algorithm,order), y = total_count, color = quality)) +
      facet_grid(model ~ ., scales="free_y") +
      geom_boxplot(outlier.alpha = 0.0) +
      ylab("total tempered restarts") + 
      scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
      $custom_colours +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("actualTemperedRestarts-box.pdf", width = 4, height = 3)
    
  ess <- read.csv("${aggregated}/allEss.csv.gz")
  ess <- ess %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  
  ess\$order <- 
    ifelse(
      ess\$algorithm == "F--T F--T", 
      3,
      ifelse(
        ess\$algorithm == "V--T* F--T",
        2,
        1
      )
    )
  
  ess %>%
    filter(algorithm != "Reference") %>%
    filter(variable %in% c(${models.stream().map{m -> "\"" + m.stat + "\""}.toList().join(",")})) %>%
    group_by(quality, model, algorithm, variable, order) %>%
    summarize(total_ess = sum(value) ) %>%
    ggplot(aes(x = reorder(algorithm,order), y = total_ess, color = quality)) +
      facet_grid(model ~ ., scales="free_y") +
      geom_boxplot() +
      ylab("ESS") + 
      scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
      $custom_colours +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("ess-box.pdf", width = 4, height = 6)
  
  global <- read.csv("${aggregated}/globalLambda.csv.gz")
  global <- global %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  global %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = round, y = value, color = quality, group = factor(seed))) +
      facet_grid(model ~ algorithm, scales="free_y") +
      $custom_colours +
      ylab("Lambda") + 
      geom_line() +
      theme_classic()
  ggsave("globalLambda.pdf", $big_dims)  
  
  lambda <- read.csv("${aggregated}/lambdaInstantaneous.csv.gz")
  lambda <- lambda %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  lambda %>% 
    filter(isAdapt == "false") %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = beta, y = value, color = quality, group = factor(seed))) +
      facet_grid(model ~ algorithm, scales="free_y") +
      $custom_colours +
      geom_line(alpha = 0.5) +
      xlab("beta") + 
      ylab("intensity") + 
      theme_bw()
  ggsave("lambda.pdf", $big_dims)
  """  
}


