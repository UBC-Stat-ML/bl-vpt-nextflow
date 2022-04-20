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

nScans = 10000
nScans_ref = 20000
ks_threshold = 0.2

addModel('coll-rockets',   'p0',    ' --model ptbm.models.CollapsedHierarchicalRockets --model.data data/failure_counts.csv --engine.nChains 10 ')
addModel('toy-mix',        'x',     ' --model ptbm.models.ToyMix --engine.nChains 10 ')
addModel('mrna-no-transf', 'beta',  ' --model ptbm.models.MRNATransfectionNoTransform --model.data data/m_rna_transfection/processed.csv --engine.nChains 30 ')
addModel('sparse-car',     'alpha', ' --model ptbm.models.SparseCAR --model.data data/scotland_lip_cancer/data.csv --model.spatialData.adjacency data/scotland_lip_cancer/adj.csv --engine.nChains 20 ')
addModel('8-schools',      'tau',   ' --model ptbm.models.EightSchools --model.data data/eight-schools.csv --model.useTPrior true --engine.nChains 10 ')
addModel('titanic',        'sigma', ' --model ptbm.models.LogisticRegression --model.data data/titanic/titanic-covariates-original.csv --model.instances.name Name --model.instances.maxSize 20 --model.labels.dataSource data/titanic/titanic.csv --model.labels.name Survived --model.useTPrior false --engine.nChains 10 ')
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
  nScans = 250
  nScans_ref = 400
  seeds = (1..2)
  models = models.subList(0, 2)
  algos.retainAll{k, v -> ['V--T* F--T', reference].contains(k)}
}

process runMatching {

  input:
    each model from models
    each seed from seeds
    each algo from algos.entrySet()
    file code
    file data
    
  time '5h'
  errorStrategy 'ignore'
    
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
  
  cp results/latest/executionInfo/stdout.txt output
  cp results/latest/executionInfo/stdout.txt fixedRefOutput
  cp results/latest/executionInfo/stderr.txt output
  cp results/latest/executionInfo/stderr.txt fixedRefOutput
  mv results/latest/samples/${model.stat}.csv output/statistic.csv
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/fixedReferencePT/monitoring/*.csv fixedRefOutput
  mv results/latest/fixedReferencePT/samples/${model.stat}.csv fixedRefOutput/statistic.csv
  
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

results_all = results.concat(fixedRefResults)

process aggregate {
  time '1h'
  echo false
  scratch false
  input:
    file 'exec_*' from results_all.toList()
  output:
    file 'results/aggregated/' into aggregated
    file 'statistic.csv.gz' into statistic
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder lambdaInstantaneous.csv actualTemperedRestarts.csv globalLambda.csv logNormalizationConstantProgress.csv statistic.csv \
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

big_w = 17
big_h = 8
big_dims = "width = ${big_w}, height = ${big_h}"
colours = 'values = c(  "good" = "blue", "poor" = "red")'
custom_colours = "scale_fill_manual($colours) + scale_colour_manual($colours)"


process computeKS {
  time '1h'
  input:
    file statistic
  output:
    file 'ks_distances.csv' into ks_distances
    file 'statistic.pdf'
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
      facet_grid(algorithm ~ model, scales="free") +
      geom_density_ridges(panel_scaling = TRUE) +
      xlab("statistic") + 
      ylab("density") + 
      $custom_colours +
      theme_bw()
  ggsave("statistic.pdf", $big_dims)   
  """

}

process plot {
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
  
  ks_distances <- read.csv("${ks_distances}")
  ks_distances\$quality <- ifelse(ks_distances\$ks_stat > ${ks_threshold}, "poor", "good")
  
  restarts <- read.csv("${aggregated}/actualTemperedRestarts.csv.gz")
  restarts <- restarts %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))

  restarts %>%
    filter(algorithm != "Reference") %>%
    group_by(quality, model, algorithm) %>%
    summarize(total_count = sum(count) ) %>%
    ggplot(aes(x = algorithm, y = total_count, color = quality)) +
      facet_grid(model ~ ., scales="free_y") +
      geom_boxplot() +
      ylab("total tempered restarts") + 
      scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
      $custom_colours +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("actualTemperedRestarts-box.pdf", width = 5, height = ${big_h})
  
  global <- read.csv("${aggregated}/globalLambda.csv.gz")
  global <- global %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
  global %>%
    filter(fixedRefChain == "false") %>%
    ggplot(aes(x = round, y = value, color = quality, group = factor(seed))) +
      facet_grid(model ~ algorithm, scales="free_y") +
      $custom_colours +
      ylab("Lambda") + 
      geom_line() +
      theme_bw()
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


