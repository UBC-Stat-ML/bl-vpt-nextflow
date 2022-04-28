#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')

process buildCode {
  executor 'local'
  cache true 
  input:
    val gitRepoName from 'ptanalysis'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'ed74dca2bf28d35b9e514ab2d02012c6bec7750c' 
    val snapshotPath from "${System.getProperty('user.home')}/w/ptanalysis"
  output:
    file 'code' into code
    file 'ptanalysis/data' into data
  script:
    template 'buildRepo.sh'
}

seeds = (1..10)
sizes = (2..13).collect{Math.pow(2, it)}
nScans = 2000

models = []

class Model {
  String name
  String sizeArg
  String args
  def int hashCode() { return [name, sizeArg, args].hashCode() }
}

def addModel(String n, String s, String a) {
  m = new Model(name: n, sizeArg: s, args: a)
  models.add(m)
}


addModel('coll-rockets',    ' --model.rocketTypes.maxSize ', ' --model ptbm.models.CollapsedHierarchicalRockets --model.data data/failure_counts.csv --engine.nChains 10 ')
addModel('unidentiable',    ' --model.nTrials ',             ' --model ptbm.models.UnidentifiableProduct ')
addModel('Cauchy-Cauchy',   ' --model.obs.maxSize ',         ' --model ptbm.models.CauchyCauchy --model.data data/cc-100k/ys.csv ')

postprocessor = ' --postProcessor ptgrad.VariationalPostprocessor '

params.dryRun = false

if (params.dryRun) {
  seeds = seeds.subList(0, 2)
  sizes = sizes.subList(0, 2)
  models = models.subList(0, 1)
  nScans = 300
}

process runMatching {

  input:
    each model from models
    each seed from seeds
    each size from sizes
    each isVariational from true, false
    file code
    file data
    
  time '1h'
  errorStrategy 'ignore'
    
  output:
    file 'output' into results
  """
  java -Xmx5g -cp code/lib/\\*  blang.runtime.Runner \
    --experimentConfigs.resultsHTMLPage false \
    $postprocessor \
    --engine ptbm.OptPT \
    --engine.random $seed \
    ${model.sizeArg} $size \
    --engine.nScans $nScans \
    --engine.scmInit.nParticles 10 \
    --engine.scmInit.temperatureSchedule.threshold 0.9 \
    --engine.nPassesPerScan 1 \
    --engine.useFixedRefPT true \
    --engine.minSamplesForVariational ${if (isVariational) "100" else "INF"} \
    ${model.args} \
    --engine.nThreads single
    
  mkdir output
  
  cp results/latest/executionInfo/stdout.txt output
  cp results/latest/executionInfo/stderr.txt output
  mv results/latest/*.csv output
  mv results/latest/monitoring/*.csv output
  mv results/latest/ess/allEss.csv output
  
  echo "\nmodelDescription\t${model.name}" >> results/latest/arguments.tsv
  echo "isVariational\t${isVariational}" >> results/latest/arguments.tsv
  echo "size\t${size}" >> results/latest/arguments.tsv
  echo "path\t\$(pwd)" >> results/latest/arguments.tsv
  cp results/latest/*.tsv output  
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
  """
  aggregate \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder \
        globalLambda.csv \
        roundTimings.csv \
        allEss.csv \
    --experimentConfigs.tabularWriter.compressed true \
    --keys \
      modelDescription as model \
      engine.random as seed \
      isVariational \
      size \
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
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("dplyr")
  
  timing <- read.csv("${aggregated}/roundTimings.csv.gz") %>% rename(time = value)
  
  
  global <- read.csv("${aggregated}/globalLambda.csv.gz")
  global <- global %>% inner_join(timing, by = c("round", "model", "isVariational", "seed", "size"))
  global %>% 
    filter(isAdapt == "false") %>%
    filter(model != "coll-rockets" | size < 369) %>%
    group_by(size, isVariational, model, round) %>%
    summarize(
      mean_gcb = mean(value),
      se_gcb = sd(value)/sqrt(n())) %>%
    ggplot(aes(x = size, y = mean_gcb, colour = factor(isVariational))) +
      geom_line() + 
      geom_errorbar(aes(ymin=mean_gcb-se_gcb, ymax=mean_gcb+se_gcb), width=.1) +
      scale_x_log10() +
      scale_y_log10() + 
      facet_grid(. ~ model) +
      theme_bw()
      
  ggsave("scaling.pdf")
  """  
}


