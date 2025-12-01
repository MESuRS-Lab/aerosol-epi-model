nextflow.enable.dsl=2

// parametres
params.simu_count = 500
params.folder_result = 'results'
params.folder_simu = 'simu'

// Paths to R codes
params.rlib= "cpp/dev-sensibility-analysis.cpp"
params.rlib2= "R/nodscov2/helper-functions-simulations.R"

// Path to input data
params.rinput = "out/"



process generation_param {
  publishDir "${params.folder_result}", mode: 'copy'

  output:
	path 'param_grid_interventions.txt' 

  script:
  """
  echo 'Generate parameter grid'
	Rscript ${workflow.projectDir}/bin/generation_param_interventions.R 
  """
}

process intervention_simu {
  scratch true
	
	publishDir "${params.folder_result}/${params.folder_simu}", mode: 'copy'
  
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  
	input:
	path rlib
	path rlib2
	path rinput
  val simu_count
  val line
	
	output:
	path 'summary_stat*.csv', emit: resu

  script:
  """
	Rscript ${workflow.projectDir}/bin/interventions.R $rlib $rlib2 $rinput ${line}
  """
}


workflow {
  param_grid = generation_param() 
	lines = param_grid
		.splitText(keepHeader: true)

  rlib=file(params.rlib)
  rlib2=file(params.rlib2)
  rinput=file(params.rinput)
  intervention_simu=intervention_simu(rlib, rlib2, rinput, params.simu_count, lines)

  intervention_simu.resu
    .collectFile(name: 'resu_interventions_all.txt', newLine: true, keepHeader: true, skip: 1)
    .subscribe { file -> 
      def cleanedFile = file.readLines().findAll { it.trim() }
        
      def outputFile = new File("${params.folder_result}/resu_interventions_all.txt")
        
      if (outputFile.exists()) {
          outputFile.append("\n" + cleanedFile.join('\n'))
      } else {
        outputFile.text = cleanedFile.join('\n')
      }
    }
}

