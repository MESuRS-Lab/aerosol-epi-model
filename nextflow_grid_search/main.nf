nextflow.enable.dsl=2

// parametres
params.simu_count = 100
params.folder_result = 'results'
params.folder_simu = 'simu'

// chemin entier vers le code R du modele
// sur le cluster
params.rlib= "cpp/dev-sensibility-analysis.cpp"
params.rlib2= "R/nodscov2/helper-functions-simulations.R"

// chemin vers le dossier ou se trouve les fichiers .RData input
params.rinput = "out/"


process generation_param {
  publishDir "${params.folder_result}", mode: 'copy' 

  output:
	path 'param_grid.txt' 

  script:
  """
  echo 'Generate parameter grid'
	Rscript ${workflow.projectDir}/bin/generation_param.R
  """
}

process simulation {
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
	Rscript ${workflow.projectDir}/bin/simu_nextflow.R $rlib $rlib2 $rinput ${simu_count} ${line}
  """
}


workflow {
  param_grid = generation_param() 
	lines = param_grid
		.splitText(keepHeader: true)

	rlib=file(params.rlib)
	rlib2=file(params.rlib2)
  rinput=file(params.rinput)
	simulation=simulation(rlib, rlib2, rinput, params.simu_count, lines)
	
  simulation.resu
		.collectFile(name: 'resu_simu_all.txt', newLine: true, keepHeader: true, skip: 1)
		.subscribe { file -> 
      // Read the file, filter out empty lines, and append to the existing file
      def cleanedFile = file.readLines().findAll { it.trim() }
        
      // Append the cleaned content to the output file
      def outputFile = new File("${params.folder_result}/resu_simu_all.txt")
        
      // Check if the file exists
      if (outputFile.exists()) {
        // Append the cleaned content to the file, adding a newline before appending
          outputFile.append("\n" + cleanedFile.join('\n'))
      } else {
        // If the file doesn't exist, create it and write the cleaned content
        outputFile.text = cleanedFile.join('\n')
      }
		}
}
