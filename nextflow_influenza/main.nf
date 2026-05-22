// Permet de lancer des simulations en parallele a partir d une grille de parametre qui va etre genere
// TUTO
// pour faire tourner la pipeline nextflow, il faut se placer dans un dossier et avoir dedans ce script, le script nextflow.config, et un dossier bin avec les codes R
// il faut également le fichier code_model.R et renseigner son path dans ce fichier
// on lance la ligne de commande suivante : nextflow run main.nf --simu_count 2
// parametre simu_count : nombre de simulations effectuees avec un jeu de parametre

// pour modifier les valeurs des parametres de la grille, il faut modifier le fichier generation_param.R
// pour modifier les parametres à faire varier, il faut modifier les fichiers generation_param.R, simu_nextflow.R et stats_simu.R
// pour modifier les indicateurs des simulations, il faut modifier le fichier simu_nextflow.R et stats_simu.R

// pour utiliser sur le cluster de Pasteur, il faut d abord charger les modules avec : module load apptainer, module load graalvm et module load nextflow


nextflow.enable.dsl=2

// parametres
params.simu_count = 500
params.folder_result = 'results'
params.folder_simu = 'simu'

// chemin entier vers le code R du modele
params.rlib= "/pasteur/helix/users/maylayan/Pacri/Hospitals_modes_transmission/github_codes/cpp/dev-sensibility-analysis.cpp"
params.rlib2= "/pasteur/helix/users/maylayan/Pacri/Hospitals_modes_transmission/github_codes/R/nodscov2/helper-functions-simulations.R"

// chemin vers le dossier ou se trouve les fichiers .RData input
params.rinput = "/pasteur/helix/users/maylayan/Pacri/Hospitals_modes_transmission/github_codes/out/"



process generation_param {
  // permet de generer la grille de parametre a partir du nombre de simulation que l on veut effectuer pour chaque jeu de parametre
  publishDir "${params.folder_result}", mode: 'copy' 

  output:
	path 'param_grid_influenza.txt' 

  script:
  """
  echo 'Generate parameter grid'
	Rscript ${workflow.projectDir}/bin/generation_param_influenza.R 
  """
}

process sensitivity_simu {
  scratch true
	// realise une simulation a partir d un jeu de parametre et calcul des indicateurs interessants
	publishDir "${params.folder_result}/${params.folder_simu}", mode: 'copy'
  
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  
	input:
	path rlib
	path rlib2
	path rinput
  val simu_count
  val line 
	
	output:
	path 'summary_stat*.csv', emit: resu //fichier .txt avec une seule ligne avec tous les indicateurs de la simulation

  script:
  """
	Rscript ${workflow.projectDir}/bin/influenza.R $rlib $rlib2 $rinput ${line}
  """
}


workflow {
	// genere la grille de set de parametre 
  param_grid = generation_param() 
	
	// permet de split chaque ligne de param_grid 
	lines = param_grid
		.splitText(keepHeader: true)

  // fait tourner la simulation sur chaque set de parametre
  rlib=file(params.rlib)
  rlib2=file(params.rlib2)
  rinput=file(params.rinput)
  sensitivity_simu=sensitivity_simu(rlib, rlib2, rinput, params.simu_count, lines)

  // collecter tous les indicateurs de chaque simu dans un seul fichier que l on enregistre dans le dossier de resultat
  sensitivity_simu.resu
    .collectFile(name: 'resu_influenza_all.txt', newLine: true, keepHeader: true, skip: 1)
    .subscribe { file -> 
      // Read the file, filter out empty lines, and append to the existing file
      def cleanedFile = file.readLines().findAll { it.trim() }
        
      // Append the cleaned content to the output file
      def outputFile = new File("${params.folder_result}/resu_influenza_all.txt")
        
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

