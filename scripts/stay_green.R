#
#	TRANSFORMS A DARTSEQ DATA TABLE TO HAPMAP FORMAT
#	TO PERFORM A GWAS ANALYSIS USING GAPIT SCRIPT.
#
#	AUTHOR: JAVIER MONTALVO-ARREDONDO
#	CONTACT: buitrejma@gmail.com
#	UNIVERSIDAD AUTONOMA AGRARIA ANTONIO NARRO
#	DEPARTAMENTO DE BIOTECNOLOGIA
#	CALZADA ANTONIO NARRO, 1923 BUENAVISTA SALTILLO, COAHUILA
#

# LOAD DARTSEQ DATA TABLE
#	HERE GOES THE FILE ----�
#		file = ?   # <------|
		data <- read.delim(file, header=T, sep=",")