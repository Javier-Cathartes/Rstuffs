#
#	CALCULATES THE FREQUENCY OF EACH CODON
#	BY TOTAL CODONS FOR EVERY AMINO ACID
#	
#	AUTHOR: JAVIER MONTALVO
#	CONTACT: buitrejma@gmail.com
#	UNIVERSIDAD AUTONOMA AGRARIA ANTONIO NARRO
#	CALZADA ANTONIO NARRO 1923 BUENAVISTA SALTILLO COAHUILA
#	DEPARTAMENTO DE BIOTECNOLOGIA
#

# WRITE-FQTABLE-TO-FILE METHOD								#################
# NOTE THAT, IT DEPENDS ON CODON_TABLE & GCLIDX OBJECTS,	#   HARDCODED   #
# DO NOT ERASE THEM OR DO NOT CHANGE THE NAME				#################
		write.FQTable <- function(fqt, outfile) {
			gene_id <- fqt@ID
			indices <- 1:length(gene_id)
			for (idx in indices) {
				write(paste(">", gene_id[idx], sep=""), file=outfile, append=T)
				table <- cbind(
								attr(fqtable@counts, "dimnames")[[2]],
								codon_table@counts[
													idx,
													as.vector(unlist(GCLidx))
													],
								fqtable@counts[idx,])
				table <- cbind(rownames(table), table)
				colnames(table) <- c(
										"Amino acid",
										"codon",
										"counts",
										"frequency")
				
				write.table(
								x = table,
								file = outfile,
								append=T,
								sep="\t",
								col.names=T,
								row.names=F,
								quote=F
							 	)
				}
		}

# VALIDITY DATA FUNCTIONS
		check_data <- function(object) {
			if (length(object@counts) == 0) {
				stop("FQTable matrix is emtpy!")
			}
		}

# DECLARE A FQTABLE OBJECT
		create_a_FQTable <- setClass("FQTable", representation(
															ID="character",
															counts="matrix",
															len="numeric",
															#ImTheCat
															KO="character",
															COG="character"),
										validity = check_data)

# LOAD SEQUENCES
# CALCULATE CODON USAGE TABLE
		library(coRdon)
#		Here goes the INPUT filename! UNMASK IT  #---¬
#										                V
#		dna <- readSet(file="https://raw.githubusercontent.com/BioinfoHR/coRdon-examples/master/LD94.fasta")
		if (!"dna" %in% ls()){
			stop("DNA sequences file unavailable!")
		}
		codon_table <- codonTable(dna)

# LOAD GENETIC CODE
		library(Biostrings)
		GC <- getGeneticCode()
		GC <- sort(GC)

# PREPARE GENETIC CODE HASH TABLE
		GCL <- sapply(names(table(GC)), function (x) {
											names(GC[GC == x])
											})

# GET INDICES OF COUNTS MATRIX WHERE COLNAMES
# MATCH WITH CODONS BY AMINO ACID
		GCLidx <- sapply(names(GCL), function(x){
			which(attr(codon_table@counts, "dimnames")[[2]] %in% GCL[[x]])
		})

# CALCULATE THE FRECUENCY OF CODONS BY AMINO ACIDS
# Fac = (Oac*100) / SUM(Oac)
		m <- c()
		for (name in names(GCLidx)) {
			tmp_m <- codon_table@counts[,GCLidx[[name]]]
			if (is.matrix(tmp_m)){
				rs <- rowSums(tmp_m)
			}
			else {
				rs <- tmp_m
			}
			tmp_m <- ((tmp_m * 100) / rs)
			m <- cbind(m, round(tmp_m, 2)) #ImTheCat
		}
		colnames(m) <- unlist(GCL)
		fqtable <- create_a_FQTable(
										ID=codon_table@ID,
										counts=m,
										len=codon_table@len,
										#ImTheCat
										KO=codon_table@KO,
										COG=codon_table@COG
										)
# FORMAT RESULTS TO OUTPUT TABLE
		attr(fqtable, "class") <- "codonTable"

# WRITE TO FILE THE FQTABLE
#    		HERE GOES THE OUTPUT FILENAME	-------¬
#                                               V
		write.FQTable(fqt=fqtable, outfile="out.fqt")

# THAT'S ALL FOLKS

#//~CAT
