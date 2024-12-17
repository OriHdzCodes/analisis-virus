#Evidencia 2
#Se uso la versión 4.4

install.packages("seqinr")
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("msa",force = TRUE)
library("msa")
library("Biostrings")
library(seqinr)
library(rentrez)
library(ape)
library(msa)
library(phylogram)
library(Biostrings)


#Función para obtener las secuencias complementarias
complementarias <- function(virus) {
  nucleotidos <- c("a", "c", "g", "t")
  complementarios <- c("t", "g", "c", "a")
  complemento <- chartr(paste(nucleotidos, collapse = ""), paste(complementarios, collapse = ""), virus)
  return(complemento)
}

secuencia90 <- read.fasta(file = "secuencia_90.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia90
secuencia89 <- read.fasta(file = "secuencia_89.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia89
secuencia88 <- read.fasta(file = "secuencia_88.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia88
secuencia87 <- read.fasta(file = "secuencia_87.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia87
secuencia86 <- read.fasta(file = "secuencia_86.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia86
secuencia85 <- read.fasta(file = "secuencia_85.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia85
secuencia84 <- read.fasta(file = "secuencia_84.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia84
secuencia83 <- read.fasta(file = "secuencia_83.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia83
secuencia82 <- read.fasta(file = "secuencia_82.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia82
secuencia81 <- read.fasta(file = "secuencia_81.fasta", as.string = FALSE, set.attributes = FALSE)
secuencia81

secuencia90_vector <- secuencia90$LQ980785.1
secuencia90_vector
secuencia89_vector <- secuencia89$LQ980784.1
secuencia89_vector
secuencia88_vector <- secuencia88$LQ980783.1
secuencia88_vector
secuencia87_vector <- secuencia87$LQ980782.1
secuencia87_vector
secuencia86_vector <- secuencia86$LQ980781.1
secuencia86_vector
secuencia85_vector <- secuencia85$LQ980780.1
secuencia85_vector
secuencia84_vector <- secuencia84$LQ980779.1
secuencia84_vector
secuencia83_vector <- secuencia83$LQ980778.1
secuencia83_vector
secuencia82_vector <- secuencia82$LQ980777.1
secuencia82_vector
secuencia81_vector <- secuencia81$LQ980776.1
secuencia81_vector

#Tamaño de cada secuencia 
longitud90 <- length(secuencia90_vector)
longitud90 
longitud89 <- length(secuencia89_vector)
longitud89 
longitud88 <- length(secuencia88_vector)
longitud88
longitud87 <- length(secuencia87_vector)
longitud87
longitud86 <- length(secuencia86_vector)
longitud86
longitud85 <- length(secuencia85_vector)
longitud85
longitud84 <- length(secuencia84_vector)
longitud84
longitud83 <- length(secuencia83_vector)
longitud83
longitud82 <- length(secuencia82_vector)
longitud82
longitud81 <- length(secuencia81_vector)
longitud81

secuencia90table <- seqinr::count(secuencia90_vector, 1)
secuencia89table <- seqinr::count(secuencia89_vector, 1)
secuencia88table <- seqinr::count(secuencia88_vector, 1)
secuencia87table <- seqinr::count(secuencia87_vector, 1)
secuencia86table <- seqinr::count(secuencia86_vector, 1)
secuencia85table <- seqinr::count(secuencia85_vector, 1)
secuencia84table <- seqinr::count(secuencia84_vector, 1)
secuencia83table <- seqinr::count(secuencia83_vector, 1)
secuencia82table <- seqinr::count(secuencia82_vector, 1)
secuencia81table <- seqinr::count(secuencia81_vector, 1)

#Crear grafica de frecuencia para comparar las secuencias de viruses
nucleotides <- matrix(0, nrow = 4, ncol = 10, dimnames = list(c("A", "C", "T", "G"), c("S90","S89","S88","S87","S86","S85","S84","S83","S82","S81")))

nucleotides[,"S90"] <- table(secuencia90_vector)
nucleotides[,"S89"] <- table(secuencia89_vector)
nucleotides[,"S88"] <- table(secuencia88_vector)
nucleotides[,"S87"] <- table(secuencia87_vector)
nucleotides[,"S86"] <- table(secuencia86_vector)
nucleotides[,"S85"] <- table(secuencia85_vector)
nucleotides[,"S84"] <- table(secuencia84_vector)
nucleotides[,"S83"] <- table(secuencia83_vector)
nucleotides[,"S82"] <- table(secuencia82_vector)
nucleotides[,"S81"] <- table(secuencia81_vector)

nucleotides_df <- as.data.frame(nucleotides)

nucleotides_df <- t(nucleotides_df)

barplot(nucleotides_df, beside = TRUE, col = c("lightblue", "red", "green", "orange","yellow","purple","pink","gold","blue","lightgreen"),
        main = "Tabla de Frecuencia",
        xlab = "Virus", ylab = "Frecuencia")
legend(x = "bottomleft", legend = c("S90", "S89","S88","S87","S86","S85","S84","S83","S82","S81"), fill = c("royalblue", "grey"), 
       title = "Loan")

#Sacar los alineamientos
variantes <- c(secuencia90_vector, secuencia89_vector, secuencia88_vector, secuencia87_vector, secuencia86_vector, 
               secuencia85_vector, secuencia84_vector, secuencia83_vector, secuencia82_vector, secuencia81_vector)
secuencias = DNAStringSet(variantes)
secuencias
AlineamientoSecuencias = msa(secuencias, method = "ClustalW")
AlineamientoSecuencias

distancias <- dist.hamming(as.character(consensusSequence(AlineamientoSecuencias)))

# Crear el árbol filogenético utilizando Neighbor-Joining
arbol <- nj(dist.dna(AlineamientoSecuencias))

# Graficar el árbol filogenético
plot(arbol, main = "Árbol Filogenético de las Secuencias")
