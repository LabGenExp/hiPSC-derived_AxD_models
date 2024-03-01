setwd("d:/20211221_Lund_AxD/")

#################CREATION WORKING FOLDER###############################################
cat(paste0("mkdir qa;\n",
             "mkdir qa/raw;\n",
             "mkdir qa/fastq_screen;\n",
             "mkdir log;\n",
             "mkdir STAR;\n"))

samples<-read.table("samples.csv",header=T,sep=",")
##################QUALITY CONTROL BEFORE ANALYSIS######################################

for(i in 1:nrow(samples)) {
  cat(paste0("zcat ",
              samples$LibraryName[i],"_","L001_R1_001.fastq.gz ",
              samples$LibraryName[i],"_","L002_R1_001.fastq.gz > ",
              samples$LibraryName[i],"_1.fastq;\n",
              "pigz -p 12 ",samples$LibraryName[i],"_1.fastq;\n"))
  
  cat(paste0("zcat ",
              samples$LibraryName[i],"_","L001_R2_001.fastq.gz ",
              samples$LibraryName[i],"_","L002_R2_001.fastq.gz > ",
              samples$LibraryName[i],"_2.fastq;\n",
              "pigz -p 12 ",samples$LibraryName[i],"_2.fastq;\n"))
}



cat(paste0("fastqc -t 10 -o qa/raw raw/*;\n"))

#################scRNA-Seq#################################################
################DEFINITION OF VARIABLES################################################

#UMI and cell barcode length
CBlength = c(16)
UMIlength = c(12)
expectedCells = c(10000)
#10x barocde list
BarcodeList_sc = c("/mnt/d/10x_barcodes/3M-february-2018.txt")
#GenomeDir
genomeDir1=c("/mnt/f/Genomes/Homo_sapiens_GRCh38_single_cell")
#NumberOfThreads
threads=c(10)

#############ANALYSIS SAMPLE AFTER SAMPLE##############################################
for(i in 1:nrow(samples)) {
  #control of contaminations
  cat(paste0("~/fastq_screen_v0.11.1/fastq_screen --aligner bwa --threads ",threads,
              " --conf /mnt/f/Genomes/fastq_screen/fastq_screen.conf",
              " raw/",samples$LibraryName[i],"_2.fastq.gz;\n"))
  
  cat(paste0("mv ",samples$LibraryName[i],"_2_screen* qa/fastq_screen;\n"))
  
  cat(paste0("unpigz -p ",threads," raw/",samples$LibraryName[i],"_2.fastq.gz;\n"))
  cat(paste0("unpigz -p ",threads," raw/",samples$LibraryName[i],"_1.fastq.gz;\n"))
  
  
  #mapping
  cat(paste0("mkdir STAR/",samples$SampleName[i],";\n"))
  
  cat(paste0("~/STAR-2.7.9a/bin/Linux_x86_64/STAR",
              " --genomeDir ",genomeDir1,"/",
              " --readFilesIn",
              " raw/",samples$LibraryName[i],"_2.fastq",
              " raw/",samples$LibraryName[i],"_1.fastq",
              " --soloType CB_UMI_Simple",
              " --runThreadN ",threads,
              " --soloCBwhitelist ",BarcodeList_sc,
              " --soloCBstart 1",
              " --soloCBlen ",CBlength,
              " --soloUMIstart ",CBlength + 1,
              " --soloUMIlen ",UMIlength,
              " --soloBarcodeReadLength 0",
              " --soloFeatures Gene GeneFull SJ Velocyto",
              " --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts",
              " --soloUMIdedup 1MM_Directional_UMItools",
              " --soloUMIfiltering MultiGeneUMI",
              " --soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx",
              " --soloCellFilter EmptyDrops_CR ",expectedCells," 0.99 10 45000 90000 100 0.01 20000 0.001 10000",
              " --outFileNamePrefix STAR/",samples$SampleName[i],"/",samples$SampleName[i],
              " --outFilterMultimapNmax 1",
              " --outSAMtype BAM SortedByCoordinate",
              " --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM",
              " --soloMultiMappers Uniform PropUnique EM Rescue",
              " --clipAdapterType CellRanger4;\n"))
  
  cat(paste0("pigz -p ",threads," raw/",samples$LibraryName[i],"_2.fastq;\n"))
  cat(paste0("pigz -p ",threads," raw/",samples$LibraryName[i],"_1.fastq;\n"))
  cat(paste0("cp STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],"Log.final.out log/",as.vector(samples$SampleName)[i],"_STAR.log;\n"))
}

