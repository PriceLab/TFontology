####################################################################
########### Transcription Factor Ontologies Reference ##############
########### Alison Paquette 1/14/2018                 ##############
####################################################################

####Load Packages
library(MotifDb)
library(GO.db)
library(org.Hs.eg.db)

####################################################
####### Part 1:Motif DB Gene Motifs ################
####################################################

Motifs<-query(MotifDb,TFs_GO)

mdb.human=query(MotifDb,"hsapiens",c("HOCOMOCOv10","jaspar2018","SwissRegulon"))

#Note: All otpions ( cisbp_1.02:  313
      # HOCOMOCOv10:  640
            #hPDI:  437
# JASPAR_2014:  117
 # JASPAR_CORE:   66
  # jaspar2016:  442
  #jaspar2018:  537
   #jolma2013:  710
     # stamlab:  683
   #SwissRegulon:  684
     #  UniPROBE:2
genesymbols=unique(mcols(mdb.human)$geneSymbol)
Motif_Table<-as.data.frame(mcols(mdb.human))
save(Motif_Table,file="/Users/alisonpaquette/TrenaProjectPlacenta/Motif_Table.Rdata")

#Make Gene Specific Data Frame
Motif_GeneDF<-as.data.frame(c(genesymbols))
colnames(Motif_GeneDF)<-"GeneSymbol"
rownames(Motif_GeneDF)<-genesymbols
Motif_GeneDF$HOCOMOCOv10<-rep(NA,length(genesymbols))
Motif_GeneDF$SwissRegulon<-rep(NA,length(genesymbols))
Motif_GeneDF$jaspar2018<-rep(NA,length(genesymbols))

#Test to see if Motifs from individual Databases are in Gene DF
for(i in 1:length(Motif_GeneDF$GeneSymbol)){
  temp<-subset(Motif_Table,geneSymbol==as.character(Motif_GeneDF$GeneSymbol[i]))
    if(length(grep("HOCOMOCOv10",temp$dataSource)>0)){Motif_GeneDF$HOCOMOCOv10[i]=T} else {Motif_GeneDF$HOCOMOCOv10[i]=F}
    if(length(grep("SwissRegulon",temp$dataSource))>0){Motif_GeneDF$SwissRegulon[i]=T} else {Motif_GeneDF$SwissRegulon[i]=F}
    if(length(grep("jaspar2018",temp$dataSource))>0){Motif_GeneDF$jaspar2018[i]=T} else {Motif_GeneDF$jaspar2018[i]=F}
  }

save(Motif_GeneDF,file="/Users/alisonpaquette/TFontology/Motif_DF.RData")

####################################################
####### Part 2:GO BP Transcription Factors##########
####################################################
go.term <- "GO:0003700"
columns(GO.db)  #  "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"
keytypes(GO.db) # "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"
tbl.go <- select(GO.db, keys=go.term, keytype="GOID", columns=columns(GO.db))
as.data.frame(t(select(GO.db, keys=go.term, keytype="GOID", columns=columns(GO.db))))

grep("GO", keytypes(org.Hs.eg.db), v=TRUE) # [1] "GO"    "GOALL"
unique( select(org.Hs.eg.db, "GO:0003700", "SYMBOL", "GOALL")$SYMBOL)
tbl.go <- select(org.Hs.eg.db, keys=go.term, columns="SYMBOL", keytype="GOALL")
TFs_GO<-as.character(unique(tbl.go$SYMBOL))
length(TFs_GO)
#Number of Feasible TFs based on GO_BP & Create Table

#Evidence Code Guidance
#http://www.geneontology.org/page/guide-go-evidence-codes
EvidenceCodes<-c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP","ISS","ISO","ISA","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","IC","ND","IEA")


#Make Gene Specific Data Frame for Evidence Codes
GO_GeneDF<-data.frame(matrix(NA, nrow = length(TFs_GO), ncol = 27))
colnames(GO_GeneDF)<-c("GeneSymbol",EvidenceCodes)
rownames(GO_GeneDF)<-TFs_GO
GO_GeneDF$GeneSymbol<-TFs_GO

#Test to see if Motifs from individual Databases are in Gene DF
for(j in 1:length(EvidenceCodes)){
for(i in 1:length(TFs_GO)){
  temp<-subset(tbl.go,SYMBOL==as.character(GO_GeneDF$GeneSymbol[i]))
  if(length(grep(EvidenceCodes[j],temp$EVIDENCEALL)>0)){GO_GeneDF[j+1]=T} else {GO_GeneDF[i,j+1]=F}
}
}

save(GO_GeneDF,file="/Users/alisonpaquette/TFontology/GO_EvidenceDF.RData")
#######################################################
###########Create Master File & Analyze ###############
#######################################################
library(dplyr)
AllTFs<-merge(GO_GeneDF,Motif_GeneDF,by="GeneSymbol",all=T)
rownames(AllTFs)<-AllTFs$GeneSymbol

#######################################################
########### Evaluate all TFS       ###############
#######################################################
length(intersect(GO_GeneDF$GeneSymbol,Motif_GeneDF$GeneSymbol))

length(GO_GeneDF$GeneSymbol)-length(intersect(GO_GeneDF$GeneSymbol,Motif_GeneDF$GeneSymbol))

length(Motif_GeneDF$GeneSymbol)-length(intersect(GO_GeneDF$GeneSymbol,Motif_GeneDF$GeneSymbol))

A<-setdiff(Motif_GeneDF$GeneSymbol,GO_GeneDF$GeneSymbol)
B<-intersect(Motif_GeneDF$GeneSymbol,GO_GeneDF$GeneSymbol)
C<-setdiff(GO_GeneDF$GeneSymbol, Motif_GeneDF$GeneSymbol)

A.TFs<-AllTFs[A,]
B.TFs<-AllTFs[B,]
C.TFs<-AllTFs[C,]


par(mfrow=c(1,1), mai = c(0.75, 2, 0.75, 0.75))
##Group A
A_Table<-data.frame(matrix(NA, nrow = (length(A.TFs)), ncol = 2))
for(i in 2:length(A.TFs)){
  A_Table[i,2]<-sum(A.TFs[,i])
}
rownames(A_Table)<-colnames(A.TFs)
A_Table[,1]<-colnames(A.TFs)
A_Table[A_Table==0] <- NA
A_Table<-na.omit(A_Table)
colors<-c(rep("dodgerblue4",3))
barplot(A_Table$X2,xlim=c(0,300),col=colors,las=2,names.arg=rownames(A_Table),horiz=T)
abline(v=length(rownames(A.TFs)),col="red",lwd=2,lty=2)

#Group B
B_Table<-data.frame(matrix(NA, nrow = (length(B.TFs)), ncol = 2))
for(i in 2:length(B.TFs)){
B_Table[i,2]<-sum(B.TFs[,i])
}
rownames(B_Table)<-colnames(B.TFs)
B_Table[,1]<-colnames(B.TFs)
B_Table[B_Table==0] <- NA
B_Table<-na.omit(B_Table)
colors<-c(rep("dodgerblue",11),rep("dodgerblue4",3))
barplot(B_Table$X2,xlim=c(0,800),col=colors,las=2,names.arg=rownames(B_Table),horiz=T)
abline(v=length(rownames(B.TFs)),col="red",lwd=2,lty=2)

#Group C
C_Table<-data.frame(matrix(NA, nrow = (length(C.TFs)), ncol = 2))
for(i in 2:length(C.TFs)){
  C_Table[i,2]<-sum(C.TFs[,i])
}
rownames(C_Table)<-colnames(C.TFs)
C_Table[,1]<-colnames(C.TFs)
C_Table[C_Table==0] <- NA
C_Table<-na.omit(C_Table)
colors<-c(rep("dodgerblue",11),rep("dodgerblue4",3))
barplot(C_Table$X2,xlim=c(0,1200),col=colors,las=2,names.arg=rownames(C_Table),horiz=T)
abline(v=length(rownames(C.TFs)),col="red",lwd=2,lty=2)
