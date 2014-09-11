
##' @title Regional Plot Association Results of Pig
##' 
##' @description This is a plotting tool to simultaneously display the association p-value,
##'              the LD, recombination rate and annotated gene in region of interesting in pig.
##'              The interpretation and visualization tool of association results will
##'              contribute to further identification of causal mutations.
##'
##' @details nothing
##'  
##' @param pval A R object.the input data of GWAS result constructed by at least 2 columns,
##'             P-value and SNP identity in Illumina system. The two colnames of pval must be
##'             "MkName" and "Pval".
##' @param plotChr a character:the chromosome to plot. such as "1", "2", "X" and "Y".
##' @param plotRegion a numeric vector containing the start position and the end position of plot. 
##'                   The unit is "Mb".
##' @param geneRegion a numeric vector containing the start position and the end position of gene structure in bottom panel.
##'                   The unit is "Mb".
##' @param isLD a logical value indicating whether plots the LD information
##' @param LDresult the data of linkage disequilibrium (LD) consisting of the maker name of other SNPs in plotting Region and 
##'                 the R2 associated with the top SNP. The two colnames are "SNP_B" and "R2", respectively.
##' @param geno the corresponding path of genotype data if we want to use own LD pattern by
##'                 inputting a plink formatted genotype. Then, RPARP will automatically call PLINK software (version 1.07) and calculated LD
##' @param isRecomb a logical value indicating whether plots the recombination information
##' @param recomb A dataframe contains the SNPs name and the recombination information, 
##'                of which colnames are "MkName" and "recomb", respectively.
##' @param isAnno a logical value indicating whether plots every SNP annotaion information
##' @param main The main title (on top) using font and size (character expansion) par("font.main") and color par("col.main").
##' @param ref  the name of the reference genome. There are two different genome: 
##'             the widespread commercial breed: Duroc (Sscrofa10.2) or Chinese indigenous breed: Wuzhishan pig (WZSP)
##' @param xlab a title for the x axis.
##' @return a figure
##' @examples png("MCV.png", width=3200, height=2200,res=300) # this is the most appropriate size and resolution
##' plotRegion(MCV_pval)
##' dev.off()

##' @author Feng Zhang and Zhiyan Zhang
 

plotRegion <- function(pval,plotChr,plotRegion,geneRegion,isLD=TRUE,LDresult,geno,isRecomb=TRUE,recomb,isAnno=TRUE,main="",ref='Wzs',xlab) {

  # mapping P-value to reference genome
  pval <- merge(get(paste(ref,"Mk",sep="")),pval,by="MkName")
  pval <- pval[pval$Chromosome!=0,]
  
  if(missing(plotChr)) plotChr <- pval[which.min(pval$Pval),]$Chromosome
  pval <- pval[pval$Chromosome==plotChr,]
  pval <- pval[order(pval$Position),]
  HstPvalIdx <- which.min(pval$Pval)
  HstSNP <- pval$MkName[HstPvalIdx]
  
  #select the region to plot
  if(missing(plotRegion))    plotRegion <- plot_region(pval,HstPvalIdx)
  
  pval$Position <- pval$Position/1000000
  pval <- pval[pval$Position>=min(plotRegion)& pval$Position <= max(plotRegion),]
  
  #merge linkage disequilibrium result with other SNPs
  if(isLD){
    if(missing(LDresult)) {
      if(missing(geno)){
        LDresult <- MeanLD[MeanLD$SNP_A==HstSNP,][,c('SNP_B','R2')]
        pval <- merge(pval,LDresult,by.x='MkName',by.y='SNP_B',all.x=T)
      }else{
        ind <- which(colnames(geno) %in% pval$MkName)
        tmp <- LD(makeGenotypes(geno[,ind]))$R^2
        LDresult <- data.frame(R2=tmp[rownames(tmp)==HstSNP,],
                               SNP_B=names(tmp[rownames(tmp)==HstSNP,]))
        pval<- merge(pval,LDresult,by.x='MkName',by.y='SNP_B',all.x=T)
      }
    }else{
      pval<- merge(pval,LDresult,by.x='MkName',by.y='SNP_B',all.x=T)
    }
  }
  


  #merge recombination with other SNPs
  if(isRecomb){
    if(missing(recomb)){
      pval <- merge(pval,SutaiRecomb,by='MkName',all.x=T)
    }else{
      if(recomb=='F2Recomb'){
        pval<- merge(pval,F2Recomb,by='MkName',all.x=T)
      }else{
        pval <- merge(pval,recomb,by='MkName',all.x=T)
      }
    }
  }
 
  
  #set attribution of pval
  pval <- pval_set(pval,HstSNP,isAnno)
  
  #draw the regional figure
  layout(matrix(c(1, 2), 2), 1, c(4, 1))
  layout.show(2)
  plot.first(pval,main,isAnno,isRecomb)  

  #merge the gene annotation 
  pval <- pval[order(pval$Position),] 
  gff <- get(paste(ref,'Gff',sep=''))
  if(missing(geneRegion)) geneRegion <- gene_region(pval,HstSNP)
  gff <- gff[gff$chr==plotChr & gff$leftPos>=geneRegion[1]*1000000 & gff$rightPos <= geneRegion[2]*1000000,]  #set gene structual attribution in allGene
  allGene <- allGene_set(gff)
  if (missing(xlab)) xlab=paste(paste("SSC",plotChr,sep = ""), "(Mb)",sep = " ")
  plot.second(allGene,gff,geneRegion,xlab) #draw gene structure 
}



gene_region <- function(pval,HstSNP){
  #StartIdx <- which(pval$Position < (pval[pval$MkName==HstSNP,'Position']-0.2))
  #EndIdx <- which(pval$Position > (pval[pval$MkName==HstSNP,'Position']+0.2))
  geneStart <- pval[pval$MkName==HstSNP,'Position']-0.2
  geneEnd <- pval[pval$MkName==HstSNP,'Position']+0.2
  if(geneStart <= 0) geneStart <- min(pval$Position)
  if(geneEnd>max(pval$Position)) geneEnd <- max(pval$Position)
  c(geneStart,geneEnd)
}
plot_region <- function(pval,HstPvalIdx){
  CI.left <- which(pval[1:(HstPvalIdx-1),]$Pval/pval[HstPvalIdx,]$Pval <100) #select left snp threshold of which pvalue 1000 smaller than the top snp
  if(length(CI.left)>0) {CI.left <- pval$Position[min(CI.left)-1]
  } else {
    cat("all left hand SNPs' pvalue 100 times big than the top SNPs, extend to 1000 kb of the left hand\n") 
    indCI <- max(which(pval$Position < (pval[HstPvalIdx,'Position']-1000000)))
    CI.left <- pval$Position[indCI]
  }
  
  CI.Right <- which(pval[(HstPvalIdx+1):nrow(pval),]$Pval/pval[HstPvalIdx,]$Pval <100)
  if(length(CI.Right)>0) {CI.Right <- pval$Position[max(CI.Right)+HstPvalIdx+1]
  } else {
    cat("all Right hand SNPs' pvalue 100 times big than the top SNPs, extend to 1000 kb of the Right hand\n") 
    indCI <- min(which(pval$Position > (pval[HstPvalIdx,'Position']+1000000)))
    CI.Right <- pval$Position[indCI]
  }  
  
  c(CI.left/1000000,CI.Right/1000000)
  #pval <<- pval[CI.left:CI.Right,]
}


pval_set <- function(pval,HstSNP,isAnno){
  #set color
  pval$color <- "black"
  pval$color[pval$R2 >= 0.8] = "red"
  pval$color[pval$MkName==HstSNP]="red"
  pval$color[pval$R2 < 0.8 & pval$R2 >= 0.6] = "blue" 
  pval$color[pval$R2 < 0.6 & pval$R2 >= 0.4] = "orange"
  pval$color[pval$R2 < 0.4 & pval$R2 >= 0.2] = "purple" 
  
  #set cex
   pval$cex <- 1
# #   pval$cex[pval$R2 >= 0.8] = 1.25
# #   pval$cex[pval$R2 < 0.8 & pval$R2 >= 0.5] = 1.25
# #   pval$cex[pval$R2 < 0.5 & pval$R2 >= 0.2] = 1.25
   pval$cex[pval$MkName==HstSNP] <-2
  
  #set pch
  pval$pch=21
  if(isAnno){
    pval$pch[pval$Anno=='non-synonymous'] <- 25
    pval$pch[pval$Anno=='synonymous'] <- 24
    pval$pch[pval$Anno=="3'UTR"|pval$Anno=="5'UTR"]<- 23
    pval$pch[pval$Anno=='Intron'] <- 22
  }
  
  pval
}

plot.first <- function(pval,main,isAnno,isRecomb){
  pval <- pval[order(pval$Position),]
  par(mar=c(1,4,1,4),mgp=c(2,0.5,0))
  plot(pval$Position,-log10(pval$Pval),lwd = 2, col = as.character(pval$color),bg= as.character(pval$color),type = "p", pch=pval$pch,cex = pval$cex,
       ylab = expression(-log[10](p)),xlab ='' ,
       main = main,  tck=-0.01,
       ylim = range(0, ceiling(max(-log10(pval$Pval)))), cex.lab=1.25)
 
  legend("topright",legend=c("0.8~1","0.6~0.8","0.4~0.6","0.2~0.4","0~0.2"), 
         pch=15,pt.cex=2,col=c("red","blue","orange","purple","black"),
         title=expression(r^2),title.adj=0.2, 
         bty="o",cex=0.8,y.intersp=0.9)
  if(isAnno){
    legend("topleft",legend=c('non-synony','synonym','UTR','Intron','InterGe/Un'), col = "black",
           bty="o",cex=0.8,pt.cex=1,pch=c(25,24,23,22,21)) 
  }

  if(isRecomb==TRUE){
    par(new=TRUE)
    pval_RF <- pval[,c("Position","recomb")]
    pval_RF <- rbind(data.frame(Position=pval_RF$Position-0.00001,recomb=0),pval_RF,
                     data.frame(Position=pval_RF$Position+0.00001,recomb=0))
    pval_RF <- pval_RF[order(pval_RF$Position),]
    pval_RF$recomb[is.na(pval_RF$recomb)] <- 0
    plot(pval_RF$Position,pval_RF$recomb*100,type='l',lty=1,xaxt="n",yaxt="n",xlab="",ylab='')
    axis(4)
    mtext("recombination ratio (%)", side = 4, cex=1.25,line=2)
  }
}		   

allGene_set <- function(gff){
  allGene <- gff[gff$structure=='mRNA',]
   
  #set y axis position
  allGene$strand_ypos <-  rep(c(0.15,-0.15),len=nrow(allGene))
  allGene$gene_ypos <- allGene$strand_ypos+rep(c(0.06,0.06,0.13,0.13),len=nrow(allGene))
  
  
  
  #set arrow
  allGene$geneLabel <- paste('<--',allGene$geneName,sep=' ')
  allGene$geneLabel[allGene$strand=='+'] <- paste(allGene$geneName[allGene$strand=='+'],'-->',sep=' ')
  
  allGene$leftPos <- allGene$leftPos/1000000
  allGene$rightPos <- allGene$rightPos/1000000
  allGene
}

plot.second <- function(allGene,gff,geneRegion,xlab){
  par(mar=c(3,4,0.5,4),mgp=c(1.5,0.5,0))
  plot(0,type='n',xlim=c(geneRegion[1], geneRegion[2]),ylim=c(-0.19,0.33),
       xlab=xlab,ylab='',yaxt='n')
  box()
  segments(allGene$leftPos ,allGene$strand_ypos ,allGene$rightPos ,allGene$strand_ypos,col='black',lwd=2)
  midPos <- (allGene$rightPos +allGene$leftPos )/2  
  gff$leftPos <- gff$leftPos/1000000
  gff$rightPos <- gff$rightPos/1000000
  for(x in 1:nrow(allGene)){
    geneGff <- gff[gff$leftPos>=allGene$leftPos[x] & gff$rightPos <= allGene$rightPos[x] & gff$structure != 'mRNA',]
    #set UTR colours
    geneGff$col <- 'black'
    geneGff$col[grep('UTR',geneGff$structure)] <- 'red'
    
    rect(xleft=geneGff[,'leftPos'],ybottom= allGene$strand_ypos[x]-0.02,
         xright=geneGff[,'rightPos'],ytop=allGene$strand_ypos[x]+0.02,col=geneGff$col,border=geneGff$col)
    text( midPos[x], y=allGene$gene_ypos[x],labels=allGene$geneLabel[x],cex=0.7,font=2)
  }
}
  
  
