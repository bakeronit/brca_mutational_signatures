# functions from scarHRD (https://github.com/sztup/scarHRD/tree/master/R)

preprocess.hrd<-function(seg){
  seg <- seg[!seg[,2] %in% c(paste('chr',c('X','Y','x','y',23,24),sep=''),c('X','Y','x','y',23,24)),]
  seg[,1] <- as.character(seg[,1])
  
  if(! all(seg[,8] <= seg[,7]) ){
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  seg <- shrink.seg.ai.wrapper(seg)
  return(seg)
  
}

shrink.seg.ai<-function(chr.seg){
  new.chr <- chr.seg
  if(nrow(chr.seg) > 1){
    new.chr <- matrix(0,0,ncol(chr.seg))
    colnames(new.chr) <- colnames(chr.seg)
    new.chr <- chr.seg
    seg.class <- c(1)
    for(j in 2:nrow(new.chr)){
      seg_test <- new.chr[(j-1),7] == new.chr[j,7] & new.chr[(j-1),8] == new.chr[j,8]
      if(seg_test){
        seg.class <- c(seg.class, seg.class[j-1])
      }
      if(!seg_test){
        seg.class <- c(seg.class, seg.class[j-1]+1)
      }
    }
    for(j in unique(seg.class)){
      new.chr[seg.class %in% j,4] <- max(new.chr[seg.class %in% j,4])
      new.chr[seg.class %in% j,5] <- sum(new.chr[seg.class %in% j,5])
    }
    new.chr<- new.chr[!duplicated(seg.class),]
  }
  if(nrow(chr.seg) == 1){
    new.chr <- chr.seg
  }
  return(new.chr)
}

shrink.seg.ai.wrapper<-function(seg){
  new.seg <- seg[1,]
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
      }
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}

calc.hrd<-function(seg, nA=7, return.loc=FALSE,sizelimit1){
  nB <- nA+1
  output <- rep(0, length(unique(seg[,1])))
  names(output) <- unique(seg[,1])
  if(return.loc) {
    out.seg <- matrix(0,0,9)
    colnames(out.seg) <- c(colnames(seg)[1:8],'HRD breakpoint')
  }
  #For multiple patients
  for(i in unique(seg[,1])){
    segSamp <- seg[seg[,1] %in% i,]
    chrDel <-vector()
    for(j in unique(segSamp[,2])){
      if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
        chrDel <- c(chrDel, j)
      }
    }
    segSamp[segSamp[,nA] > 1,nA] <- 1
    segSamp <- shrink.seg.ai.wrapper(segSamp)
    segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
    segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,,drop=F]
    segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
    output[i] <- nrow(segLOH)
    if(return.loc){
      if(nrow(segLOH) < 1){next}
      segLOH <- cbind(segLOH[,1:8], 1)
      colnames(segLOH)[9] <- 'HRD breakpoint'
      out.seg <- rbind(out.seg, segLOH)
    }
  }
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
}


calc.ai_new<-function(seg, chrominfo, min.size=1e6, cont = 0,ploidyByChromosome=TRUE, shrink=TRUE){
  seg <- seg[seg[,4]- seg[,3] >= min.size,]
  seg <- seg[seg[,10] >= cont,]
  if(shrink){
    seg <- shrink.seg.ai.wrapper(seg)
  }
  AI <- rep(NA, nrow(seg))
  seg <- cbind(seg, AI)
  samples <- as.character(unique(seg[,1]))
  
  ascat.ploidy <- setNames(seg[!duplicated(seg[,1]),9], seg[!duplicated(seg[,1]),1])
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    if(!ploidyByChromosome){
      ploidy <- vector()
      for(k in unique(sample.seg[,6])){
        tmp <- sample.seg[sample.seg[,6] %in% k,]
        ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
      }
      ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
      sample.seg[,9] <- ploidy
      if(ploidy %in% c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] == sample.seg[,8], c('TRUE', 'FALSE'))]
      }
      if(!ploidy %in%  c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] + sample.seg[,8] == ploidy & sample.seg[,7] != ploidy, c('TRUE', 'FALSE'))]
      }
    }
    new.sample.seg<- sample.seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) == 0){ next}
      if(ploidyByChromosome){
        ploidy <- vector()
        for(k in unique(sample.seg[,6])){
          tmp <- sample.chrom.seg[sample.chrom.seg[,6] %in% k,]
          ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
          ploidy <- ploidy[!names(ploidy) %in% 0]
        }
        ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
        sample.chrom.seg[,9] <- ploidy # update "ploidy" column, so the new calculated value can be returned
        if(ploidy %in% c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] == sample.chrom.seg[,8], c('TRUE', 'FALSE'))]
        }
        if(!ploidy %in%  c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] + sample.chrom.seg[,8] == ploidy & sample.chrom.seg[,8] != 0, c('TRUE', 'FALSE'))]
        }
        sample.seg[sample.seg[,2] %in% i,9] <-ploidy
        sample.seg[sample.seg[,2] %in% i,'AI'] <-sample.chrom.seg[,'AI']
      }
      
      if(class(chrominfo) != 'logical'){# Here we consider the centromere
        if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[1,4] < (chrominfo[i,2])){
          sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
        }
        if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[nrow(sample.chrom.seg),3] > (chrominfo[i,3])){
          sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
        }
      }
      if(nrow(sample.seg[sample.seg[,2]==i,]) == 1 & sample.seg[sample.seg[,2]==i,'AI'][1] != 0){
        sample.seg[sample.seg[,2]==i,'AI'][1] <- 3
      }
    }
    seg[seg[,1] %in% j,] <- sample.seg
  }
  samples <- as.character(unique(seg[,1]))
  #0 = no AI, 1=telomeric AI, 2=interstitial AI, 3= whole chromosome AI
  no.events <- matrix(0, nrow=length(samples), ncol=12)
  rownames(no.events) <- samples
  colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", "Telomeric LOH",  "Mean size", "Interstitial LOH", "Mean Size", "Whole chr LOH", "Ploidy", "Aberrant cell fraction")
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    no.events[j,1] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,2] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,3] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,4] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,5] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
    no.events[j,11] <- ascat.ploidy[j]
    no.events[j,12] <- unique(sample.seg[,10]) # aberrant cell fraction
    #Here we restrict ourselves to real LOH
    sample.seg <- sample.seg[sample.seg[,8] == 0,]
    no.events[j,6] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,7] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,8] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,9] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,10] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
  }
  return(no.events)
}


calc.lst<-function(seg, chrominfo=chrominfo,nA=7,chr.arm='no'){
  nB <- nA+1
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- c()
    chroms <- unique(sample.seg[,2])
    chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]
    for(i in chroms){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(chr.arm !='no'){
        p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
        q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
      }
      if(nrow(sample.chrom.seg) < 2) {next}
      sample.chrom.seg.new <- sample.chrom.seg
      if(chr.arm == 'no'){
        p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,2],,drop=F] # split into chromosome arms
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,3],,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- chrominfo[i,3]
        if(nrow(p.arm) > 0){
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- chrominfo[i,2]
        }
      }
      if(chr.arm != 'no'){
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- q.min
        if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
          p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] # split into chromosome arms
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- p.max
        }
      }
      n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        p.arm <- p.arm[-(n.3mb[1]),]
        p.arm <- shrink.seg.ai(p.arm)
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      }
      if(nrow(p.arm) >= 2){
        p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(p.arm)){
          if(p.arm[k,9] == 1 & p.arm[(k-1),9]==1 & (p.arm[k,3]-p.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
          }
        }
      }
      n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        q.arm <- q.arm[-(n.3mb[1]),]
        q.arm <- shrink.seg.ai(q.arm)
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      }
      if(nrow(q.arm) >= 2){
        q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(q.arm)){
          if(q.arm[k,9] == 1 & q.arm[(k-1),9]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
            
          }
        }
      }
    }
    output[j] <- sum(sample.lst)
  }
  return(output)
}

scar_score<-function(seg,reference = "grch38", chr.in.names=TRUE, m,seqz=FALSE, ploidy=NULL, sizelimitLOH=15e6, outputdir=NULL){
  
  if (is.null(outputdir)){
    outputdir=getwd()
  }
  if (reference == "grch38"){
    chrominfo = chrominfo_grch38
    if(!chr.in.names){
      chrominfo$chrom = gsub("chr","",chrominfo$chr)
      rownames(chrominfo) = chrominfo$chr
    }
  } else if(reference == "grch37"){
    chrominfo = chrominfo_grch37
    if(!chr.in.names){
      chrominfo$chrom = gsub("chr","",chrominfo$chr)
      rownames(chrominfo) = chrominfo$chr
    }
  } else if(reference == "mouse"){
    chrominfo = chrominfo_mouse
    if(!chr.in.names){
      chrominfo$chrom = gsub("chr","",chrominfo$chr)
      rownames(chrominfo) = chrominfo$chr
    }
  } else {
    stop()
  }
  

  #seg<-read.table(seg,header=T, check.names = F, stringsAsFactors = F, sep="\t")
  
  seg[,9]<-seg[,8]
  seg[,8]<-seg[,7]
  seg[,7]<-seg[,6]
  seg[,10]<-rep(1,dim(seg)[1])
  
    
  #prep
  cat('Determining HRD-LOH, LST, TAI \n')
  seg<-preprocess.hrd(seg)
  #Calculating the hrd score:
  res_hrd <- calc.hrd(seg,sizelimit1=sizelimitLOH)
  #Calculating the telomeric allelic imbalance score:
  res_ai<- calc.ai_new(seg = seg, chrominfo = chrominfo) #<-- the first column is what I need
  #Calculating the large scale transition scores:
  res_lst <- calc.lst(seg = seg, chrominfo = chrominfo) #<-- You need to use the chrominfo.snp6 file! Nicolai sent it to you!
  sum_HRD0<-res_lst+res_hrd+res_ai[1]
  
  if (is.null(ploidy)){
    sum_HRDc<-NA
  } else {
    sum_HRDc<-res_lst-15.5*ploidy+res_hrd+res_ai[1]
  }
  
  #HRDresulst<-c(res_hrd,res_ai,res_lst,sum_HRD0,sum_HRDc)
  #names(HRDresulst)<-c("HRD",colnames(res_ai),"LST", "HRD-sum","adjusted-HRDsum")
  HRDresulst<-c(res_hrd,res_ai[1],res_lst,sum_HRD0)
  names(HRDresulst)<-c("HRD",colnames(res_ai)[1],"LST", "HRD-sum")
  run_name<-names(sum_HRD0)
  #write.table(t(HRDresulst),paste0(outputdir,"/",run_name,"_HRDsumresults.txt"),col.names=NA,sep="\t",row.names=unique(seg[,1]))
  return(t(HRDresulst))
}