
#dependency

library(dplyr,quietly = T,warn.conflicts = T)
library(ggplot2,quietly = T)


mytime <- function(expression) {
  t1 <- proc.time()
  res <- eval(expression)
  print(proc.time() - t1)
  res
}

runparallel <- function(objlist, FUN, cores = 7, ...) {
  
  cl <- parallel::makeCluster(cores)
  res <- parallel::parLapply(cl, X = objlist, fun = FUN, ...)
  parallel::stopCluster(cl)
  
  invisible(res)
}

# from hyb file to contact matrix core function
hyb2CM <- function(hyblines, binsize = 10) {
  library(dplyr)
  hyblines <- floor(hyblines / binsize) # as.integer only supports vector
  Arm5=(0:3000)*10
  Arm3=(0:3000)*10
  matrix <- matrix(data = 0,nrow = 3001,ncol = 3001,dimnames=list(Arm5,Arm3))
  for(i in 1:nrow(hyblines)){
    matrix[hyblines[i,1]:hyblines[i,2],hyblines[i,3]:hyblines[i,4]] <- matrix[hyblines[i,1]:hyblines[i,2],hyblines[i,3]:hyblines[i,4]] +1
  }
  data <- reshape2::melt(matrix, varnames = c("Arm5", "Arm3"),value.name="n") %>% dplyr::filter(n>0)
  
}


# frome hyb to contact matrix file function
HybToMatrix <- function(Hyb,filename,...){
  cores <- parallel::detectCores() -1
  idxes <- sort(sample(letters[1:cores], nrow(Hyb), replace = TRUE))
  HybList <- split(Hyb, idxes)
  result <- mytime({
    runparallel(HybList, FUN = hyb2CM, cores = cores, binsize = 10)
  })%>% bind_rows()  %>% group_by(Arm5,Arm3) %>% summarise(n = sum(n)) %>% as.data.frame()
  data.table::fwrite(result,file = filename,row.names = F,sep = "\t")
 # return(result)
}

call_intra_matrix <- function(hybfile,removeOverlap=F,Bait_RNA="NC_045512-2_SARS-CoV-2_virusRNA",matrxfile){
  hyb <- data.table::fread(hybfile,sep = "\t",select = 1:15)
  library(dplyr)
  intra.hyb <- hyb %>% filter(V4 == V10) %>% filter(V4 == Bait_RNA) %>% select(c(7,8,13,14))
  if(removeOverlap){
    intra.hyb <-  intra.hyb %>% mutate(overlap = 1+pmin(V8,V14)-pmax(V7,V13)) %>% dplyr::filter(overlap < 0)
  }
  # intra 
  Contact_matrix <- HybToMatrix(Hyb = intra.hyb,filename = matrxfile)
}



class_stats_hyb <- function(hybfile,Bait_RNA="NC_045512-2_SARS-CoV-2_virusRNA",breaks,...){
  library(dplyr)
  library(ggplot2)
  hyb <- data.table::fread(hybfile,sep = "\t",select = 1:15)
  intra.hyb <- hyb %>% filter(V4 == V10) %>% filter(V4 == Bait_RNA) %>% select(c(7,8,13,14)) %>% mutate(overlap = 1+pmin(V8,V14)-pmax(V7,V13)) 
  overlap_hist <- intra.hyb %>% count(overlap)
  overlap_hist <- overlap_hist %>% mutate(class = cut(overlap,breaks = breaks,labels = c("long range","distal","proximal","overlap"))) %>% group_by(class) %>% summarise(n=sum(n))
  
}







stats_range_oritation_hyb <- function(hybfile,Bait_RNA="NC_045512-2_SARS-CoV-2_virusRNA",breaks=c(-30000,-1000,-100,0,200),labels=c("long range","distal","proximal","overlap"),...){
  library(dplyr)
  library(ggplot2)
  hyb <- data.table::fread(hybfile,sep = "\t",select = 1:15)

  intra.hyb <- hyb %>% filter(V4 == V10) %>% filter(V4 == Bait_RNA) %>% select(c(7,8,13,14)) %>% mutate(overlap = 1+pmin(V8,V14)-pmax(V7,V13)) 
 # intra.hyb <- hyb %>% filter(V4 == V10) %>% filter(V4 == Bait_RNA) %>% mutate(overlap = 1+pmin(V8,V14)-pmax(V7,V13)) 
  intra.hyb <- intra.hyb %>% mutate(orientation= ifelse(overlap > 0,"self",ifelse(V7< V13, "5'-3'","3'-5'")))
  
  intra.all <- nrow(intra.hyb)
  Count_orientation <- intra.hyb %>% count(orientation)
  overlap_hist <- intra.hyb %>% mutate(class = cut(overlap,breaks = breaks,labels = labels)) %>% group_by(class,orientation) %>% summarise(n=n())
}


# for a given samplelist, calculate distance and orientation distribution, then ggplot 
stat_distance_orientation <- function(samplelist,hybdir="../hyb_R2/hyb/",suffix= "_comp_SARS-CoV-2_no_polyA_hybrids_ua.hyb.gz",breaks = breaks,labels = labels,strait=T,...){
  
  stats_range_df <- lapply(names(samplelist), function(s){
    
    x <- samplelist[[s]]
    hybfile <- paste0(hybdir,x,suffix)
    
    stats_range_oritation_x <- stats_range_oritation_hyb2(hybfile = hybfile,Bait_RNA = "NC_045512-2_SARS-CoV-2_virusRNA",breaks = breaks,labels = labels) %>% as.data.frame() %>% mutate(sample=s)
  }) %>% bind_rows() %>% as.data.frame()
    maxcount <- stats_range_df %>% group_by(class,sample) %>% summarise(n = sum(n)) 
    maxcount <- max(maxcount$n)
  
  p_list <- lapply(names(samplelist), function(x){
    gdata <- stats_range_df %>% filter(sample == x)
    
    p <- ggplot(gdata,aes(class,n))+
      geom_bar(aes(fill=orientation),stat = "identity")+
      scale_fill_manual(values =col_orientation )+
    #  scale_y_continuous(limits = c(0,maxcount))+
      labs(title = x,y="chimeric reads counts")+
      theme_classic()+theme(axis.text.x  = element_text(angle = 60,hjust = 1))
    
    if(strait){
      p<- p+scale_y_continuous(limits = c(0,maxcount))
    }
    return(p)
  })
  
  
}



# export contact matrix to java treeview cdt files
# you need total mapped count for each sample/group

CPM2cdt <- function(samplelist,indir="matrix_all/",matxsufix="_R2_all.matrix.gz",outdir="cdt/",range=c(0,500),cutoff=0,cdtFileName,contrast=1,readcounts,...){
  
  # read raw matrix and total mapped readcount to make CPM matrix
  # then export to cdt
  #readcounts is a named vector, recording total read count of the sample
  suppressWarnings(library(ggplot2,quietly = T,warn.conflicts = F))
  suppressWarnings(library(dplyr,quietly = T,warn.conflicts = F))
  
  levels = 10*(floor(range[2]/10) : floor(range[1]/10))
  
  data <- lapply(samplelist, function(x){
    
    data <- data.table::fread(paste0(indir,x,matxsufix),sep = "\t",col.names = c("Arm5","Arm3","n"))
    data <- data %>% filter((Arm5 <= range[2] & Arm5 >= range[1]) &(Arm3 <= range[2] & Arm3 >= range[1]))
    data <- data %>% filter(n > cutoff)  %>% mutate(n=1e6*n/readcounts[x])
    
    matrix2cdt(matrix = data,levels =levels,contrast = contrast,cdtFileName = paste0(outdir,x,cdtFileName) )
    
  }) %>% bind_rows () %>% as.data.frame()
  
}




plot_local_mat <- function(samplelist,indir="matrix_rm_SL/",matxsufix=".merged_R2_remove_selfligation.matrix.gz",range=c(0,500),cutoff=0,ceiling=0.95,annotation=NULL,...){
  suppressWarnings(library(ggplot2,quietly = T,warn.conflicts = F))
  suppressWarnings(library(dplyr,quietly = T,warn.conflicts = F))
  data <- lapply(samplelist, function(x){

    data <- data.table::fread(paste0(indir,x,matxsufix),sep = "\t",col.names = c("Arm5","Arm3","n"))
    data <- data %>% filter((Arm5 <= range[2] & Arm5 >= range[1]) &(Arm3 <= range[2] & Arm3 >= range[1]))
    data <- data %>% filter(n > cutoff) %>% mutate(sample=x)
  }) %>% bind_rows () %>% as.data.frame()
  
  
  data$n[data$n > quantile(data$n,ceiling)] <- quantile(data$n,ceiling)
  
  maxpoint <- max(data$n)
  regionMax <- max(data$Arm5,data$Arm3)
  regionMin <- min(data$Arm5,data$Arm3)
  
  p.list <- lapply(samplelist,function(x){
    gdata <- data %>% filter(sample == x)
    p <- ggplot(gdata,aes(Arm5,Arm3))+
      geom_rect(aes(fill=n,xmin=Arm5,xmax=Arm5+10,ymin=Arm3,ymax=Arm3+10))+
      annotate(geom = "rect",xmin = range[1]-5, xmax = range[2]-5, ymin = range[1]-5, ymax = range[2]-5,fill=NA,color="black")+ 
      
      annotate(geom = "rect",xmin = floor(annotation$V2/10)*10, xmax = floor(1+annotation$V3/10)*10, ymin = floor(annotation$V4/10)*10, ymax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      annotate(geom = "rect",ymin = floor(annotation$V2/10)*10, ymax = floor(1+annotation$V3/10)*10, xmin = floor(annotation$V4/10)*10, xmax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      scale_fill_gradient(low = "#ffffff",high="#000000",limits=c(0,maxpoint))+
      scale_x_continuous(breaks = c(regionMin,regionMax))+
      labs(fill="chimeras",title = x)+
      theme_bw()+
      theme(aspect.ratio = 1,panel.grid = element_blank())
  })
}

#plot local normalized contact matrix
plot_local_CPM <- function(samplelist,indir="matrix_all/",matxsufix=".merged_R2_all.matrix.gz",range=c(0,500),cutoff=0,ceiling=0.95,annotation=NULL,readcounts,...){
  suppressWarnings(library(ggplot2,quietly = T,warn.conflicts = F))
  suppressWarnings(library(dplyr,quietly = T,warn.conflicts = F))
  data <- lapply(samplelist, function(x){
    
    data <- data.table::fread(paste0(indir,x,matxsufix),sep = "\t",col.names = c("Arm5","Arm3","n"))
    data <- data %>% filter((Arm5 <= range[2] & Arm5 >= range[1]) &(Arm3 <= range[2] & Arm3 >= range[1]))
    data <- data %>% filter(n > cutoff) %>% mutate(sample=x) %>% mutate(n=1e6*n/readcounts[x])
  }) %>% bind_rows () %>% as.data.frame()
  
  
  data$n[data$n > quantile(data$n,ceiling)] <- quantile(data$n,ceiling)
  
  maxpoint <- max(data$n)
  regionMax <- max(data$Arm5,data$Arm3)
  regionMin <- min(data$Arm5,data$Arm3)
  
  p.list <- lapply(samplelist,function(x){
    gdata <- data %>% filter(sample == x)
    p <- ggplot(gdata,aes(Arm5,Arm3))+
      geom_rect(aes(fill=n,xmin=Arm5,xmax=Arm5+10,ymin=Arm3,ymax=Arm3+10))+
      annotate(geom = "rect",xmin = range[1], xmax = range[2]+10, ymin = range[1], ymax = range[2]+10,fill=NA,color="black")+ 
      annotate(geom = "rect",xmin = floor(annotation$V2/10)*10, xmax = floor(1+annotation$V3/10)*10, ymin = floor(annotation$V4/10)*10, ymax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      annotate(geom = "rect",ymin = floor(annotation$V2/10)*10, ymax = floor(1+annotation$V3/10)*10, xmin = floor(annotation$V4/10)*10, xmax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      annotate(geom = "text",x=(annotation$V2+annotation$V3)/2,y=(annotation$V4+annotation$V5)/2,label=annotation$V1,color="#0071bc",size=2)+
      scale_fill_gradient(low = "#ffffff",high="#000000",limits=c(0,maxpoint))+
      scale_x_continuous(breaks = c(regionMin,regionMax),name = NULL)+
      labs(fill="chimeras per million reads",title = x)+
      theme_bw()+
      theme(aspect.ratio = 1,panel.grid = element_blank(),panel.border =  element_blank())
  })
}

# plot interaction of two regions from contact matrix
plot_2region_mat <- function(samplelist,indir="matrix_rm_SL/",matxsufix=".merged_R2_remove_selfligation.matrix.gz",symetry=T,range1=c(0,500),range2=c(29500,29870),cutoff=0,ceiling=0.95,annotation=NULL,...){
  suppressWarnings(library(ggplot2,quietly = T,warn.conflicts = F))
  suppressWarnings(library(dplyr,quietly = T,warn.conflicts = F))
  data <- lapply(samplelist, function(x){
    
    data <- data.table::fread(paste0(indir,x,matxsufix),sep = "\t",col.names = c("Arm5","Arm3","n"))
    data.1 <- data %>% filter((Arm5 <= range1[2] & Arm5 >= range1[1]) &(Arm3 <= range2[2] & Arm3 >= range2[1]))
    data.2 <- data %>% filter((Arm5 <= range2[2] & Arm5 >= range2[1]) &(Arm3 <= range1[2] & Arm3 >= range1[1]))
    
    
    data <- rbind(data.1,data.2) %>% filter(n > cutoff) %>% mutate(sample=x)
  }) %>% bind_rows () %>% as.data.frame()
  
  #make symetry
  data.1 <- data.frame(Arm5=data$Arm5,Arm3=data$Arm3,n=data$n,sample=data$sample) 
  data.2 <- data.frame(Arm5=data$Arm3,Arm3=data$Arm5,n=data$n,sample=data$sample)
  data <- rbind(data.1,data.2) %>% mutate(bin1=pmin(Arm5,Arm3),bin2=pmax(Arm5,Arm3)) %>% group_by(bin1,bin2,sample) %>% summarise(n=sum(n)) %>%
    mutate(Arm5=bin1,Arm3=bin2)
  
  data$n[data$n > quantile(data$n,ceiling)] <- quantile(data$n,ceiling)
  
  
  
  maxpoint <- max(data$n)
  regionMax <- max(data$Arm5,data$Arm3)
  regionMin <- min(data$Arm5,data$Arm3)
  
  p.list <- lapply(samplelist,function(x){
    gdata <- data %>% filter(sample == x)
    p <- ggplot(gdata,aes(Arm5,Arm3))+
      geom_rect(aes(fill=n,xmin=Arm5,xmax=Arm5+10,ymin=Arm3,ymax=Arm3+10))+
      annotate(geom = "rect",xmin = range1[1]-5, xmax = range1[2]-5, ymin = range2[1]-5, ymax = range2[2]-5,fill=NA,color="black")+ 
      
      annotate(geom = "rect",xmin = floor(annotation$V2/10)*10, xmax = floor(1+annotation$V3/10)*10, ymin = floor(annotation$V4/10)*10, ymax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      annotate(geom = "rect",ymin = floor(annotation$V2/10)*10, ymax = floor(1+annotation$V3/10)*10, xmin = floor(annotation$V4/10)*10, xmax = floor(1+annotation$V5/10)*10,fill=NA,color="black")+
      scale_fill_gradient(low = "#ffffff",high="#000000",limits=c(0,maxpoint))+
      
      labs(fill="chimeras",title = x)+
      theme_bw()+
      theme(aspect.ratio = 1,panel.grid = element_blank())
  })
}



# contact matrix to zoom-in heatmaps
CPM2cdt_zoom <- function(samplelist,indir="matrix_all/",matxsufix="_R2_all.matrix.gz",outdir="cdt/",range_Arm5=c(0,500),range_Arm3=c(0,500),cutoff=0,cdtFileName,contrast=1,readcounts,...){
  
  # read raw matrix and total mapped readcount to make CPM matrix
  # then export to cdt
  #readcounts is a named vector, recording total read count of the sample
  suppressWarnings(library(ggplot2,quietly = T,warn.conflicts = F))
  suppressWarnings(library(dplyr,quietly = T,warn.conflicts = F))
  
  levels_5 = 10*(floor(range_Arm5[2]/10) : floor(range_Arm5[1]/10))
  levels_3 = 10*(floor(range_Arm3[1]/10) : floor(range_Arm3[2]/10))
  data <- lapply(samplelist, function(x){
    
    data <- data.table::fread(paste0(indir,x,matxsufix),sep = "\t",col.names = c("Arm5","Arm3","n"))
    data <- data %>% filter((Arm5 <= range_Arm5[2] & Arm5 >= range_Arm5[1]) &(Arm3 <= range_Arm3[2] & Arm3 >= range_Arm3[1]))
    matrix <- data %>% filter(n > cutoff)  %>% mutate(n=1e6*n/readcounts[x])
    
    
    
    data.table <- data.table::data.table(bin1=factor(matrix$Arm5,levels = levels_5),
                                         bin2=factor(matrix$Arm3,levels = levels_3),
                                         n   =matrix$n/contrast) 
    
    d <- data.table::dcast.data.table(data =data.table,formula = bin2~bin1,fill = 0,drop = F,value.var="n")
    
    d1 <- cbind(data.frame(COL_1 =d$bin2,COL_2=d$bin2,GWEIGHT=1),d[,2:ncol(d)]) %>% type.convert()
    d2 <- rbind(c("EWEIGHT",rep(1,ncol(d1)-1)),d1)
    
    data.table::fwrite(d2,file=paste0(outdir,x,cdtFileName),sep = "\t",quote = F)    
    
    
    
  })
  
}

# correaltion between two matrix
cor_matrix <- function(x1,x2,range=c(0,30000)){
  mat1 <- data.table::fread(file = x1,sep = "\t",col.names = c("Arm5","Arm3","n")) %>% 
    filter(Arm5 <= range[2] & Arm5 >= range[1] &Arm3 <= range[2] & Arm3 >= range[1]) %>% as.data.frame()
  mat2 <- data.table::fread(file = x2,sep = "\t",col.names = c("Arm5","Arm3","n")) %>%
    filter(Arm5 <= range[2] & Arm5 >= range[1] &Arm3 <= range[2] & Arm3 >= range[1]) %>% as.data.frame()
  levels <- 10*(floor(range[1]/10):floor(range[2]/10) )
  mat1 <- data.table::data.table(Arm5=factor(mat1$Arm5,levels = levels),Arm3=factor(mat1$Arm3,levels = levels),n=mat1$n) %>%
    data.table::dcast.data.table(formula = Arm5~Arm3,fill=0,drop = F,value.var = "n") %>% data.table::melt.data.table(id.vars = "Arm5",variable.name = "Arm3",value.name = "n")
  
  mat2 <- data.table::data.table(Arm5=factor(mat2$Arm5,levels = levels),Arm3=factor(mat2$Arm3,levels = levels),n=mat2$n) %>%
    data.table::dcast.data.table(formula = Arm5~Arm3,fill=0,drop = F,value.var = "n") %>% data.table::melt.data.table(id.vars = "Arm5",variable.name = "Arm3",value.name = "n")
  
  data <- data.frame(mat1=mat1$n,mat2=mat2$n) %>% filter(mat1+mat2 >0)
 # t1 <- cor.test(data$mat1,data$mat2)
  t2 <- cor.test(mat1$n,mat2$n)
  return(t2)
}


# turn COMRADE cluster to bin-bin pair

cluster2CM <- function(sheet = "SARS-CoV-2_cluster_Sample1",binsize = 10){
  Cluster_gRNA_S1 <- read_excel("/Volumes/HiR/SARS2-HiR/publicData/1-s2.0-S1097276520307826-mmc3(clusters).xlsx", 
                                sheet = sheet)
  hyblines <- Cluster_gRNA_S1[,2:5]
  CM <- hyb2CM(hyblines =hyblines, binsize =binsize) %>% mutate(bin1=pmin(Arm5,Arm3),bin2=pmax(Arm5,Arm3)) %>% count(bin1,bin2) 
}



# counting comrades score

count_chimera_1seq <- function(start=1,end=2,Comrades.scores.list=NULL ,indir="/Users/yan/Documents/scientific/cronovirus/G3D/nCov-hiR/analysis/counting_chimeras/", ct.file.name=NULL,outdir="secondary_structue/",name=NULL){
  
  ct.head <- readLines(con = paste0(indir,ct.file.name),n = 1)
  ct <- data.table::fread(paste0(indir,ct.file.name),skip = 1)
  ct$V6 <- c(start:end)
  ct$V2 <- sub(pattern = "[Tt]",replacement = "u",perl = T,x = ct$V2)
  ct <- ct %>% mutate(paired = ifelse(V5==0,0,start+V5 - 1))
  
  writeLines(text =ct.head, con =paste0(outdir,"modified.",ct.file.name), sep = "\n")
  data.table::fwrite(x = ct,file =paste0(outdir,"modified.",ct.file.name),append = T,row.names = F,col.names = F,quote = F,sep = "\t" )
  
  len1 <- end-start+1
  
  
  color.maps <- lapply((names(Comrades.scores.list)),function(x){
    score.sample <- Comrades.scores.list[[x]] %>% as.data.frame() %>% filter(V1 <= end & V1 >= start & V2 <= end & V2 >= start) %>% type.convert() %>%
      mutate(P1 = pmin(V1,V2),P2=pmax(V1,V2)) %>% group_by(P1,P2) %>% summarise(n=sum(V3))
    
    
    color.map <- lapply (1:nrow(ct), function(i){
      
      bin1 <- ct$V6[i]
      paired.i <- ct$paired[i]
      
      if(paired.i == 0){ #not paired, don't count chimeras
        res <- 0
      }
      else{
        # bin2 <- start-1+paired.i
        
        found <- score.sample %>% filter((P1 == bin1 & P2 == paired.i )|(P1 == paired.i & P2 == bin1 ))
        res <- ifelse(nrow(found)>0,sum(found$n) ,0)
      }
      
      data <- data.frame(i=i,value=res) %>% mutate(SPM = log2(1+1e7*value/mapped_RC_group[x]))
    }) %>% bind_rows() 
    
    data.table::fwrite(color.map,file=paste0(outdir,x,"_",name,"_",start,"_",end,".color.map.txt"),sep = "\t",col.names = F,row.names = F)
    max(color.map$SPM) %>% print()
  })
  
}

# count basepairing score for two interacted fragments
count_chimera_2seq <- function(start1=1,end1=2,start2=1,end2=2,Comrades.scores.list=NULL ,indir="/Users/yan/Documents/scientific/cronovirus/G3D/nCov-hiR/analysis/counting_chimeras/", ct.file.name=NULL,outdir="secondary_structue/",name=NULL){
  
  ct.head <- readLines(con = paste0(indir,ct.file.name),n = 1)
  ct <- data.table::fread(paste0(indir,ct.file.name),skip = 1)
  ct$V6 <- c(start1:end1,start2:end2)
  ct$V2 <- sub(pattern = "[Tt]",replacement = "u",perl = T,x = ct$V2)
  
  len1 <- end1-start1+1
  len2 <- end2-start2+1
  # V5 是配对的碱基的相对位置
  
  ct <- ct %>% mutate(paired = ifelse(V5==0,0,ifelse(V5 <= len1,start1-1+V5,V5-len1+start2-1)))
  
  writeLines(text =ct.head, con =paste0(outdir,"modified.",ct.file.name), sep = "\n")
  data.table::fwrite(x = ct,file =paste0(outdir,"modified.",ct.file.name),append = T,row.names = F,col.names = F,quote = F,sep = "\t" )
  
  
  color.maps <- lapply((names(Comrades.scores.list)),function(x){
    score.sample <- Comrades.scores.list[[x]] %>% as.data.frame() %>% filter((V1 <= end1 & V1 >= start1 & V2 <= end2 & V2 >= start2)|(V1 <= end2 & V1 >= start2 & V2 <= end1 & V2 >= start1)) %>%
      mutate(P1 = pmin(V1,V2),P2=pmax(V1,V2)) %>% group_by(P1,P2) %>% summarise(n=sum(V3))
    
    
    color.map <- lapply (1:nrow(ct), function(i){
      
      bin1 <- ct$V6[i]
      paired.i <- ct$paired[i]
      
      if(paired.i == 0){ #not paired, don't count chimeras
        res <- 0
      }
      else{
        # bin2 <- start-1+paired.i
        
        found <- score.sample %>% filter((P1 == bin1 & P2 == paired.i )|(P1 == paired.i & P2 == bin1 ))
        res <- ifelse(nrow(found)>0,sum(found$n) ,0)
      }
      
      data <- data.frame(i=i,value=res) %>% mutate(SPM = log2(1+1e7*value/mapped_RC_group[x]))
    }) %>% bind_rows() 
    
    data.table::fwrite(color.map,file=paste0(outdir,x,"_",name,"_",start1,"_",end1,"_",start2,"_",end2,".color.map.txt"),sep = "\t",col.names = F,row.names = F)
    max(color.map$SPM) %>% print()
  })
  
}

