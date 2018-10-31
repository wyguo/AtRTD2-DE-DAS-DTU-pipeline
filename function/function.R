#############################################################################################
boxplot.normalised <- function(data.before,data.after,condition,sample.name){
  data2plot <- data.before
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g1 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) + 
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x='samples',y='read counts',title='Data distribution before normalisation')
  
  data2plot <- data.after
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g2 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x='samples',y='log2(CPM)',title='Data distribution after normalisation')
  
  if(length(sample.name)>20){
    g1 <- g1+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
    g2 <- g2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
  }
  return(list(g1=g1,g2=g2))
}

#############################################################################################
condition2design <- function(condition,batch.effect=NULL){
  design.data <- data.frame(condition=condition)
  if(!is.null(batch.effect)){
    colnames(batch.effect) <- paste0('batch',1:ncol(batch.effect))
    design.data <- data.frame(design.data,batch.effect)
  }
  design.fomula <- as.formula(paste0('~0+',paste0(colnames(design.data),collapse = '+')))
  design <- model.matrix(design.fomula, data=design.data)
  idx <- grep('condition',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('condition')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design
}

#############################################################################################
counts2CPM <- function (obj, lib.size = NULL,Log=F){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size)) 
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  if(Log){
    t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  } else {
    t(t(counts)/(lib.size) * 1e+06)
  }
  
}

#############################################################################################
###filter cpm
cpm.filter<-function(counts,sample.n=45,cpm.cut=1,...){
  rowSums(counts2CPM(counts,...)>=cpm.cut)>=sample.n
}

TPM.filter<-function(TPM,sample.n=45,tpm.cut=1){
  rowSums(TPM>=tpm.cut)>=sample.n
}


#' Filter low expressed genes
#' @param counts read counts at transcript level.
#' @param mapping transcript-gene name mapping, with first column of transcript list and second column of gene list.
#' @param cpm.cut abundance cut-off. 
#' @param sample.n number of samples cut across all samples.
#' @param Log if TRUE, the log2-CPM is used for filters.
counts.filter <- function(counts,mapping,cpm.cut=1,sample.n=3,Log=F){
  colnames(mapping) <- c('trans','genes')
  rownames(mapping) <- mapping$trans
  filter.idx<-cpm.filter(counts=counts,sample.n = sample.n,cpm.cut =cpm.cut,Log=Log)
  trans_high <- names(filter.idx)[filter.idx==T]
  genes_high <- unique(mapping[trans_high,]$genes)
  list(trans_high=trans_high,genes_high=genes_high,mapping_high=mapping[trans_high,])
}

#############################################################################################
data.error <- function (x, error.type = "stderr"){
  if (error.type == "stderr") 
    error <- sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
  if (error.type == "sd") 
    error <- sd(na.omit(x))
  return(error)
}

#############################################################################################
# pval.cutoff=0.01
# lfc.cutoff=1
# deltaPS.cutoff=0.1
# expressed <- target_high$genes_high
# x <- unique(DE_genes$target)
# y <- unique(DAS_genes$target)
# type='genes'
# 
# load('data/target_high.RData')
# load('data/DE_genes.RData')
# load('data/DAS_genes.RData')
# load('data/DE_trans.RData')
# load('data/DTU_trans.RData')
# plot.flow.chart(expressed = target_high$genes_high,
#                 x = unique(DE_genes$target),
#                 y = unique(DAS_genes$target),
#                 type='genes')
# 
# plot.flow.chart(expressed = target_high$trans_high,
#                 x = unique(DE_trans$target),
#                 y = unique(DTU_trans$target),
#                 type='transcripts')

plot.flow.chart <- function(expressed,x,y,
                            type=c('genes','transcripts'),
                            pval.cutoff=0.01,lfc.cutoff=1,deltaPS.cutoff=0.1){
  type <- match.arg(type,c('genes','transcripts'))
  n.expressed <- length(expressed)
  n.x <- length(x)
  n.notx <- n.expressed-n.x
  n.onlyx <- length(setdiff(x,y))
  n.xy <- length(intersect(x,y))
  n.onlyy <- length(setdiff(y,x))
  n.no<- n.notx-n.onlyy
  n.y <- length(y)
  
  grid.newpage()
  width <- 0.2
  ##---Expressed
  exp.bx <- boxGrob(
    label = paste0('Expressed\n',length(expressed)),
    x = 0.1,
    y = 0.5,
    width = 0.14,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---DE genes
  de.genes.bx <- boxGrob(
    label = paste0('DE ',type,'\n',n.x),
    x = 0.35,
    y = 0.7,
    width = 0.16,
    box_gp=gpar(fill = 'gold'),
    just = 'center')
  
  ##---not DE genes
  notde.genes.bx <- boxGrob(
    label = paste0('Not\n DE ',type,'\n',n.notx),
    x = 0.35,
    y = 0.3,
    width = 0.16,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---Only transcription regulation
  de.only.bx <- boxGrob(
    label = paste0('Only transcription\nregulation\n',n.onlyx),
    x = 0.65,
    y = 0.8,
    width = width,
    box_gp=gpar(fill = 'gold'),
    just = 'center')
  
  ##---Transcription and AS regulation
  de.das.bx <- boxGrob(
    label = paste0('Transcription\nand AS regulation\n',n.xy),
    x = 0.65,
    y = 0.6,
    width = width,
    box_gp=gpar(fill = 'seagreen3'),
    just = 'center')
  
  ##---Only AS regulation
  das.only.bx <- boxGrob(
    label = paste0('Only AS\nregulation\n',n.onlyy),
    x = 0.65,
    y = 0.4,
    width = width,
    box_gp=gpar(fill = 'deepskyblue3'),
    just = 'center')
  
  ##---No regulation
  no.regulation.bx <- boxGrob(
    label = paste0('No\nregulation\n',n.no),
    x = 0.65,
    y = 0.2,
    width = width,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---DAS genes
  das.bx <- boxGrob(
    label = ifelse(type=='genes',paste0('DAS genes\n',n.y),paste0('DTU transcripts\n',n.y)),
    x = 0.9,
    y = 0.5,
    width = 0.15,
    box_gp=gpar(fill = 'deepskyblue3'),
    just = 'center')
  
  
  arrow_obj <- getOption('connectGrobArrow', 
                         default = arrow(angle = 0, ends = 'last', type = 'closed'))
  lty_gp <- getOption('connectGrob', default = gpar(fill = 'black',lwd=4))
  plot(connectGrob(exp.bx, de.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(exp.bx, notde.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(de.genes.bx, de.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(de.genes.bx, de.das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(notde.genes.bx, das.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(notde.genes.bx, no.regulation.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  
  lty_gp <- getOption('connectGrob', default = gpar(col = 'red',lwd=4))
  plot(connectGrob(de.das.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(das.only.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  
  plot(exp.bx)
  plot(de.genes.bx)
  plot(notde.genes.bx)
  plot(de.only.bx)
  plot(de.das.bx)
  plot(das.only.bx)
  plot(no.regulation.bx)
  plot(das.bx)
  grid.text(ifelse(type=='genes','DE and DAS genes','DE and DTU transcripts'), x=0.5, y=0.95, 
            gp=gpar(fontsize=16, col='black'))
  grid.text(paste0('P-value cut-off: ',pval.cutoff,
                   '\nLog2 fold change cut-off: ',lfc.cutoff,
                   '\ndeltaPS cut-off: ', deltaPS.cutoff),
            x=0.05, y=0.1, just = 'left',
            gp=gpar(fontsize=12, col='black'))
}



#############################################################################################
##-----------------------------------------------------------------------------------------
edgeR.pipeline <- function(dge,
                           method=c('glm','glmQL'),
                           design,
                           contrast,
                           diffAS=F,
                           deltaPS=NULL,
                           adjust.method='BH'){
  
  start.time <- Sys.time()
  results <- list()
  method <- match.arg(method,c('glm','glmQL'))
  targets <- rownames(dge$counts)
  message('Estimate dispersion ...')
  Disp <- estimateGLMCommonDisp(dge, design)
  Disp <- estimateGLMTagwiseDisp(Disp, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  
  switch(method, 
         glm = {
           message('Fit Genewise Negative Binomial Generalized Linear Models...')
           fit <- glmFit(Disp, design = design)
           results$fit.glm <- fit
           
           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmLRT(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           ###overall test
           x <- glmLRT(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         },
         glmQL = {
           message('Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests...')
           fit <- glmQLFit(Disp,design = design)
           results$fit.glmQL <- fit
           
           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmQLFTest(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           
           ###overall test
           x <- glmQLFTest(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         }
  )
  names(DE.pval.list) <- contrast
  results$DE.pval.list <- DE.pval.list
  
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$FDR))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast
  
  DE.stat <- summary.stat(x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  ##################################################################
  ###---DAS
  if(diffAS){
    fit.splice <- lapply(contrast,function(i){
      message(paste0('\ndiffSpliceDGE of contrast: ', i))
      diffSpliceDGE(fit, contrast=contrast.matrix[,i,drop=F], geneid="GENEID")
    })
    names(fit.splice) <- contrast
    results$fit.splice<-fit.splice
    
    genes.idx <- unique(fit.splice[[1]]$genes$GENEID)
    trans.idx <- unique(fit.splice[[1]]$genes$TXNAME)
    
    DTU.pval.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="exon", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast
    
    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast
    
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),]
    DTU.stat <- summary.stat(x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))
    
    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat
    
    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice[[1]]$genes$TXNAME,],
                     INDICES = fit.splice[[1]]$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     })
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,]
    results$maxdeltaPS<-maxdeltaPS
    
    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="gene", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.F.list) <- contrast
    
    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,]
    colnames(DAS.pval.F) <- contrast
    
    DAS.F.stat <- summary.stat(x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat
    
    ##---Simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="Simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.simes.list) <- contrast
    
    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,]
    colnames(DAS.pval.simes) <- contrast
    
    DAS.simes.stat <- summary.stat(x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))
    
    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}



#############################################################################################
gg.color.hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

distinct.color <- function(n){
  col.lib <- c("#C0D57E","#BF6CF4","#D435BB","#C2BBCB","#CAE678","#BEC3AD","#EC75F6","#F8DAD6","#D9FBD7","#46C3CF","#5FECDE","#AEE1F5","#AA67CC",
               "#81C1F1","#D6A381","#81DADB","#99AD6D","#B8ECF3","#828D6F","#F832A9","#8DD092","#7FC0D0","#6EDDD0","#D4BE95","#B8F7F2","#9752C4",
               "#AF9D7C","#468CE2","#9AA843","#4F9F7A","#8C7DB6","#CA44D8","#93D0B9","#4A718E","#F259C8","#E1C358","#BB924E","#5682F1","#B3E18E",
               "#59C1A6","#3BA7EB","#4ED06A","#A28FF0","#F4F2E7","#F7C2DA","#6921ED","#E8AAE6","#C2F2DF","#A3F7A1","#F0A85C","#BFF399","#E93B87",
               "#D7C7F0","#BCF7B6","#39B778","#DD714B","#B7D9A6","#E19892","#BF70C4","#B882F4","#74F38B","#D6E3C2","#9E23B2","#65F55F","#97F3EB",
               "#6BCD37","#8A7B8E","#DEEDA7","#F03E4F","#8FE77B","#ACF165","#8FCAC4","#CA1FF5","#E854EE","#43E11C","#6EEDA7","#B584B3","#BAAAF2",
               "#DE88C9","#415B96","#4EF4F3","#81F8CA","#C1B4AC","#BBA4CF","#AB55E5","#F789EC","#CEDFF8","#A9BEF0","#658DC6","#94ACEA","#91CD6F",
               "#BAA431","#8FB2C5","#F8F1A3","#D1CCCE","#F533D8","#844C56","#F3DA7D","#BDE5C0","#F1F69C","#F17C85","#C5E743","#B8C8ED","#3ED3EC",
               "#E6C5D9","#AD7967","#EAF4C3","#3F4DC1","#F4DDEC","#F1B136","#C247F3","#6FD3A1","#DCF7A3","#B8D1C2","#545CED","#2C7A29","#E4B1BB",
               "#785095","#8373F1","#D0F2C5","#B53757","#49F79F","#F0BE8B","#9630DD","#D6D89F","#79EED3","#EEEFD7","#7FA41F","#D1BF71","#6AEFF6",
               "#D661D8","#CFEEF8","#9ECB59","#67D382","#A8945E","#BB9B9E","#F3E4B9","#9964F3","#835BD4","#A5E2BF","#9CA6C8","#CD9BE9","#8FF038",
               "#B6C08A","#D19BB5","#8893EB","#51746E","#5920AC","#F3A1BB","#D7E1DA","#D668C3","#EC9EF2","#C2738A","#EDE979","#615AB4","#E9EF68",
               "#F3DECA","#90A6A0","#F2B8E7","#36DEA3","#EED89A","#5AADA5","#92EFAD","#A3288A","#ECC9F4","#91C9E7","#ECC2B8","#EB9BCC","#2735BB",
               "#5E7A2D","#C3BDF7","#DCF8F4","#52E0BC","#C2F8D6","#B3CFDA","#6E93A2","#9AECCA","#E3D8F0","#D97EE7","#9CBE9A","#60A858","#55F9D7",
               "#F47AA2","#DCF5E3","#B55398","#AC81D3","#644AF2","#546DEA","#E7EB37","#EA7E2A","#A1517E","#E131E4","#DE92F4","#EC66B1","#82DDF3",
               "#E99475","#E4EAF3","#E8D144","#40C0EE","#559ABA")
  col.lib[1:n]
}

gg_legend <- function(g){ 
  tmp <- ggplot_gtable(ggplot_build(g)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 



#############################################################################################
limma.pipeline <- function(dge,
                           contrast,
                           span=0.5,
                           design,
                           deltaPS=NULL,
                           diffAS=F,
                           adjust.method='BH'){
  start.time <- Sys.time()
  results <- list()
  ##########################################################
  ##--limma voom
  message('Limma-voon to estimate mean-vriance trend ...')
  voom.oboject<-voom(dge,design,plot=F,span = span)
  results$voom.oboject <- voom.oboject
  targets <- rownames(voom.oboject$E)
  
  ##########################################################
  ##--Fit a basic linear model
  message('Fit a basic linear model ...')
  fit.lmFit <- lmFit(voom.oboject, design)
  results$fit.lmFit <- fit.lmFit
  
  ##########################################################
  ##--Fit a basic linear model
  message('Fit the contrast model ...')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
  fit.contrast<-contrasts.fit(fit.lmFit, contrast.matrix)
  results$fit.contrast<-fit.contrast
  
  ##########################################################
  ##--Fit a eBayes model
  message('Fit a eBayes model ...')
  fit.eBayes<-eBayes(fit.contrast)
  results$fit.eBayes<-fit.eBayes
  
  ##########################################################
  ##--Testing statistics for each contrast group
  message('Testing for each contrast group ...')
  DE.pval.list <- lapply(contrast,function(i){
    x <- topTable(fit.eBayes, 
                  coef=i, 
                  adjust.method =adjust.method,
                  number = Inf)
    x <- x[targets,]
    x
  })
  names(DE.pval.list) <- contrast
  results$DE.pval.list<-DE.pval.list
  
  ###---DE pval and lfc
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast
  
  DE.stat <- summary.stat(x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  
  # DE.stat <- cbind(DE.pval,DE.lfc)
  # colnames(DE.stat) <- c(paste0('pval:',contrast),paste0('lfc:',contrast))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  
  ##########################################################
  ##--Testing statistics for across all contrast groups
  message('Testing across all contrast groups ...')
  DE.stat.overalltest<-topTable(fit.eBayes,number = Inf,
                                coef = contrast,adjust.method =adjust.method )
  # DE.stat.overalltest <- DE.stat.overalltest[targets,]
  col.idx <- gsub('-','.',contrast)
  col.idx <- grep(paste0(col.idx,collapse = '|'),colnames(DE.stat.overalltest))
  DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                    contrast='overall',DE.stat.overalltest,row.names = NULL)
  colnames(DE.stat.overalltest)[col.idx] <- gsub('[.]','-',colnames(DE.stat.overalltest)[col.idx])
  results$DE.stat.overalltest<-DE.stat.overalltest
  
  if(diffAS){
    message('Fit a splicing model ...')
    if(is.null(deltaPS))
      stop('Please provide deltaPS for DAS analysis...')
    fit.splice<-diffSplice(fit.contrast, geneid = 'GENEID')
    results$fit.splice<-fit.splice
    
    # ##########################################################
    ##---DTU transcripts
    ##DTU transcript pval list
    genes.idx <- unique(fit.splice$genes$GENEID)
    trans.idx <- unique(fit.splice$genes$TXNAME)
    
    DTU.pval.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="t", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast
    
    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast
    
    ##DTU transcript deltaPS
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),]
    
    DTU.stat <- summary.stat(x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))
    
    # DTU.stat <- cbind(DTU.pval,DTU.deltaPS)
    # colnames(DTU.stat) <- c(paste0('pval:',contrast),paste0('deltaPS:',contrast))
    
    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat
    
    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice$genes$TXNAME,],
                     INDICES = fit.splice$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     })
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,]
    results$maxdeltaPS<-maxdeltaPS
    
    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="F", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.F.list) <- contrast
    
    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,]
    colnames(DAS.pval.F) <- contrast
    
    DAS.F.stat <- summary.stat(x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    
    # DAS.F.stat <- cbind(DAS.pval.F,maxdeltaPS)
    # colnames(DAS.F.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    # 
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat
    
    ##---simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.simes.list) <- contrast
    
    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,]
    colnames(DAS.pval.simes) <- contrast
    
    DAS.simes.stat <- summary.stat(x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))
    # 
    # DAS.simes.stat <- cbind(DAS.pval.simes,maxdeltaPS)
    # colnames(DAS.simes.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    
    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}


#############################################################################################
mean.variance.trend <- function (obj, design = NULL,lib.size = NULL, 
                                 normalize.method = "none",                              
                                 span = 0.5, ...){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size)) 
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  list(fit=fit,sx=sx,sy=sy,l=l)
} 


plot.mean.variance <- function(x,y,l,fit.line.col='red',
                               xlab = "log2( count size + 0.5 )", 
                               ylab = "Sqrt(standard deviation)",
                               main="Mean-variance trend",lwd=1.5,...){
  plot(x, y, pch = 16, cex = 0.25,xlab=xlab,ylab=ylab,main=main,...)
  lines(l,col=fit.line.col,lwd=lwd)
}


check.mean.variance <- function(counts.raw,
                                counts.filtered,
                                condition,
                                span = 0.5,
                                points.col1='black',
                                points.col2='black',make.plot=T,...){
  ###---design matrix
  condition <- factor(condition,levels = unique(condition))
  design<-model.matrix(~0+condition)
  colnames(design) <- gsub('condition','',colnames(design))
  
  ###---generate DGEList object
  message('=> Generate DGEList object')
  dge.raw <- DGEList(counts=counts.raw)
  dge.raw <- calcNormFactors(dge.raw)
  dge.filtered<- DGEList(counts=counts.filtered)
  dge.filtered <- calcNormFactors(dge.filtered)
  
  ###---fit mean-variance trend
  message('=> Fit mean-variance trend')
  #before filter
  fit.raw <- mean.variance.trend(obj = dge.raw,design = design,span = span)
  #after filter
  fit.filtered <- mean.variance.trend(obj = dge.filtered,design = design,span = span)
  
  ###---plot mean-variance trend
  if(make.plot){
    message('=> Plot mean-variance trend')
    par(mfrow=c(1,2))
    plot.mean.variance(x = fit.raw$sx,y = fit.raw$sy,
                       l = fit.raw$l,lwd=2,fit.line.col ='gold',col=points.col1,...)
    title('\n\nRaw counts')
    plot.mean.variance(x = fit.filtered$sx,y = fit.filtered$sy,
                       l = fit.filtered$l,lwd=2,col=points.col2,...)
    title('\n\nFiltered counts')
    lines(fit.raw$l, col = "gold",lty=4,lwd=2)
    legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
           legend = c('low-exp removed','low-exp kept'))
    message('Done!!!')
  } else {
    message('Done!!!')
    return(list(fit.raw=fit.raw,fit.filtered=fit.filtered))
  }
  
}

#############################################################################################
plot.DDD.bar <- function(targets,contrast,title=NULL){
  data2plot <- lapply(contrast,function(i){
    if(nrow(targets[[i]])==0){
      x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
    } else {
      x <- data.frame(contrast=i,table(targets[[i]]$up.down))
      colnames(x) <- c('contrast','regulation','number')
    }
    x
  })
  data2plot <- do.call(rbind,data2plot)
  data2plot$contrast <- factor(data2plot$contrast,levels = unique(data2plot$contrast))
  data2plot$regulation <- factor(data2plot$regulation,levels = c('up_regulate','down_regulate'))
  g <- ggplot(data2plot,aes(x=contrast,y=number,group=regulation,
                            fill=regulation,label=number))+
    geom_bar(stat='identity')+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme_bw()+
    labs(title=title)
  if(length(data2plot$contrast)>10)
    g <- g+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g 
}

plot.euler.diagram <- function(x,
                               fill=gg.color.hue(length(x)),
                               shape = "ellipse"){
  fit <- euler(x)
  # grid.newpage()
  g <- plot(fit,quantities = TRUE,
            labels = list(font =1),
            fill=fill,shape = shape)
  g
}

# targets <- lapply(contrast,function(i){
#   DE_genes[[i]]$target
# })
# names(targets) <- contrast
# g <- plot.euler.diagram(x = targets)
# grid.arrange(g,top=textGrob('DE genes', gp=gpar(cex=1.2)))


# 
# load('data/trans_TPM.RData')
# gene=DE_genes$target[1]
# DE.genes <- unique(DE_genes$target)
# DAS.genes <- unique(DAS_genes$target)
# DE.trans <- unique(DE_trans$target)
# DTU.trans <- unique(DTU_trans$target)
# 
# genes.ann <- venn2(DE.genes,DAS.genes)
# names(genes.ann) <- c('DEonly','DE&DAS','DASonly')
# genes.ann <- plyr::ldply(genes.ann,cbind)[,c(2,1)]
# genes.ann[1:10,]
# colnames(genes.ann) <- c('target','annotation')
# head(genes.ann)
# 
# trans.ann <- venn2(DE.trans,DTU.trans)
# names(trans.ann) <- c('DEonly','DE&DTU','DTUonly')
# trans.ann <- plyr::ldply(trans.ann,cbind)[,c(2,1)]
# trans.ann[1:10,]
# colnames(trans.ann) <- c('target','annotation')
# head(trans.ann)
# 


plot.abundance<- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  trans.expressed=NULL,
  reps,
  x.lab='Conditions',
  y.lab='TPM',
  marker.size=3,
  legend.text.size=11,
  error.type='stderr',
  show.annotation=T,
  plot.gene=T,
  error.bar=T
){
  #####################################################################
  ##prepare plot data                                                
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  plot.title <- paste0('Gene: ',gene)
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]
  
  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }
  
  
  data2plot<-data.exp[trans,,drop=F]
  
  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)
  
  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(!is.null(genes.ann)){
      idx <- genes.ann[gene,'annotation']
      if(length(idx)>0)
        genes.idx <- paste0(gene,' (',idx,')') 
    }
    
    if(!is.null(trans.ann)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      if(length(idx)>0)
        trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)
  legend.ncol <- ceiling(length(trans.idx)/15)
  
  ##--mean and error
  mean2plot <- t(rowmean(t(data2plot),group = reps))
  sd2plot <- by(t(data2plot),INDICES = reps,FUN = function(x){
    apply(x,2,function(y){
      data.error(y,error.type = error.type)
    })
  })
  sd2plot <- do.call(cbind,sd2plot)
  sd2plot <- sd2plot[rownames(mean2plot),]
  
  mean2plot <- reshape2::melt(mean2plot)
  sd2plot <- reshape2::melt(sd2plot)
  colnames(mean2plot) <- c('Targets','Conditions','mean')
  colnames(sd2plot) <- c('Targets','Conditions','error')
  data2plot <- merge(mean2plot,sd2plot)
  
  ##--plot gene or not
  if(!plot.gene){
    data2plot <- data2plot[-which(data2plot$Targets %in% genes.idx),]
    plot.title <- paste0('Gene: ',genes.idx)
  }
  
  
  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))
  
  profiles <- ggplot(data2plot,aes(x=Conditions,y=mean,group=Targets,color=Targets))+
    geom_line(size=1)+
    geom_point(size=marker.size,aes(fill=Targets,shape=Targets))+
    scale_shape_manual(name="Targets",values=c(25:0,25:0,rep(25,500)))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  if(error.bar)
    profiles <- profiles+
    geom_errorbar(aes(ymin=mean-error,ymax=mean+error),width=0.1,color='black',size=0.3)
  profiles
}



plot.PS <- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  show.annotation=T,
  trans.expressed=NULL,
  reps,
  y.lab='TPM',
  x.lab='Conditions',
  marker.size=3,
  legend.text.size=11
){
  #####################################################################
  ##prepare plot data                                                
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]
  
  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }
  data2plot<-data.exp[trans,,drop=F]
  data2plot <- t(rowmean(t(data2plot),reps))
  
  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)
  
  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(!is.null(genes.ann)){
      idx <- genes.ann[gene,'annotation']
      if(length(idx)>0)
        genes.idx <- paste0(gene,' (',idx,')') 
    }
    
    if(!is.null(trans.ann)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      if(length(idx)>0)
        trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)
  legend.ncol <- ceiling(length(trans.idx)/15)
  plot.title <- paste0('Gene: ',genes.idx)
  ##PS value
  data2plot[trans.idx,] <- t(t(data2plot[trans.idx,])/data2plot[genes.idx,])
  data2plot <- data2plot[-which(rownames(data2plot) %in% genes.idx),,drop=F]
  data2plot <- reshape2::melt(data2plot)
  colnames(data2plot) <- c('Targets','Conditions','PS')
  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))
  
  profiles <- ggplot(data2plot,aes(x=Conditions,y=PS,group=Targets,color=Targets))+
    geom_bar(stat='identity',aes(fill=Targets))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  profiles
}


#############################################################################################

#' counts.filtered = trans_counts[targets_high$trans_high,]
#' dge <- DGEList(counts=counts.filtered) %>% calcNormFactors()
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#' rownames(data2pca) <- paste0(samples$baseline,'.',samples$brep)
#' groups<-samples$brep
#' dim1 <- 'PC1'
#' dim2 <- 'PC2'
#' plot_PCA_ind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = rownames(data2pca))


plot.PCA.ind <- function(data2pca,dim1='PC1',dim2='PC2',
                         groups,
                         title='PCA plot',
                         ellipse.type=c('none','ellipse','polygon'),
                         add.label=T,adj.label=F){
  
  ellipse.type <- match.arg(ellipse.type,c('none','ellipse','polygon'))
  dim1 <- toupper(dim1)
  dim2 <- toupper(dim2)
  
  fit <- prcomp(data2pca,scale = T)
  fit.stat <- summary(fit)$importance
  dim1.p <- round(fit.stat[2,dim1]*100,2)
  dim2.p <- round(fit.stat[2,dim2]*100,2)
  data2plot <- data.frame(groups=groups,dim1=fit$x[,dim1],dim2=fit$x[,dim2])
  data2plot$groups <- factor(data2plot$groups,levels = unique(data2plot$groups))
  data2plot$labels <- rownames(data2plot)
  
  g <- ggplot(data2plot,aes(x=dim1,y=dim2),frame=T)+
    geom_point(aes(colour=groups,shape=groups))+
    geom_hline(yintercept = 0,linetype=2)+
    geom_vline(xintercept=0,linetype=2)+
    theme_bw()+
    theme(panel.border = element_blank())+
    labs(x=paste0(dim1,' (',dim1.p,'%)'),
         y=paste0(dim2,' (',dim2.p,'%)'),
         title=title)+
    scale_shape_manual(values = rep(c(18:0,19:25),3)[1:length(unique(groups))])+
    coord_cartesian(xlim = c(min(data2plot$dim1)*1.1,max(data2plot$dim1)*1.1),
                    ylim = c(min(data2plot$dim2)*1.1,max(data2plot$dim2)*1.1))
  
  if(ellipse.type=='ellipse')
    g <- g+stat_ellipse(geom = "polygon",aes(fill=groups,colour=groups),alpha=0.2)
  if(ellipse.type=='polygon'){
    hulls <- plyr::ddply(data2plot, "groups", function(x){
      x[chull(x$dim1,x$dim2),]
    })
    g <- g+geom_polygon(data = hulls,aes(x=dim1,y=dim2,fill=groups),alpha = 0.3)
  }
  if(add.label & !adj.label)
    g <- g+geom_text(aes(label=labels,colour=groups),hjust=0, vjust=0)
  if(add.label & adj.label)
    g <- g+ggrepel::geom_text_repel(label=rownames(data2pca),aes(colour=groups))
  g
  
}


#############################################################################################

# go.table <- readr::read_csv('data/DAS genes GO terms.csv')
# col.lab <- '-log10(FDR)'
plotGO <- function(go.table,col.lab,plot.title='GO annotation'){
  data2plot <- go.table[,c(1,2,which(colnames(go.table)==col.lab))]
  colnames(data2plot) <- c('Category','Term','Value')
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  data2plot$Term <- factor(data2plot$Term,levels = rev(data2plot$Term))
  g <- ggplot(data2plot,aes(x=Term,y=Value))+
    geom_bar(stat='identity',aes(fill=Category))+
    coord_flip() +
    theme_bw()+
    theme(axis.title.y = element_blank())+
    facet_grid(Category~.,scales = 'free_y',space='free_y')+
    labs(y=col.lab,title=plot.title)
  g
}


#############################################################################################


remove_batch<-function(read.counts,design,group){
  start.time <- Sys.time()
  read.counts<-round(as.matrix(read.counts),0)
  message('Estimate norm factor...')
  y <- DGEList(counts=read.counts, group=group)
  y <- calcNormFactors(y)
  message('Estimate Common Dispersion for Negative Binomial GLMs ...')
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  # seqUQ <- betweenLaneNormalization(read.counts, which="upper")
  controls = rep(TRUE,dim(read.counts)[1])
  message('Remove Unwanted Variation Using Residuals...')
  # batch_ruv_res = RUVr(seqUQ,controls,k=1,res) 
  batch_ruv_res = RUVr(read.counts,controls,k=1,res)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
  message('Done!!!')
  return(c(batch_ruv_res,res=list(res)))
}


# load('data/samples.new.RData')
# load('data/targets_high.RData')
# load('data/trans_counts.RData')
# 
# read.counts <- trans_counts[targets_high$trans_high,]
# condition<-paste0(samples.new$condition)
# condition<-factor(condition,levels = unique(condition))
# design<-model.matrix(~0+condition)
# colnames(design)<-gsub('condition','',colnames(design))
# group <- samples.new$condition
# contrast <- c('T10-T2','T19-T2')
# contrast <- makeContrasts(contrast,levels = design)
# 
# 
# cIdx = rep(TRUE,dim(read.counts)[1])
# 
# xar <- remove.batch(read.counts = read.counts,method = 'RUVr',condition = condition,group = condition,cIdx = cIdx)
# xag <- remove.batch(read.counts = read.counts,method = 'RUVg',condition = condition,group = condition,cIdx = cIdx)
# xas <- remove.batch(read.counts = read.counts,method = 'RUVs',condition = condition,group = condition,cIdx = cIdx)
# 
# x1 <- t(counts2CPM(read.counts,Log = T))
# xx1 <- cbind(brep=samples.new$brep,x1)
# g1 <- autoplot(prcomp(x1,scale. = T),label=T,data = xx1,colour='brep')
# 
# x2 <- t(counts2CPM(xar$normalizedCounts,Log = T))
# xx2 <- cbind(brep=samples.new$brep,x2)
# g2 <- autoplot(prcomp(x2,scale. = T),label=T,data = xx2,colour='brep')
# 
# x3 <- t(counts2CPM(xag$normalizedCounts,Log = T))
# xx3 <- cbind(brep=samples.new$brep,x3)
# g3 <- autoplot(prcomp(x3,scale. = T),label=T,data = xx3,colour='brep')
# 
# x4 <- t(counts2CPM(xas$normalizedCounts,Log = T))
# xx4 <- cbind(brep=samples.new$brep,x4)
# g4 <- autoplot(prcomp(x4,scale. = T),label=T,data = xx4,colour='brep')
# 
# gridExtra::grid.arrange(g1,g2,g3,g4,ncol=2)

# results <- remove.batch(read.counts = read.counts,
#                         condition = condition,
#                         design = design,
#                         contrast=contrast,
#                         group = condition,
#                         method = 'RUVs')

remove.batch<-function(read.counts,
                       condition,
                       design=NULL,
                       contrast=NULL,
                       group=NULL,
                       method=c('RUVr','RUVg','RUVs'),
                       cIdx=NULL,
                       k=1,...){
  start.time <- Sys.time()
  method <- match.arg(method,c('RUVr','RUVg','RUVs'))
  if(is.null(design)){
    condition<-factor(condition,levels = unique(condition))
    design<-model.matrix(~0+condition)
    colnames(design)<-gsub('condition','',colnames(design))
  }
  if(is.null(contrast)){
    contrast <- unique(condition)
    contrast <- paste0(contrast[-1],'-',contrast[1])
  }
  
  contrast <- makeContrasts(contrasts = contrast, levels=design)
  read.counts<-as.matrix(round(read.counts,0))
  message('Estimate norm factor...')
  y <- DGEList(counts=read.counts, group=group)
  y <- calcNormFactors(y)
  message('Estimate Common Dispersion for Negative Binomial GLMs ...')
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
  
  fit <- glmFit(y, design)
  ##########################################################################
  ##---get the negative control
  if(is.null(cIdx)){
    lrt <- glmLRT(fit, contrast = contrast)
    top <- topTags(lrt, n=nrow(read.counts))$table
    empirical <- top[order(top$PValue,decreasing = T),]
    empirical <- rownames(empirical)[empirical$PValue>0.1]
  }
  
  switch(method, 
         RUVr = {
           message('Remove Unwanted Variation Using RUVr...')
           res <- residuals(fit, type="deviance")
           results <-  RUVr(x = read.counts,cIdx = empirical,residuals = res,k=k)
         },
         RUVg = {
           message('Remove Unwanted Variation Using RUVg...')
           results <- RUVg(x = read.counts,cIdx = empirical,k=k)
         },
         RUVs = {
           message('Remove Unwanted Variation Using RUVs...')
           results <- RUVs(x = read.counts,cIdx = empirical,k=k,scIdx = makeGroups(condition))
         }
  )
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
  message('Done!!!')
  results$method <- method
  return(results)
}


remove.batch.shiny<-function(read.counts,
                             condition,
                             design=NULL,
                             contrast=NULL,
                             group=NULL,
                             method=c('RUVr','RUVg','RUVs'),
                             cIdx=NULL,
                             k=1,...){
  withProgress(message = 'Estimate batch effects...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0.1)
                 start.time <- Sys.time()
                 method <- match.arg(method,c('RUVr','RUVg','RUVs'))
                 if(is.null(design)){
                   condition<-factor(condition,levels = unique(condition))
                   design<-model.matrix(~0+condition)
                   colnames(design)<-gsub('condition','',colnames(design))
                 }
                 if(is.null(contrast)){
                   contrast <- unique(condition)
                   contrast <- paste0(contrast[-1],'-',contrast[1])
                 }
                 
                 contrast <- makeContrasts(contrasts = contrast, levels=design)
                 read.counts<-as.matrix(round(read.counts,0))
                 
                 message('Estimate norm factor...')
                 incProgress(0.2)
                 y <- DGEList(counts=read.counts, group=group)
                 y <- calcNormFactors(y)
                 
                 message('Estimate Common Dispersion for Negative Binomial GLMs ...')
                 incProgress(0.3)
                 y <- estimateGLMCommonDisp(y, design)
                 y <- estimateGLMTagwiseDisp(y, design)
                 
                 message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
                 incProgress(0.4)
                 fit <- glmFit(y, design)
                 ##########################################################################
                 ##---get the negative control
                 # if(is.null(cIdx)){
                 #   lrt <- glmLRT(fit, contrast = contrast)
                 #   top <- topTags(lrt, n=nrow(read.counts))$table
                 #   empirical <- top[order(top$PValue,decreasing = T),]
                 #   empirical <- rownames(empirical)[empirical$PValue>0.1]
                 # }
                 empirical <- rep(TRUE,dim(read.counts)[1])
                 
                 switch(method, 
                        RUVr = {
                          message('Remove Unwanted Variation Using RUVr...')
                          res <- residuals(fit, type="deviance")
                          results <-  RUVr(x = read.counts,cIdx = empirical,residuals = res,k=k)
                        },
                        RUVg = {
                          message('Remove Unwanted Variation Using RUVg...')
                          results <- RUVg(x = read.counts,cIdx = empirical,k=k)
                        },
                        RUVs = {
                          message('Remove Unwanted Variation Using RUVs...')
                          results <- RUVs(x = read.counts,cIdx = empirical,k=k,scIdx = makeGroups(condition))
                        }
                 )
                 incProgress(0.7)
                 end.time <- Sys.time()
                 time.taken <- end.time - start.time
                 message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
                 showNotification("Done!!!")
                 message('Done!!!')
                 results$method <- method
                 return(results)
               })
}



#############################################################################################
reorder.clusters <- function(clusters,dat){
  targets <- rownames(dat)
  names(clusters) <- targets
  idx0 <- rowmean(dat,group = clusters,reorder = T)
  para <- lapply(1:ncol(idx0),function(i) idx0[,i])
  idx0 <- idx0[do.call(order,para),]
  idx0 <- idx0[rev(rownames(idx0)),]
  order.idx <- rownames(idx0)
  idx <- 1:nrow(idx0)
  names(idx) <- order.idx
  idx <- idx[as.character(clusters)]
  names(idx) <- names(clusters)
  idx
}

#############################################################################################
rowmean <- function (x, group, reorder = T, na.rm = T) {
  order.idx <- as.character(unique(group))
  if (reorder) 
    order.idx <- gtools::mixedsort(order.idx)
  counts <- table(group)[order.idx]
  sums <- rowsum(x, group = group)[order.idx, ]
  means <- (diag(1/counts)) %*% as.matrix(sums)
  rownames(means) <- order.idx
  if (na.rm) 
    means[is.na(means)] <- 0
  return(means)
}

#############################################################################################
rowratio <- function (x, group, reorder = T, na.rm = T){
  y <- rowsum(x, group = group)
  y <- y[group, ]
  ratio = x/y
  rownames(ratio) <- rownames(x)
  if (na.rm) 
    ratio[is.na(ratio)] <- 0
  if (reorder & !is.null(rownames(ratio))) 
    ratio <- ratio[gtools::mixedsort(rownames(ratio)), ]
  return(ratio)
}


#############################################################################################
#' Average Over Replicate Arrays
#' 
#' @param x a matrix-like object.
#' @param group grouping of sample identifier.
#' @return A data object of the same class as x with a column for sums according grouping.
#' @examples 
#' set.seed(100)
#' x<- matrix(rnorm(8*4),8,4)
#' colnames(x) <- c("a","a","b","b")
#' sumarrays(x)
#' data.frame(a=x[,1]+x[,2],b=x[,3]+x[,4])

sumarrays<-function(x,group=NULL){
  if(is.null(group))
    group<-colnames(x)
  colnames(x)<-group
  
  x<-rowsum(t(x),group = group)
  #order the columns as the input group odering
  x<-data.frame(t(x)[,unique(group)])
  return(x)
}


#############################################################################################
summary.stat <- function(x,y,target,contrast,stat.type=c('adj.pval','lfc'),srot.by=stat.type[1]){
  x <- x[target,]
  y <- y[target,]
  stat <- lapply(contrast,function(i){
    z <- data.frame(targets =target,contrast=i, x[,i],y[,i],row.names = NULL)
    colnames(z) <- c('target','contrast',stat.type)
    z <- z[order(z[,srot.by]),]
  })
  names(stat) <- contrast
  stat <- do.call(rbind,stat)
  rownames(stat) <- NULL
  stat
}



##-----------------------------------------------------------------------------------------
# lfc <- trans_DDD_stat$DE.lfc
# lfc <- reshape2::melt(as.matrix(lfc))
# colnames(lfc) <- c('target','contrast','log2FC')
# 
# dim(lfc)
# dim(stat)
# 
# x <- trans_DDD_stat$DTU.pval
# y <- trans_DDD_stat$DTU.deltaPS
# 
# stat <- genes_DDD_stat$DE.stat

summary.DE.target <- function(stat,cutoff=c(adj.pval=0.01,log2FC=1)){
  names(cutoff) <- c('adj.pval','log2FC')
  stat$up.down <- 'up_regulate'
  stat$up.down[stat[,'log2FC']<0] <- 'down_regulate'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,'log2FC'])>=cutoff['log2FC'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}


##-----------------------------------------------------------------------------------------
# stat <- trans_DDD_stat$DTU.stat
# lfc <- trans_DDD_stat$DE.lfc
# 
# stat <- trans_DDD_stat$DAS.F.stat
# lfc <- genes_DDD_stat$DE.lfc
# 
# 
# lfc <- reshape2::melt(as.matrix(lfc))
# colnames(lfc) <- c('target','contrast','log2FC')

summary.DAS.target <- function(stat,lfc,cutoff=c(adj.pval=0.01,deltaPS=0.1)){
  names(cutoff) <- c('adj.pval','deltaPS')
  lfc <- lfc[which(lfc$target %in% stat$target),]
  stat <- merge(stat,lfc)
  stat$up.down <- 'up_regulate'
  stat$up.down[stat[,'log2FC']<0] <- 'down_regulate'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,grep('deltaPS',colnames(stat))])>=cutoff['deltaPS'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}



summary.target <- function(stat,lfc=NULL,
                           stat.col=c('adj.pval','log2FC'),
                           cutoff=c(0.01,1)){
  if(!is.null(lfc)){
    
  }
  
  stat$up.down <- 'up_regulate'
  stat$up.down[stat[,stat.col[2]]<0] <- 'down_regulate'
  idx <- (abs(stat[,stat.col[1]])<=cutoff[1]) & (abs(stat[,stat.col[2]])>=cutoff[2])
  stat <- stat[idx,]
  split(stat , f = stat$contrast )
}

##-----------------------------------------------------------------------------------------
# x <- DE_genes
# y <- DAS_genes
# summary.vs <- function(x,y,contrast,type='DE vs DAS'){
#   idx <- unlist(strsplit(type,' vs '))
#   idx <- c(idx[1],paste0(idx[1],'&',idx[2]),idx[2])
#   vs <- lapply(contrast,function(i){
#     x <- venn2(x[[i]]$target,y[[i]]$target)
#     names(x) <- idx
#     x
#   })
#   names(vs) <- contrast
#   n1 <- do.call(rbind,lapply(vs,function(i) sapply(i,length)))
#   vs.expand <- unlist(vs,recursive=FALSE)
#   n2 <- sapply(idx,function(i){
#     length(Reduce(intersect,vs.expand[grep(paste0('\\.',i,'$'),names(vs.expand))]))
#   })
#   n2 <- data.frame(t(n2),row.names = 'Intersection')
#   colnames(n2) <- colnames(n1)
#   x <- rbind(n1,n2)
#   x <- data.frame(contrast=rownames(x),x)
#   colnames(x) <- c('contrast',idx)
#   x
# }

# x <- split(DE_genes$target,DE_genes$contrast)
# y <- split(DAS_genes$target,DAS_genes$contrast)
summary.DDD.vs <- function(x,y,contrast,
                           idx = c('DEonly','DE&DAS','DASonly')){
  num <- lapply(contrast,function(i){
    x <- venn2(x[[i]],y[[i]])
    names(x) <- idx
    x
  })
  names(num) <- contrast
  
  n1 <- lapply(contrast,function(i){
    data.frame(contrast=i,t(sapply(num[[i]],length)))
  })
  n1 <- do.call(rbind,n1)
  colnames(n1) <- c('Contrast',idx)
  
  num <- unlist(num,recursive=FALSE)
  n2 <- sapply(idx,function(i){
    length(Reduce(intersect,num[grep(paste0('\\.',i,'$'),names(num))]))
  })
  n2 <- data.frame(Contrast='Intersection',t(n2))
  names(n2) <- c('Contrast',idx)
  rbind(n1,n2)
}


summary.DDD.number <- function(DE_genes,DAS_genes,DE_trans,DTU_trans,contrast){
  # n1 <- lapply(DE_genes,function(x) x$target)
  idx <- factor(DE_genes$contrast,levels = contrast)
  n1 <- split(DE_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x<- data.frame(`DE genes`=c(sapply(n1,length),`Intersection`=n2))
  
  # n1 <- lapply(DAS_genes,function(x) x$target)
  idx <- factor(DAS_genes$contrast,levels = contrast)
  n1 <- split(DAS_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DAS genes`=c(sapply(n1,length),`Intersection`=n2)))
  
  # n1 <- lapply(DE_trans,function(x) x$target)
  idx <- factor(DE_trans$contrast,levels = contrast)
  n1 <- split(DE_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DE transcripts`=c(sapply(n1,length),`Intersection`=n2)))
  
  
  # n1 <- lapply(DTU_trans,function(x) x$target)
  idx <- factor(DTU_trans$contrast,levels = contrast)
  n1 <- split(DTU_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DTU transcripts`=c(sapply(n1,length),`Intersection`=n2)))
  
  x <- data.frame(contrast=rownames(x),x,row.names = NULL)
  colnames(x) <- c('contrast','DE genes','DAS genes','DE transcripts','DTU transcripts')
  x
}



#############################################################################################
transAbundance2PS <- function(transAbundance=NULL,
                              PS=NULL,
                              contrast,
                              condition,
                              mapping){
  if(is.null(PS)){
    transAbundance.mean <- t(rowmean(t(transAbundance),group = condition,reorder = F))
    ##---PS
    PS <- rowratio(x = transAbundance.mean[as.vector(mapping$TXNAME),],
                   group = as.vector(mapping$GENEID), reorder = F)
    PS <- data.frame(PS[mapping$TXNAME,])
  }
  
  ##---deltaPSI
  deltaPS <- lapply(strsplit(contrast,'-'),function(x){
    PS[,x[1]]-PS[,x[2]]
  })
  deltaPS <- do.call(cbind,deltaPS)
  colnames(deltaPS) <- contrast
  rownames(PS) <- rownames(deltaPS) <- mapping$TXNAME
  list(PS=PS,deltaPS=deltaPS)
}

#############################################################################################
venn2 <- function(x,y){
  a1 <- setdiff(x,y)
  a2 <- intersect(x,y)
  a3 <- setdiff(y,x)
  results <- list(a1=a1,a2=a2,a3=a3)
  attributes(results) <- list(
    a1='x-x&y',
    a2='x&y',
    a3='y-x&y'
  )
  names(results) <- c('a1','a2','a3')
  return(results)
}