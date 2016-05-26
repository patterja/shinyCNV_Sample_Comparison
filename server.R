shinyServer(function(input, output) {

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Depth of Coverage 
output$plot_doc1 <- renderPlot({
  
  ggplot(doc1, aes(x=coverage)) +
    geom_histogram(binwidth=1) +
    geom_rect(data=data.frame(xmin=min(input$depth1), xmax=max(input$depth1)), 
              aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
              color="red", alpha=0.5, inherit.aes = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
    labs(x="depth", y="frequency")
})

output$plot_doc2 <- renderPlot({
  ggplot(doc2, aes(x=coverage)) +
    geom_histogram(binwidth=1) +
    geom_rect(data=data.frame(xmin=min(input$depth2), xmax=max(input$depth2)), 
              aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
              color="red", alpha=0.5, inherit.aes = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
    labs(x="depth", y="frequency")
})

output$plot_doc1_ec <- renderPlot({
  
  ggplot(doc1, aes(x=coverage)) +
    geom_histogram(binwidth=1) +
    geom_rect(data=data.frame(xmin=min(input$depth1_ec), xmax=max(input$depth1_ec)), 
              aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
              color="red", alpha=0.5, inherit.aes = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
    labs(x="genomic location", y="depth")
})

output$plot_doc2_ec <- renderPlot({
  ggplot(doc2, aes(x=coverage)) +
    geom_histogram(binwidth=1) +
    geom_rect(data=data.frame(xmin=min(input$depth2_ec), xmax=max(input$depth2)), 
              aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
              color="red", alpha=0.5, inherit.aes = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
    labs(x="genomic location", y="depth")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ED abbrevation text for help
output$text1 <- renderUI({
  HTML(paste("ED= ExomeDepth","EC=ExomeCNV",sep="<br/>"))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Data for between Tumors
doc1_filt<- reactive({
  makeGRangesFromDataFrame(filter(
    doc1, coverage >= min(input$depth1) & coverage <= max(input$depth1)), 
    keep.extra.columns = TRUE)})

doc2_filt<- reactive({
  makeGRangesFromDataFrame(filter(
    doc2, coverage >= min(input$depth2) & coverage <= max(input$depth2)), 
    keep.extra.columns = TRUE)})

doc1_filt_ec<- reactive({
  makeGRangesFromDataFrame(filter(
    doc1, coverage >= min(input$depth1_ec) & coverage <= max(input$depth1_ec)), 
    keep.extra.columns = TRUE)})

doc2_filt_ec<- reactive({
  makeGRangesFromDataFrame(filter(
    doc2, coverage >= min(input$depth2_ec) & coverage <= max(input$depth2_ec)), 
    keep.extra.columns = TRUE)})

ed1_filt<- reactive({ 
  ed1_int <- ed1  %>% 
    filter(reads.ratio > (1+as.numeric(input$sigfilter))|
                      reads.ratio <= (1-as.numeric(input$sigfilter))) %>%
    filter(!grepl('chrX|chrY', chromosome)) 
  subsetByOverlaps(query=makeGRangesFromDataFrame(ed1_int,keep.extra.columns = TRUE),subject=doc1_filt())})

ed2_filt<- reactive({
  ed2_int<-ed2  %>% 
    filter(reads.ratio > (1+as.numeric(input$sigfilter))|reads.ratio <= (1-as.numeric(input$sigfilter))) %>% 
    filter(!grepl('chrX|chrY', chromosome)) 
  
  subsetByOverlaps(query=makeGRangesFromDataFrame(ed2_int,keep.extra.columns = TRUE),subject=doc2_filt())})

ec1_filt<- reactive({
  ec1_int<-ec1  %>% 
    filter(ratio > (1+as.numeric(input$sigfilter_ec))|ratio <= (1-as.numeric(input$sigfilter_ec))) %>%
    filter(!grepl('chrX|chrY', chr)) 
  
  subsetByOverlaps(query=makeGRangesFromDataFrame(ec1_int,keep.extra.columns = TRUE),subject=doc1_filt_ec())})

ec2_filt<- reactive({
  ec2_int<-ec2  %>% 
    filter(ratio > (1+as.numeric(input$sigfilter_ec))|ratio <= (1-as.numeric(input$sigfilter_ec))) %>%
    filter(!grepl('chrX|chrY', chr)) 
  subsetByOverlaps(query=makeGRangesFromDataFrame(ec2_int,keep.extra.columns = TRUE),subject=doc2_filt_ec())})

ovfilt_tumor_ed <- reactive({
  ov<-olRanges(query=ed1_filt(),subject=ed2_filt(), output="df")
  filter(ov, OLtype %in% input$ovfilter_tumor_ed)})

ovfilt_tumor_ec <- reactive({
  ov<-olRanges(query=ec1_filt(),subject=ec2_filt(), output="df")
  filter(ov, OLtype %in% input$ovfilter_tumor_ec)})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Circos plots

output$circos_tumor_ed <- renderPlot({
  if(input$circdata_tumor_ed=="overlap"){
    ed1gr <- ed1_filt()[ovfilt_tumor_ed()$Qindex]
    seqlevels(ed1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed1gr,force=TRUE)<-seqinfo(hg19sub)
    ed2gr <- ed2_filt()[ovfilt_tumor_ed()$Sindex]
    seqlevels(ed2gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed2gr,force=TRUE)<-seqinfo(hg19sub)
    
  }else{
    ed1gr <- ed1_filt()
    seqlevels(ed1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed1gr,force=TRUE)<-seqinfo(hg19sub)
    ed2gr <- ed2_filt()
    seqlevels(ed2gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ed2gr,force=TRUE)<-seqinfo(hg19sub)}
 
  ggbio(buffer = 5, radius = 120) + 
    circle(hg19sub,geom ="text",aes(label=seqnames,angle=90, hjust=1), vjust=1,size=4, trackwidth=10) +
    circle(hg19sub,geom ="ideo",fill="gray70",trackWidth=2) +
    circle(ed1gr,geom="rect",stat="identity",aes(y="type"),fill="forestgreen",color="forestgreen",trackWidth=30, radius=140) +
    circle(ed2gr,geom="rect",stat="identity",aes(y="type"),fill="dark magenta",color="darkmagenta", trackWidth=30,radius=170)

},height = 600,width = 600)

#ExomeCNV
output$circos_tumor_ec <- renderPlot({
  if(input$circdata_tumor_ec=="overlap"){
    ec1gr <- ec1_filt()[ovfilt_tumor_ec()$Qindex]
    seqlevels(ec1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec1gr,force=TRUE)<-seqinfo(hg19sub)
    ec2gr <- ec2_filt()[ovfilt_tumor_ec()$Sindex]
    seqlevels(ec2gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec2gr,force=TRUE)<-seqinfo(hg19sub)
    
  }else{
    ec1gr <- ec1_filt()
    seqlevels(ec1gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec1gr,force=TRUE)<-seqinfo(hg19sub)
    ec2gr <- ec2_filt()
    seqlevels(ec2gr,force=TRUE)<-seqlevels(hg19sub)
    seqinfo(ec2gr,force=TRUE)<-seqinfo(hg19sub)}
  
  ggbio(buffer = 5, radius = 120) + 
    circle(hg19sub,geom ="text",aes(label=seqnames,angle=90, hjust=1), vjust=1,size=4) +
    circle(hg19sub,geom ="ideo",fill="gray70",trackWidth=2) +
    circle(ec1gr,geom="rect",color="forestgreen", fill="forestgreen",trackWidth=30, radius=140) +
    circle(hg19sub,geom ="ideo",fill="dark green",colour="dark green",trackWidth=0.5,radius=155) +
    circle(ec2gr,geom="rect",color="darkmagenta", fill="darkmagenta", trackWidth=30,radius=170)
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Venn Diagrams
output$depth_venn<-renderD3vennR({
  ov<-ovfilt_tumor_ed()
    venn_tooltip2(
      d3vennR(
        data = list(
          list( sets = list(0), 
                label=paste("Tumor1","ExomeDepth", length(ed1_filt()), sep=",\n"), 
                size=length(ed1_filt()), gene=""),
          list( sets = list(1), 
                label=paste("Tumor2","ExomeDepth", length(ed2_filt()), sep=",\n"), 
                size=length(ed2_filt()), gene=""), 
          list( sets = list(0,1), 
                size=length(unique(ov$Qindex)), 
                gene="")
        ),layoutFunction = '
function(d) { return venn.venn(d, { initialLayout: venn.classicMDSLayout });}
  '
      )
    )
})

output$depth_venn_ec<-renderD3vennR({
  ov<-ovfilt_tumor_ec()
    venn_tooltip2(
      d3vennR(width = 450,
        data = list(
          list( sets = list(0), 
                label=paste("Tumor1", "ExomeCNV", length(ec1_filt()),sep=",\n"), 
                size=length(ec1_filt()), gene=""),
          list( sets = list(1), 
                label=paste("Tumor2","ExomeCNV", length(ec2_filt()),sep=",\n"), 
                size=length(ec2_filt()), gene=""), 
          list( sets = list(0,1), 
                size=length(unique(ov$Qindex)), 
                gene="")
        )
      )
    )
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Stats Tables
output$ed_depthvalues <-renderTable({
  ov<- olRanges(ed1_filt(), ed2_filt(), output = "gr") 
  
  data.frame(        
    Metric = c("CNV within depth range",
               "LOF within depth range",
               "GOF within depth range",
               "Average size within depth range",
               "Std. deviation within depth range",
               "Min size within depth range",
               "Median size within depth range",
               "Max size within depth range",
               "Average Number of bases overlap",
               "Average % overlap"),
    
    Tumor1 = c(paste0(length(ed1_filt())," / ",nrow(ed1)),
               paste0(length(ed1_filt()$reads.ratio[ed1_filt()$type=="deletion"])," / ", 
                      length(ed1$reads.ratio[ed1$type=="deletion"])),
               paste0(length(ed1_filt()$reads.ratio[ed1_filt()$type=="duplication"])," / ",
                      length(ed1$reads.ratio[ed1$type=="duplication"])),
               mean(width(ed1_filt())),
               sd(width(ed1_filt())),
               min(width(ed1_filt())),
               median(width(ed1_filt())),
               max(width(ed1_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercQ)),
    Tumor2 = c(paste0(length(ed2_filt())," / ",nrow(ed2)),
               paste0(length(ed2_filt()$reads.ratio[ed2_filt()$type=="deletion"])," / ", 
                      length(ed2$reads.ratio[ed2$type=="deletion"])),
               paste0(length(ed2_filt()$reads.ratio[ed2_filt()$type=="duplication"])," / ",
                      length(ed2$reads.ratio[ed2$type=="duplication"])),
               mean(width(ed2_filt())),
               sd(width(ed2_filt())),
               min(width(ed2_filt())),
               median(width(ed2_filt())),
               max(width(ed2_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercS))
  )}, digits=3)

output$ec_depthvalues <-renderTable({
  ov<- olRanges(ec1_filt(), ec2_filt(), output = "gr") 
  
  data.frame(        
    Metric = c("CNV within depth range",
               "LOF within depth range",
               "GOF within depth range",
               "Average size within depth range",
               "Std. deviation within depth range",
               "Min size within depth range",
               "Median size within depth range",
               "Max size within depth range",
               "Average Number of bases overlap",
               "Average % overlap"),
    
    Tumor1 = c(paste0(length(ec1_filt())," / ",nrow(ec1)),
               paste0(length(ec1_filt()$ratio[ec1_filt()$ratio<1])," / ", 
                      length(ec1$ratio[ec1$ratio<1])),
               paste0(length(ec1_filt()$ratio[ec1_filt()$ratio>=1])," / ",
                      length(ec1$ratio[ec1$ratio>=1])),
               mean(width(ec1_filt())),
               sd(width(ec1_filt())),
               min(width(ec1_filt())),
               median(width(ec1_filt())),
               max(width(ec1_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercQ)),
    Tumor2 = c(paste0(length(ec2_filt())," / ",nrow(ec2)),
               paste0(length(ec2_filt()$ratio[ec2_filt()$ratio<1])," / ", 
                      length(ec2$ratio[ec2$ratio<1])),
               paste0(length(ec2_filt()$ratio[ec2_filt()$ratio>=1])," / ",
                      length(ec2$ratio[ec2$ratio>=1])),
               mean(width(ec2_filt())),
               sd(width(ec2_filt())),
               min(width(ec2_filt())),
               median(width(ec2_filt())),
               max(width(ec2_filt())),
               mean(ov$OLlength),
               mean(ov$OLpercS))
  )}, digits=3)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Data tables

output$ed_overlap_tumortable <- renderDataTable({
  options(scipen=6, digits=2)
  
  ed1 <- select(as.data.frame(ed1_filt()[ovfilt_tumor_ed()$Qindex,]),seqnames,start,end,type,BF,reads.ratio,gene,cytoband)
  ed2 <- select(as.data.frame(ed2_filt()[ovfilt_tumor_ed()$Sindex,]),seqnames,start,end,type,BF,reads.ratio,gene,cytoband)     
  dafr<-cbind.data.frame(overlapType=ovfilt_tumor_ed()$OLtype,
                         overlapLength=ovfilt_tumor_ed()$OLlength,
                         Tumor1_overlap=as.numeric(ovfilt_tumor_ed()$OLpercQ),
                         Tumor2_overlap=as.numeric(ovfilt_tumor_ed()$OLpercS),
                         Tumor1_chrom=ovfilt_tumor_ed()$space,Tumor1_start=ovfilt_tumor_ed()$Qstart,Tumor1_end=ovfilt_tumor_ed()$Qend,
                         Tumor2_chrom=ovfilt_tumor_ed()$space,Tumor2_start=ovfilt_tumor_ed()$Sstart,Tumor2_end=ovfilt_tumor_ed()$Send,
                         Tumor1_type=ed1_filt()$type[ovfilt_tumor_ed()$Qindex],Tumor2_type=ed2_filt()$type[ovfilt_tumor_ed()$Sindex],
                         Tumor1_reads.ratio=ed1_filt()$reads.ratio[ovfilt_tumor_ed()$Qindex],
                         Tumor2_reads.ratio=ed2_filt()$reads.ratio[ovfilt_tumor_ed()$Sindex], 
                         gene=ed1$gene[ovfilt_tumor_ed()$Qindex],
                         cytoband=ed1_filt()$cytoband[ovfilt_tumor_ed()$Qindex])
})


output$ed_setdiff_tab12<-renderDataTable({
  as.data.frame(GenomicRanges::setdiff(ed1_filt(), ed2_filt()))
})


output$ed_setdiff_tab21<-renderDataTable({
  as.data.frame(GenomicRanges::setdiff(ed2_filt(),ed1_filt()))
})

output$ed1_tumortable <- renderDataTable({
  select(as.data.frame(ed1_filt()),"Tumor1_chrom"=seqnames,"Tumor1_start"=start,
         "Tumor1_end"=end,"Tumor1_type"=type,"Tumor1_nexons"=nexons,"Tumor1_BF"=BF,
         "Tumor1_reads.ratio"=reads.ratio,"Tumor1_gene"=gene, "Tumor1_cytoband"=cytoband)})
output$ed2_tumortable <- renderDataTable({
  select(as.data.frame(ed2_filt()),"Tumor2_chrom"=seqnames,"Tumor2_start"=start,
         "Tumor2_end"=end,"Tumor2_type"=type,"Tumor2_nexons"=nexons,"Tumor2_BF"=BF,
         "Tumor2_reads.ratio"=reads.ratio,"Tumor2_gene"=gene, "Tumor2_cytoband"=cytoband)})

output$ec_overlap_tumortable <- renderDataTable({
  options(scipen=6, digits=2)
  
  dafr<-cbind.data.frame(overlapType=ovfilt_tumor_ec()$OLtype,
                         overlapLength=ovfilt_tumor_ec()$OLlength,
                         Tumor1_overlap=as.numeric(ovfilt_tumor_ec()$OLpercQ),
                         Tumor2_overlap=as.numeric(ovfilt_tumor_ec()$OLpercS),
                         Tumor1_chrom=ovfilt_tumor_ec()$space,Tumor1_start=ovfilt_tumor_ec()$Qstart,Tumor1_end=ovfilt_tumor_ec()$Qend,
                         Tumor2_chrom=ovfilt_tumor_ec()$space,Tumor2_start=ovfilt_tumor_ec()$Sstart,Tumor2_end=ovfilt_tumor_ec()$Send,
                         Tumor1_type=ec1_filt()$type[ovfilt_tumor_ec()$Qindex],Tumor2_type=ec2_filt()$type[ovfilt_tumor_ec()$Sindex],
                         Tumor1_ratio=ec1_filt()$ratio[ovfilt_tumor_ec()$Qindex],
                         Tumor2_ratio=ec2_filt()$ratio[ovfilt_tumor_ec()$Sindex], 
                         #gene=ec1_filt()$gene[ovfilt_tumor_ec()$Qindex],
                         cytoband=ec1_filt()$ovCNVec.cytoband[ovfilt_tumor_ec()$Qindex])
})

output$ec_setdiff_tab12<-renderDataTable({
  as.data.frame(GenomicRanges::setdiff(ec1_filt(), ec2_filt()))
})


output$ec_setdiff_tab21<-renderDataTable({
  as.data.frame(GenomicRanges::setdiff(ec2_filt(),ec1_filt()))
})

output$ec1_tumortable <- renderDataTable({
  select(as.data.frame(ec1_filt()),"Tumor1_chrom"=seqnames,"Tumor1_start"=start,
         "Tumor1_end"=end,"Tumor1_logR"=logR,"Tumor1_ratio"=ratio,"Tumor1_avg.coverage"=average.coverage, "cytoband"=ovCNVec.cytoband)})

output$ec2_tumortable <- renderDataTable({
  select(as.data.frame(ec2_filt()),"Tumor2_chrom"=seqnames,"Tumor2_start"=start,
         "Tumor2_end"=end,"Tumor2_logR"=logR,"Tumor2_ratio"=ratio,"Tumor2_avg.coverage"=average.coverage, "cytoband"=ovCNVec.cytoband)})
#-----------------------------------------------------IDR tab

ed_idrout<- reactive({
  sig <- cbind(((X1G6_ed$BF[ov_ed$Qindex]-min(X1G6_ed$BF[ov_ed$Qindex]))/
                  (max(X1G6_ed$BF[ov_ed$Qindex])-min(X1G6_ed$BF[ov_ed$Qindex]))),
               ((X1G8_ed$BF[ov_ed$Sindex]-min(X1G8_ed$BF[ov_ed$Sindex]))/
                  (max(X1G8_ed$BF[ov_ed$Sindex])-min(X1G8_ed$BF[ov_ed$Sindex]))))
  idr.out <- est.IDR(x=sig, mu=mean(c(sig[,1], sig[,2])),
                     sigma=sd(c(sig[,1], sig[,2])), 
                     rho=cor(sig[,1], sig[,2]), 
                     p=0.6, eps=0.001, max.ite=20)
})

idrfilt<- reactive({
  idrfilt<- data.frame(cbind(IDR=ed_idrout()$IDR, ov_ed, ed1=ed1[ov_ed$Qindex,], ed2=ed2[ov_ed$Sindex,]))
  idrfilt %>%
    filter(IDR <= input$idrthresh)
})



})