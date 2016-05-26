library(dplyr)
library(tidyr)
library(ggvis)
library(ggbio)
library(GenomicRanges)
library(d3vennR)
library(idr)
source("rangeoverlapper.R")
#credit for range overlapper to Girke Lab
#http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/rangeoverlapper.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Prep ggbio hg19 circos files
data("CRC", package = "biovizBase")
seqlevels(hg19sub) <- paste0("chr", seqlevels(hg19sub))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Depth of Coverage files

doc_samp1 <-read.delim("data/1G6_merged_sort.coverage.sample_interval_summary",header= TRUE, sep="\t", stringsAsFactors = FALSE)


doc1 <- doc_samp1 %>%
  separate(Target, into=c("chrom", "start", "end"), sep=":|-") %>%
  mutate(start=as.numeric(start), end=as.numeric(end), 
         rang = end-start+1, coverage = average_coverage) %>%
  select(chrom, start, end, rang, coverage)

#doc_samp2 <- read.delim("C:/Users/Owner/Box Sync/DepthofCoverage/raw_output_lessparam/1G8_merged_sort.coverage.sample_interval_summary", header= TRUE, sep="\t", stringsAsFactors = FALSE)
doc_samp2 <-read.delim("data/1G8_merged_sort.coverage.sample_interval_summary",header= TRUE, sep="\t", stringsAsFactors = FALSE)

doc2 <- doc_samp2 %>%
  separate(Target, into=c("chrom", "start", "end"), sep=":|-") %>%
  mutate(start=as.numeric(start), end=as.numeric(end),
         rang = end-start+1, coverage = average_coverage) %>%
  select(chrom, start, end, rang, coverage)

#------------------------------------------------------------------------------------
#ExomeDepthCNV results
# X1G6_ed <- read.delim("C:/Users/Owner/Box Sync/Analysis_CNV/ExomeDepthCNV/results_inctransthresh/1G6_exomedepth_annotated_Conrad_somaticCNV.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
# 
# X1G8_ed <- read.delim("C:/Users/Owner/Box Sync/Analysis_CNV/ExomeDepthCNV/results_inctransthresh/1G8_exomedepth_annotated_Conrad_somaticCNV.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
# X1G6_ec <- read.delim("C:/Users/Owner/Box Sync/Analysis_CNV/ExomeCNV/results/X1G6.coverage.cnv.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
# X1G8_ec <- read.delim("C:/Users/Owner/Box Sync/Analysis_CNV/ExomeCNV/results/X1G8.coverage.cnv.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
X1G6_ed <- read.delim("data/1G6_exomedepth_annotated_Conrad_somaticCNV.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)

X1G8_ed <- read.delim("data/1G8_exomedepth_annotated_Conrad_somaticCNV.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
X1G6_ec <- read.delim("data/X1G6.coverage.cnv.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)
X1G8_ec <- read.delim("data/X1G8.coverage.cnv.txt", header= TRUE, sep="\t", stringsAsFactors = FALSE)

ed1 <- X1G6_ed %>%
  select(chromosome, start, end, width=end-start+1, type, nexons, BF, reads.ratio, gene=exons.hg19, cytoband)

ed2 <- X1G8_ed %>%
  select(chromosome, start, end, width=end-start+1, type, nexons, BF, reads.ratio, gene=exons.hg19, cytoband)



#ExomeCNV results
ec1 <- ec_cnv_anno %>%
  filter(sample=="X1G6") %>%
  mutate(type=ifelse(ratio>=1,"duplication", ifelse(ratio<1,"deletion", "")))
ec2 <- ec_cnv_anno %>%
  filter(sample=="X1G8") %>%
  mutate(type=ifelse(ratio>=1,"duplication", ifelse(ratio<1,"deletion", "")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ov_ed <-olRanges(query=makeGRangesFromDataFrame(ed1,seqnames.field ="chromosome",start.field="start", end.field="end"), 
                 subject=makeGRangesFromDataFrame(ed2,seqnames.field ="chromosome",start.field="start", end.field="end"), 
                 output="df")
ov_ec <-olRanges(query=makeGRangesFromDataFrame(ec1,start.field="probe_start", end.field="probe_end"), 
                 subject=makeGRangesFromDataFrame(ec2,start.field="probe_start", end.field="probe_end"), 
                 output="df")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Functions

#------------------------------------------------------------------------------------
venn_tooltip2<- function( venn ){
  venn$x$tasks[length(venn$x$tasks)+1] <- list(
    htmlwidgets::JS('
                    function(){
                    console.log("here");
                    var div = d3.select(this);
                    
                    // add a tooltip
                    var tooltip = d3.select("body").append("div")
                    .attr("class", "venntooltip")
                    .style("position", "absolute")
                    .style("text-align", "center")
                    .style("width", 250)
                    .style("height", 500)
                    .style("background", "#333")
                    .style("color","#ddd")
                    .style("padding","2px")
                    .style("border","0px")
                    .style("border-radius","8px")
                    .style("opacity",0);
                    
                    div.selectAll("path")
                    .style("stroke-opacity", 0)
                    .style("stroke", "#fff")
                    .style("stroke-width", 0)
                    
                    // add listeners to all the groups to display tooltip on mousover
                    div.selectAll("g")
                    .on("mouseover", function(d, i) {
                    
                    // sort all the areas relative to the current item
                    venn.sortAreas(div, d);
                    
                    // Display a tooltip with the current size
                    tooltip.transition().duration(400).style("opacity", .9);
                    tooltip.html(d.size+"<br>"+d.gene);
                    
                    // highlight the current path
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 3)
                    .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                    .style("stroke-opacity", 1);
                    })
                    
                    .on("mousemove", function() {
                    tooltip.style("left", (d3.event.pageX) + "px")
                    .style("top", (d3.event.pageY - 28) + "px");
                    })
                    
                    .on("mouseout", function(d, i) {
                    tooltip.transition().duration(400).style("opacity", 0);
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 0)
                    .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                    .style("stroke-opacity", 0);
                    });
                    }
                    ')
    )
  venn
  }
