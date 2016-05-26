
shinyUI(fluidPage(
h3("Sample Comparison of Biological Replicates"),
#----------------------------------------------------------------Between Tumors
tabsetPanel(
  tabPanel("ExomeDepth",h5("ExomeDepth [Plagnol,et al, 2012]"),
         fluidRow(  
           column(6,
                  plotOutput("circos_tumor_ed")),
           column(3,
                  h4("Depth of Coverage Tumor1", style="color:darkgreen"),
                  plotOutput("plot_doc1",height = 200),
                  sliderInput('depth1', 'Minimum Depth of Coverage threshold', 
                              min = min(doc1$coverage, na.rm=TRUE), 
                              max = max(doc1$coverage, na.rm=TRUE), 
                              value = c(min(doc1$coverage, na.rm=TRUE),max(doc1$coverage, na.rm=TRUE))),
                  submitButton("UPDATE"),
                  hr(),
                  radioButtons("sigfilter", label = h4("Filter significance"),
                               choices = list("logR" = 0.5, "No Filter" = 0), 
                               selected = 0),
                  submitButton("UPDATE")),  
           column(3,
                  h4("Depth of Coverage Tumor2", style="color:darkmagenta"),
                  plotOutput("plot_doc2", height=200),
                  sliderInput('depth2', 'Minimum Depth of Coverage threshold', 
                              min = min(doc2$coverage, na.rm=TRUE), 
                              max = max(doc2$coverage, na.rm=TRUE), 
                              value = c(min(doc2$coverage, na.rm=TRUE),max(doc2$coverage, na.rm=TRUE))),
                  submitButton("UPDATE"),
                  hr(),
                  radioButtons("circdata_tumor_ed", label = h4("Overlap Filter"),
                               choices = list("All Results"="allresults", "Overlaps Only" = "overlap"), 
                               selected= "allresults"),
                  htmlOutput("text1"),
                  conditionalPanel(condition = "input.circdata_tumor_ed == 'overlap'",
                     checkboxGroupInput("ovfilter_tumor_ed",label = h4("Filter Overlaps by Type"),
                          choices=list("ED contained in EC"="contained",
                                       "EC inside ED"="inside",
                                       "overlapdown of ED"="oldown",
                                       "overlapup of ED"="olup"),
                          selected=c("contained","inside","oldown", "olup"))),
                  submitButton("UPDATE"))),
           fluidRow(
             
             column(5,h4("ExomeDepth Tumor1&2"),h5("venn overlap if any ranges overlap"),d3vennROutput("depth_venn", width=400, height=400)),
             column(6,h4("ExomeDepth Tumor 1 & 2"),tableOutput("ed_depthvalues"))
             ),
           
           
           tabsetPanel(
             tabPanel("Overlap",dataTableOutput("ed_overlap_tumortable")),
             tabPanel("Tumor1 unique", dataTableOutput("ed_setdiff_tab12")),
             tabPanel("Tumor2 unique", dataTableOutput("ed_setdiff_tab21")),
             tabPanel("ExomeDepth Tumor1", dataTableOutput("ed1_tumortable")),
             tabPanel("ExomeDepth Tumor2", dataTableOutput("ed2_tumortable"))
           )),
  

# #----------------------------------------------------------------
   tabPanel("ExomeCNV",h5("ExomeCNV [Sathirapongsasuti,etal, 2011]"),
    fluidRow(  
      column(6, plotOutput("circos_tumor_ec")),
      column(3, h4("Depth of Coverage Tumor1", style="color:darkgreen"),
            plotOutput("plot_doc1_ec",height = 200),
            sliderInput('depth1_ec', 'Minimum Depth of Coverage threshold', 
                        min = min(doc1$coverage, na.rm=TRUE), 
                        max = max(doc1$coverage, na.rm=TRUE), 
                        value = c(0, 300)),
            submitButton("UPDATE"),
              hr(),
              radioButtons("sigfilter_ec", label = h4("Filter significance"),
                         choices = list("logR" = 0.5, "No Filter" = 0), 
                         selected = 0),
             submitButton("UPDATE")),  
      column(3,h4("Depth of Coverage Tumor2", style="color:darkmagenta"),
            plotOutput("plot_doc2_ec", height=200),
            sliderInput('depth2_ec', 'Minimum Depth of Coverage threshold', 
                        min = min(doc2$coverage, na.rm=TRUE), 
                        max = max(doc2$coverage, na.rm=TRUE), 
                        value = c(min(doc2$coverage, na.rm=TRUE),max(doc2$coverage, na.rm=TRUE))),
            submitButton("UPDATE"),
              hr(),
              radioButtons("circdata_tumor_ec", label = h4("Overlap Filter"),
                         choices = list("All Results"="allresults", "Overlaps Only" = "overlap"), 
                         selected= "allresults"),
            conditionalPanel(condition = "input.circdata_tumor_ec == 'overlap'",
                             checkboxGroupInput("ovfilter_tumor_ec",label = h4("Filter Overlaps by Type"),
                                                choices=list("Tumor1 contained in Tumor2"="contained",
                                                             "Tumor1 inside"="inside",
                                                             "overlapdown"="oldown",
                                                             "overlapup"="olup"),
                                                selected=c("contained","inside","oldown", "olup"))),
             submitButton("UPDATE"))),
   fluidRow(
     column(5,h4("ExomeCNV Tumor Overlap"),h5("venn overlap if any ranges overlap"),d3vennROutput("depth_venn_ec"), offset=1),
     column(6,h4("ExomeCNV Tumor1 & 2"),tableOutput("ec_depthvalues")),
     
     hr()),
   
   
   tabsetPanel(
     tabPanel("Overlap",dataTableOutput("ec_overlap_tumortable")),
     tabPanel("Tumor1 unique", dataTableOutput("ec_setdiff_tab12")),
     tabPanel("Tumor2 unique", dataTableOutput("ec_setdiff_tab21")),
     tabPanel("ExomeCNV Tumor1", dataTableOutput("ec1_tumortable")),
     tabPanel("ExomeCNV Tumor2", dataTableOutput("ec2_tumortable")))
   )

)))