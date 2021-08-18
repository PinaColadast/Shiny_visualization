#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#check library dependencies 




list.of.packages <- c("shiny","ggplot2","ggpubr",
                      "dplyr","reshape2","tidyverse","comprehenr","ggrepel",
                      "RColorBrewer")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
# if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)


library(shiny)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(comprehenr)
library(ggrepel)
library(RColorBrewer)



project ="Rshiny-AML"
Sys.setenv(language="en")

data.anno <-  read.table( "./data/AML_GPL570_annotation_AML_healthy.txt",  
                                                    header = TRUE, sep = "\t",check.names = FALSE, row.names = 1)

data <- read.table("./data/AML_GPL570_matrix_AML_healthy.txt",
                                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

data.anno$sample_ID <- rownames(data.anno) 
data.anno <- data.anno[order(factor(data.anno$sample_ID, levels = colnames(data))), ]
data.gene <- rownames(data)
all.equal(colnames(data), rownames(data.anno))
data.anno[data.anno$FAB_score == "", ] <-"no information" 

mutation_cate <- grep("mutation", colnames(data.anno), value =TRUE)
return_mutation<- function(x1, x2){
    # print(length(x1))
    if (is.na(x1) | is.na(x2)){
        mut = "no info"
    }
    else if (x1==x2){
        mut = x1}
    
    else if (x1!=x2){
        mut = paste(x1, ", ", x2)
    }
    return(mut)
}



data.anno$mutation_status1 <- apply(data.anno[, mutation_cate],
                                    1,
                                    FUN = function(x) min(names(which(x=="positive",
                                    ))))

data.anno$mutation_status2 <- apply(data.anno[, mutation_cate],
                                    1,
                                    FUN = function(x) max(names(which(x=="positive",
                                    ))))

data.anno$mutation_status <- apply(data.anno[, c("mutation_status1",
                                                 "mutation_status2")], 1,
                                   function(y) return_mutation(y[1], y[2]))


# data.anno$FAB_score
# data_bind <- rbind(data, data.anno$Sample_type)
# rownames(data_bind) <- c(rownames(data), "Sample type")
# if (all.equal(colnames(data), rownames(data.anno))){
#     data[nrow(data) +1, ] = data.anno$Sample_type
# }

# genes <- list(data$`Gene Symbol`)
# Define UI for application that draws a histogram

give.n <- function(x){
    return(data.frame(y = 1.3*(median(x)), label = paste0("n=", length(x))))
}


ui <- fluidPage(

    # Application title
    titlePanel(
        p(
          h1("Acute Myeloid Leukemia data visualization"),
          h3("Data description:", style = "font-family: 'arial'; font-size: 20px"),
          h3("Microarray Bulk RNA-Seq data of bone marrow & PBMCs samples from 682 AML patients and 74 healthy donors",
             style = "font-family: 'arial'; font-size: 20px")
          )
        ),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width = 3,
             
            fluidRow(
                column(10, selectInput("sample_origin",
                                          h3("SAMPLE ORIGIN",
                                             style = "font-family: 'arial'; font-size: 16px"),
                                       choices = c("all", "bone marrow", "peripheral blood"),
                                       selected = "all",
                                       
                                      )
            
                                             ),
                column(10, selectizeInput("genes",
                                      h3("GENE SELECTION", 
                                         style = "font-family: 'arial'; font-size: 16px"),
                                      options = list(placeholder = 'Please select up to 5 genes',
                                                     maxOptions = 5)
                                                     ,
                                      choices = NULL,
                                      multiple = TRUE)
                      
                                       ),
                column(10, selectInput("metadata",
                                       h3("VISUALIZE DATA BY METADATA",
                                          style = "font-family: 'arial'; font-size: 16px"),
                                       choices = c("FAB score", "tumor type",
                                                   "mutations"),
                                       
                                       ))
                
                
                
            ),
    
        ),

  
        mainPanel(
          tabsetPanel(
            
            
            tabPanel("Sample Distribution",  plotOutput("Distribution",
                                                        width = "90%",
                                                        height = "420px"
            )),
            tabPanel("Gene Expression", br(), plotOutput("TargetExpre",
                                                         width = "100%",
                                                        height = "400px"))
            
            
          )
            
            

        )

        
        
    )
)


# ?mainPanel 
# Define server logic required to draw a histogram
server <-function(input, output, session) {
    updateSelectizeInput(session, 'genes', choices = data.gene,
                                                       selected = c("CD70"), server = TRUE,
                         options = list(render = I(
                             '{
    option: function(item, escape) {
    return "<div>" + escape(item.value) + "</div>";
      
    }
  }'
                             
                         )))
    
    # browser()
    data_box <- reactive({
        if (length(input$genes)!=0 ){
          inputgenes = input$genes
          
        }else{
          inputgenes = c("CD70")
        }
      
        if (input$sample_origin!="all"){
            data.anno <- filter(data.anno, Sample_tissue_of_origin == input$sample_origin)
            data <- data[,rownames(data.anno)]
        }else{
            data.anno <- data.anno
            data <- data
        }
        
        data_gene <- as.data.frame(t(data[rownames(data) %in% inputgenes, ]))
        
        if(input$metadata == "tumor type"){
            data_gene_type <- data_gene %>% add_column(sample_type = data.anno$Sample_type)
            data_box <- data_gene_type %>% 
                    melt(id.vars = c("sample_type")) %>%
                    filter(variable %in% inputgenes) %>% droplevels()
            color <- "sample_type"}
        
        else if(input$metadata == "FAB score"){
            # browser()
            
            data_gene_type <- data_gene%>% add_column(FAB_score = data.anno$FAB_score)
            data_box <- data_gene_type %>% 
                    melt(id.vars = c("FAB_score")) %>%
                    filter(variable %in% inputgenes) %>% droplevels()
            color <- "FAB_score"
            
        }
        else if(input$metadata == "mutations"){
            
            data_gene_type <- data_gene %>% add_column(mutation = data.anno$mutation_status)
            data_box <- data_gene_type %>% 
                    melt(id.vars = c("mutation")) %>%
                    filter(variable %in% inputgenes) %>% droplevels()
            color <- "mutation"
            
        }
        return (list("databox"=data_box, "color" = color,
                     "genes" = inputgenes))
    })
    
    
    #STOPIFNOT
    
    
    output$TargetExpre<- renderPlot({
        
        # stopifnot(data_box()$databox)
        data_box<- data_box()$databox
        color <- data_box()$color
        
        
        data_box <- data_box %>%
          mutate(variable = str_replace_all(variable, "_", " "))
        
        maxy <- max(data_box$value)
        
        ggplot(data_box, aes_string(x = "variable",
                  y = "value", fill = color)) +
            geom_boxplot() +
            scale_color_brewer(palette="Set3")+
        stat_compare_means(
            # group.by = color,
            aes_string(group = color),
                           label = "p.signif", label.y = maxy) + 
            stat_summary(fun.data = give.n, geom = "text",
                         position = position_dodge(width = 0.75)) +
            labs(x = "Genes", y = "RNA expression") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text=element_text(size = 11 ),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 11)
                  )
            
    })
    
    output$Distribution<- renderPlot({
      
        var_width = 5
        
    
        data_box<- data_box()$databox %>% filter(variable == data_box()$genes[1])
        category <- data_box()$color
        data_bar <- as.data.frame(table(data_box[category]))
        names(data_bar) <- c("Category", "Sample Numbers")
        
        if (max(nchar(as.vector(data_bar$Category))) > 14){
          angle = 90
        } else{
          angle = 0
        } 
        
        
        
        data_bar <- data_bar %>%
          mutate(Category = str_wrap(Category, width = var_width)) %>%
          mutate(Category = str_replace_all(Category, "_", " "))
        barplot <- ggbarplot(data_bar,x = "Category", y = "Sample Numbers",
                  fill = "Category", palette = brewer.pal(n = 8
                                                          , name = "Set3"),
                  label = TRUE, legend = "right")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text=element_text(size = 11,
                                         angle = angle),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 11)
            )
            # stat_summary(fun.data = give.n.bar, geom = "text",
            #              position = position_dodge(width = 0.75))
         
        barplot 
        
        


    })


}
options(shiny.reactlog=TRUE)

# Run the application 
shinyApp(ui = ui, server = server)
# runApp("C:/Users/jtao/work_dir/RShiny/test_with_AML_data/app.R")

