#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)


project ="Rshiny-AML"
Sys.setenv(language="en")
if (project == "Rshiny-AML"){
    working_dir<-"C:/Users/jtao/work_dir/Rshiny/test_with_AML_data"
    setwd(working_dir)
}
data.anno <-  read.table( "C:/Users/jtao/work_dir/RShiny/test_with_AML_data/data/AML_GPL570_annotation_subset.txt",  
                          header = TRUE, sep = "\t",check.names = FALSE, row.names = 1)
data <- read.table("C:/Users/jtao/work_dir/RShiny/test_with_AML_data/data/AML_GPL570_matrix_subset.txt", 
                   header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

AML_id <- rownames(data.anno[data.anno$Indication_short %in% c("healthy", "AML"), ])

data.anno <- filter(data.anno, Indication_short %in% c("AML", "healthy"))

data <- data[, colnames(data) %in% AML_id]
# length(AML_id)
# table(data.anno$Indication_short)
data.anno$sample_ID <- rownames(data.anno) 
data.anno <- data.anno[order(factor(data.anno$sample_ID, levels = colnames(data))), ]
data.gene <- rownames(data)
all.equal(colnames(data), rownames(data.anno))


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
    titlePanel("AML data visualization"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width = 3,
             
            fluidRow(
                column(10, selectizeInput("genes",
                                      h3("gene selection", 
                                         style = "font-family: 'arial'; font-si4pt"),
                                      options = list(placeholder = 'Please select up to 5 genes',
                                                     onInitialize = I('function() { this.setValue(""); }')
                                                     ),
                                      selected =TRUE,
                                      # options = c("CD70", "CD33", "IL3RA"),
                            
                                      
                                      choices = "",
                                      multiple = TRUE)
                      
                                       ),
                column(10, selectInput("metadata",
                                       h3("compare data by",
                                          style = "font-family: 'arial'; font-si4pt"),
                                       choices = c("FAB score", "tumor type",
                                                   "mutations"),
                                       
                                       ))
                
                
                
            ),
        #     fluidRow(column(10, sliderInput("bins",
        #                                    "Number of bins:",
        #                                    min = 1,
        #                                    max = 100,
        #                                    value = 30)))
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h1("interactive App for browsing target", align="center",
               style = "font-family: 'times'; font-si16pt"),
            fluidRow(
                plotOutput("boxplot")
            )
           
        )
        
    )
)


# ?mainPanel 
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    updateSelectizeInput(session, 'genes', choices = data.gene, server = TRUE)
    
    
    
    output$boxplot<- renderPlot({
        meta<- reactive({input$metadata})
        if(meta() == "tumor type"){
            
            
            
            data_gene <- reactive({as.data.frame(t(data[rownames(data) %in% input$genes, ]))})
            data_gene_type <- reactive({data_gene() %>% add_column(sample_type = data.anno$Sample_type)})
            data_box <- reactive({data_gene_type() %>% 
                    melt(id.vars = c("sample_type")) %>%
                    filter(variable %in% input$genes) %>% droplevels()})
            color <- "sample_type"}
        
        else if(meta() == "FAB score"){
            
            data_gene <- reactive({as.data.frame(t(data[rownames(data) %in% input$genes, ]))})
            data_gene_type <- reactive({data_gene() %>% add_column(FAB_score = data.anno$FAB_score)})
            data_box <- reactive({data_gene_type() %>% 
                    melt(id.vars = c("FAB_score")) %>%
                    filter(variable %in% input$genes) %>% droplevels()})
            color <- "FAB_score"
            
        }
        else if(meta() == "mutations"){
            data_gene <- reactive({as.data.frame(t(data[rownames(data) %in% input$genes, ]))})
            data_gene_type <- reactive({data_gene() %>% add_column(mutation = data.anno$mutation_status)})
            data_box <- reactive({data_gene_type() %>% 
                    melt(id.vars = c("mutation")) %>%
                    filter(variable %in% input$genes) %>% droplevels()})
            color <- "mutation"
            
        }
        
        maxy <- max(data_box()$value)
        ggplot(data_box(), aes_string(x = "variable",
                  y = "value", color = color )) +
            geom_boxplot() + 
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
                  axis.title = element_text(size = 12)
                  )
            
    })


}
options(shiny.reactlog=TRUE)

# Run the application 
shinyApp(ui = ui, server = server)
# runApp("test_trial", display.mode = "showcase")

data.anno$Sample_type

