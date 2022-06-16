library(shiny)
library(shinythemes)
library(dendextend)
library(RColorBrewer)
library(pals)
library(pheatmap)
library(phylogram)
library(DT)
library(rhandsontable)
library(msaR)
library(grid)
# TODO: Change putative effector MSA file path

# Define the UI for the app
ui <- navbarPage('FoEC2', theme = shinytheme("sandstone"),
                 tabPanel("Data",
                          sidebarLayout(
                            sidebarPanel(
                              h3("Input files"),
                              helpText("The PAV (presence absence variation) table provides information about which putative effectors are found in which genomes. TSV file."),
                              fileInput("pavFile", "PAV table",
                                        accept = c('.tsv')
                              ),
                              helpText("The genome metadata table can be used to add information about your genomes. (i.e. formae speciales). CSV file."),
                              fileInput("metaFile", "Genome metadata table",
                                        accept = c('.csv')
                              ),
                              helpText("The putative effector metadata table can be used to add information about the detected putative effectors (i.e. SIX genes). CSV file."),
                              fileInput("emetaFile", "Putative effector metadata table",
                                        accept = c('.csv')
                              ),
                            ),
                            mainPanel(
                              fluidRow(
                                column(11,
                                       tabsetPanel(
                                         tabPanel('PAV table',
                                                  fluidRow(
                                                    DT::dataTableOutput('pavTable')
                                                  )
                                         ),
                                         tabPanel('Metadata: genomes',
                                                  fluidRow(
                                                    helpText("You can edit your metadata here by following these steps:"),
                                                    helpText(HTML('&thinsp;'), "1. Right click on the table and select 'Insert column left/right'"),
                                                    helpText(HTML('&thinsp;'), "2. Fill in the new cells with data. (Note: column names cannot be modified here. To do so, an external editor such as Excel must be used.)"),
                                                    helpText(HTML('&thinsp;'), "3. Save your changes and send them to the heatmap by pressing 'Update plot'."),
                                                    helpText(HTML('&thinsp;'), "4. Download your changed file by providing a filename and then press 'Download table'.")
                                                  ),
                                                  fluidRow(
                                                    column(6,
                                                           textInput("filename", "Save as", value = paste0("metadata-", Sys.Date(),".csv"))
                                                    ),
                                                    column(3,
                                                           actionButton('updatePlot', 'Update plot')
                                                    ),
                                                    column(3,
                                                           downloadButton('downloadTable', 'Download table')
                                                    )
                                                  ),
                                                  fluidRow(
                                                    rHandsontableOutput('metaTable')
                                                  )
                                         ),
                                         tabPanel('Metadata: putative effectors',
                                                  fluidRow(
                                                    helpText("You can edit your metadata here by following these steps:"),
                                                    helpText(HTML('&thinsp;'), "1. Right click on the table and select 'Insert column left/right'"),
                                                    helpText(HTML('&thinsp;'), "2. Fill in the new cells with data. (Note: column names cannot be modified here. To do so, an external editor such as Excel must be used.)"),
                                                    helpText(HTML('&thinsp;'), "3. Save your changes and send them to the heatmap by pressing 'Update plot'."),
                                                    helpText(HTML('&thinsp;'), "4. Download your changed file by providing a filename and then press 'Download table'.")
                                                  ),
                                                  fluidRow(
                                                    column(6,
                                                           textInput("efilename", "Save as", value = paste0("putative-effector-metadata-", Sys.Date(),".csv"))
                                                    ),
                                                    column(3,
                                                           actionButton('eupdatePlot', 'Update plot')
                                                    ),
                                                    column(3,
                                                           downloadButton('edownloadTable', 'Download table')
                                                    )
                                                  ),
                                                  fluidRow(
                                                    rHandsontableOutput('emetaTable')
                                                  )
                                         )
                                       )
                                )
                              )
                            )
                          )
                 ),
                 navbarMenu('Plots',
                            tabPanel('Heatmap',
                                     sidebarLayout(
                                       sidebarPanel(id = "heatOpts", style = "overflow-y:scroll; max-height: 600px; position:relative;",
                                         h3("Options"),
                                         hr(style ='border-top: 1px solid #000000;'),
                                         textInput("heatname", "Download PDF", value = paste0("heatmap-", Sys.Date(),".pdf")),
                                         downloadButton('downloadHeat', 'Download as PDF'),
                                         textInput("reorderedName", "Download reordered CSV", value = paste0("reordered-", Sys.Date(),".csv")),
                                         downloadButton('downloadReorder', 'Download as CSV'),
                                         helpText('Download a CSV with reordered rows and columns based on clustering methods applied.'),
                                         h4("Genomes"),
                                         selectInput('dist_genomes', h5("Distance method"),
                                                     choices = list("Binary" = 'binary', "Euclidean" = 'euclidean', "Maximum" = 'maximum',
                                                                    "Manhattan" = 'manhattan', "Canberra" = 'canberra', "Minkowski" = 'minkowski'), selected = 'binary'),
                                         selectInput('clust_genomes', h5('Clustering method'),
                                                     choices = list("Average" = 'average', "Single" = 'single', "Complete" = 'complete',
                                                                    "McQuitty" = 'mcquitty', "Median" = 'median', "Centroid" = 'centroid', 'Ward.D' = 'ward.D', 'Ward.D2' = 'ward.D2'), selected = 'average'),
                                         h4("Putative effectors"),
                                         selectInput('dist_effs', h5("Distance method"),
                                                     choices = list("Binary" = 'binary', "Euclidean" = 'euclidean', "Maximum" = 'maximum',
                                                                    "Manhattan" = 'manhattan', "Canberra" = 'canberra', "Minkowski" = 'minkowski'), selected = 'binary'),
                                         selectInput('clust_effs', h5('Clustering method'),
                                                     choices = list("Average" = 'average', "Single" = 'single', "Complete" = 'complete',
                                                                    "McQuitty" = 'mcquitty', "Median" = 'median', "Centroid" = 'centroid', 'Ward.D' = 'ward.D', 'Ward.D2' = 'ward.D2'), selected = 'average'),
                                         h4("Annotation label"),
                                         helpText("To access labels, add information to the metadata table (Data > Metadata: [genomes | putative effectors])"),
                                         selectInput("catSelect", h5("Genome label"),
                                                     choices = ""),
                                         selectInput("ecatSelect", h5("Putative effector label"),
                                                     choices = ""),
                                         h4("Misc"),
                                         sliderInput('heatHeight', "Plot height (px)", 1, 2500, 1000),
                                         sliderInput('heatWidth', "Plot width (px)", 1, 2500, 1000),
                                         checkboxInput('ccheck', label = "Column labels", value = TRUE),
                                         checkboxInput('rcheck', label = "Row labels", value = TRUE),
                                         sliderInput('cfontsize', 'Column font size', 1, 20, 5),
                                         sliderInput('rfontsize', 'Row font size', 1, 20, 10),
                                         textInput('ncuts', 'Number of breaks', value = 1),
                                         selectInput("colors", h5("Colors"),
                                                     choices = list("Blues" = 'Blues', "Greens" = 'Greens', "Purples" = "Purples", "Reds" = "Reds",
                                                                    "Oranges" = "Oranges", "Greys" = "Greys", "Blue-Green" = "BuGn", "Blue-Purple" = "BuPu",
                                                                    "Green-Blue" = "GnBu"), selected = 'Blues')
                                       ),
                                       mainPanel(
                                         fluidRow(
                                           column(11,
                                                  fluidRow(
                                                    uiOutput('heatPlotsized')
                                                  )
                                           )
                                         )
                                       )
                                     )
                            ),
                            tabPanel('Dendrograms',
                                     sidebarLayout(
                                       sidebarPanel(
                                         h3('Options'),
                                         hr(style ='border-top: 1px solid #000000;'),
                                         h4('Genomes'),
                                         radioButtons(
                                           'ghoriz', 'Orientation', choices = c('Vertical' = FALSE, 'Horizontal' = TRUE)
                                         ),
                                         sliderInput('gtxtSize', 'Text size', value = 0.5, min = 0.01, max = 1),
                                         downloadButton('downloadGdend', 'Download Newick'),
                                         h4('Putative effectors'),
                                         radioButtons(
                                           'ehoriz', 'Orientation', choices = c('Vertical' = FALSE, 'Horizontal' = TRUE)
                                         ),
                                         sliderInput('etxtSize', 'Text size', value = 0.5, min = 0.01, max = 1),
                                         downloadButton('downloadEdend', 'Download Newick')
                                       ),
                                       mainPanel(
                                         fluidRow(plotOutput('gCluster')),
                                         fluidRow(plotOutput('eCluster'))
                                       )
                                     )
                            ),
                            tabPanel('MSA',
                                     fluidRow(

                                     ),
                                     fluidRow(
                                       selectInput('msaFile', h5('MSA file'), choices = list.files(path = '../output/03.presenceabsence/putative_effector_msas'))
                                     ),
                                     fluidRow(msaROutput('msaPlot')
                                     )
                            )
                 )
)
server <- function(input, output, session) {
  # Initialize metadata table
  metaData <- reactiveValues(mdata = NULL)
  emetaData <- reactiveValues(edata = NULL)
  # Get PAV file and put in table
  pav <- reactive({
    infile <- input$pavFile
    req(infile)
    d <- read.table(infile$datapath, sep = "\t", header = TRUE)
    row.names(d) <- d[,1] #rename rows to values in first column
    d[,1] <- NULL #remove the first column
    data <- as.matrix(d)
    return (data)
  })
  # Create PAV table output
  output$pavTable <- DT::renderDataTable({
    pav()
  }, rownames = TRUE, options=list(scrollX=TRUE))
  # Put metadata file in a table
  observe({
    metaInfile <- input$metaFile
    if (is.null(metaInfile)) {
      return(NULL)
    }
    metaConfig <- read.table(metaInfile$datapath, sep = ",", header = TRUE)
    row.names(metaConfig) <- metaConfig[,1] #rename rows to values in first column
    if (ncol(metaConfig) <= 1) {
      colnames(metaConfig)[1] <- 'V1'
      metaConfig[,1] <- ''
    } else {
      colnames(metaConfig)[1] <- ''
      metaConfig[,1] <- NULL
    }
    metaData$mdata <- metaConfig
  })
  # Put effector metadata file in a table
  observe({
    emetaInfile <- input$emetaFile
    if (is.null(emetaInfile)) {
      return(NULL)
    }
    emetaConfig <- read.table(emetaInfile$datapath, sep = ',', header = TRUE)
    row.names(emetaConfig) <- emetaConfig[,1] #rename rows to values in first column
    if (ncol(emetaConfig) <= 1) {
      colnames(emetaConfig)[1] <- 'V1'
      emetaConfig[,1] <- ''
    } else {
      colnames(emetaConfig)[1] <- ''
      emetaConfig[,1] <- NULL
    }
    emetaData$edata <- emetaConfig
  })
  # Genome metadata table
  output$metaTable <- renderRHandsontable({
    req(metaData$mdata)
    if (!is.null(input$metaTable)) {
      metaTableDF <- hot_to_r(input$metaTable)
    } else {
      metaTableDF <- metaData$mdata
    }
    rhandsontable(metaTableDF, rowHeaderWidth = 350, overflow = "visible", useTypes = FALSE) %>%
      hot_context_menu(allowColEdit = TRUE)
  })
  # Putative effector metadata table
  output$emetaTable <- renderRHandsontable({
    req(emetaData$edata)
    if (!is.null(input$emetaTable)) {
      emetaTableDF <- hot_to_r(input$emetaTable)
    } else {
      emetaTableDF <- emetaData$edata
    }
    rhandsontable(emetaTableDF, rowHeaderWidth = 350, overflow = "visible", useTypes = FALSE) %>%
      hot_context_menu(allowColEdit = TRUE)
  })
  # Save genome metadata table updates by pressing the 'update plot' button
  updateData <- eventReactive({input$updatePlot},{
    metaData$mdata <- hot_to_r(input$metaTable)
  })
  # Save genome metadata table and update plot when 'update plot' button is pressed
  observeEvent(input$updatePlot, updateData())
  # Save p. effector metadata table updates by pressing the 'update plot' button
  eupdateData <- eventReactive({input$eupdatePlot},{
    emetaData$edata <- hot_to_r(input$emetaTable)
  })
  # Save p. effector metadata table and update plot when 'update plot' button is pressed
  observeEvent(input$eupdatePlot, eupdateData())
  # Download changed genome metadata table
  output$downloadTable <- downloadHandler(
    filename = function(){
      input$filename
    },
    content = function(file) {
      write.csv(metaData$mdata, file)
    }
  )
  # Download changed p. effector metadata table
  output$edownloadTable <- downloadHandler(
    filename = function(){
      input$efilename
    },
    content = function(file) {
      write.csv(emetaData$edata, file)
    }
  )
  # Get genome metadata table columns to provide them as options for annotation in the plot
  observe({
    req(metaData$mdata)
    updateSelectInput(session, "catSelect", choices = names(metaData$mdata))
  })
  # Get p. effector metadata table columns to provide them as options for annotation in the plot
  observe({
    req(metaData$mdata)
    updateSelectInput(session, "ecatSelect", choices = names(emetaData$edata))
  })
  # Create genome clusters
  genomeClusters <- reactive({
    data <- pav()
    dist_met_genomes <- input$dist_genomes
    clust_met_genomes <- input$clust_genomes
    if (dist_met_genomes == 'binary'){
      data[data > 0] <- 1
    }
    gdistance <- dist(data, method = dist_met_genomes, diag = FALSE, upper = FALSE)
    gcluster <- hclust(gdistance, method = clust_met_genomes)
    return(gcluster)
  })
  # Create effector clusters
  effectorClusters <- reactive({
    data <- pav()
    dist_met_effs <- input$dist_effs
    clust_met_effs <- input$clust_effs
    if (dist_met_effs == 'binary') {
      data[data > 0] <- 1
    }
    edistance <- dist(t(data), method = dist_met_effs, diag = FALSE, upper = FALSE)
    ecluster <- hclust(edistance, method = clust_met_effs)
    return(ecluster)
  })
  # Create heatmap plot
  pheat <- reactive({
    data <- pav()
    dist_met_genomes <- input$dist_genomes
    clust_met_genomes <- input$clust_genomes
    dist_met_effs <- input$dist_effs
    clust_met_effs <- input$clust_effs
    csize <- as.numeric(input$cfontsize)
    rsize <- as.numeric(input$rfontsize)
    pwidth <- as.numeric(input$heatWidth)/96
    pheight <- as.numeric(input$heatHeight)/96
    ncutree <- strtoi(input$ncuts)
    if (ncutree < 1) {
      ncutree <- 1
    } else if (ncutree > (nrow(data) - 1)) {
      ncutree <- (nrow(data) - 1)
    }
    bigScale = append(glasbey(), c("#967ACC", "#C7B8E6", "#CCCCCC"))
    bigScale = replace(bigScale, bigScale=="#201A01", "#A5A10E")
    if (dist_met_genomes == 'binary' || dist_met_effs == 'binary') {
      data[data > 0] <- 1
    }
    req(dist_met_genomes, dist_met_effs, clust_met_genomes, clust_met_effs)
    gcluster <- genomeClusters()
    ecluster <- effectorClusters()
    metaIn = metaData$mdata
    emetaIn = emetaData$edata
    # Genome color feature
    if (!is.null(metaIn)){
      # Find which column to use for annotation
      color_feature_name = input$catSelect
      # Get input from metadata table
      color_feature_index = grep(color_feature_name, colnames(metaIn))
      color_feature = metaIn[,color_feature_index, drop = FALSE]
      # Set empty values to NA, helps visualize
      color_feature[color_feature==""]<-NA
      # If all values are NA, no annotations have been added, so no visualization is needed
      if (!all(is.na(color_feature))){
        # Not all values are NA, but some may be
        na_present = FALSE
        # Get unique annotation values
        uniqueFeat = unique(color_feature[,1])
        # Number of unique annotations
        nCats = length(uniqueFeat)
        ######## To not waste unique colors, set '-' to white and not a unique color
        # If some of the values are NA
        if (any(is.na(color_feature))){
          # Everything except for NA
          nCats = nCats - 1
          uniqueFeat = uniqueFeat[!is.na(uniqueFeat)]
          # na_present = TRUE
        }
        ########
        # Get unique colors for annotation
        colors = unname(bigScale[1:nCats])
        names(colors) <- uniqueFeat
        # If all values are NA, set everything to NA to not visualize
      } else {
        color_feature = NA
        colors = NA
        nCats = NA
      }
      # If no file was uploaded, set everything to NA to not visualize
    } else {
      color_feature = NA
      colors = NA
      nCats = NA
    }
    # Putative effector color feature
    if (!is.null(emetaIn)){
      # Find which column to use for annotation
      ecolor_feature_name = input$ecatSelect
      # Get input from metadata table
      ecolor_feature_index = grep(ecolor_feature_name, colnames(emetaIn))
      ecolor_feature = emetaIn[,ecolor_feature_index, drop = FALSE]
      # Set empty values to NA, helps visualize
      ecolor_feature[ecolor_feature==""]<-NA
      # If all values are NA, no annotations have been added, so no visualization is needed
      if (!all(is.na(ecolor_feature))){
        # Not all values are NA, but some may be
        na_present = FALSE
        # Get unique annotation values
        euniqueFeat = unique(ecolor_feature[,1])
        # Number of unique annotations
        enCats = length(euniqueFeat)
        ######## To not waste unique colors, set '-' to white and not a unique color
        # If some of the values are NA
        if (any(is.na(ecolor_feature))){
          # Everything except for NA
          enCats = enCats - 1
          euniqueFeat = euniqueFeat[!is.na(euniqueFeat)]
          na_present = TRUE
        }
        ########
        # Get unique colors for annotation
        ecolors = unname(bigScale[1:enCats])
        names(ecolors) <- euniqueFeat
        # # Set '-' to white
        # if (na_present){
        #   ecolors = c(ecolors, c('-'='white'))
        # }
      # If all values are NA, set everything to NA to not visualize
      } else {
        ecolor_feature = NA
        ecolors = NA
        enCats = NA
      }
    # If no file was uploaded, set everything to NA to not visualize
    } else {
      ecolor_feature = NA
      ecolors = NA
      enCats = NA
    }
    if (dist_met_genomes == 'binary' || dist_met_effs == 'binary') {
      legBool = FALSE
      # TODO: Change NA to default number of breaks?
      nBreaks = NA
    } else {
      legBool = TRUE
      nBreaks = c(0, 0.9, 1.9, 4.9, 9.9, 20)
    }
    if (!any(is.na(colors)) && !any(is.na(ecolors))){
      final_names = c(color_feature_name, ecolor_feature_name)
      final_colors = list(a = colors, b = ecolors)
      names(final_colors) = final_names
    } else if (!any(is.na(colors))) {
      final_colors = list(colors)
      names(final_colors)[1] <- color_feature_name
    } else if (!any(is.na(ecolors))){
      final_colors = list(ecolors)
      names(final_colors)[1] <- ecolor_feature_name
    } else {
      final_colors = NA
    }
    myHeat <- pheatmap(data,
                       cluster_rows = gcluster,
                       cluster_cols = ecluster,
                       annotation_row = color_feature,
                       annotation_col = ecolor_feature,
                       annotation_colors = final_colors,
                       legend = legBool,
                       brewer.pal(n = 5, name = input$colors),
                       breaks = nBreaks,
                       show_colnames = input$ccheck,
                       show_rownames = input$rcheck,
                       cutree_rows = ncutree,
                       angle_col = 45,
                       fontsize_col = csize,
                       fontsize_row = rsize,
                       width = pwidth,
                       height = pheight)
    return(myHeat$gtable)
  })
  # Pheatmap plot
  output$heatPlot <- renderPlot({
    pheat()
  })
  # Download plot
  observe({
    input$heatname
    output$downloadHeat <- downloadHandler(
      filename = input$heatname,
      content = function(file) {
        pdf(file = file, width = ((as.numeric(input$heatWidth)/96) + 1), height = ((as.numeric(input$heatHeight)/96) + 1))
        grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
        grid.draw(pheat())
        dev.off()
      }
    )
  })
  # Reorder for CSV
  reorderData <- reactive({
    data = pav()
    if (input$dist_genomes == 'binary' || input$dist_effs == 'binary') {
      data[data > 0] <- 1
    }
    gcluster <- genomeClusters()
    ecluster <- effectorClusters()
    reordered = data[gcluster$order,]
    reordered = reordered[, ecluster$order]
    return(reordered)
  })
  # Download reordered csv
  observe({
    input$reorderedName
    output$downloadReorder <- downloadHandler(
      filename = input$reorderedName,
      content = function(file) {
        reordered = reorderData()
        write.csv(reordered, file)
      }
    )
  })
  # Control heatmap plot size
  resizedPlot <- reactive({
    pheight =  as.numeric(input$heatHeight)
    pwidth = as.numeric(input$heatWidth)
    resized <-  plotOutput("heatPlot", height = pheight, width = paste0(pwidth, 'px'))
    return(resized)
  })
  # Output resized plot
  output$heatPlotsized <- renderUI({
    resizedPlot()
  })
  # Genome dendrogram
  output$gCluster <- renderPlot({
    gcluster <- genomeClusters()
    par(cex = input$gtxtSize)
    plot(as.dendrogram(gcluster), horiz = input$ghoriz, main = paste0("GENOMES - Distance: ", input$dist_genomes, " Clustering: ", input$clust_genomes))
  })
  # Download genome dendrogram
  output$downloadGdend <- downloadHandler(
    filename = function(){
      'genomeNewick.txt'
    },
    content = function(file) {
      gcluster <- genomeClusters()
      writeLines(paste(write.dendrogram(as.dendrogram(gcluster), file = "")), file)
    }
  )
  # Effector dendrogram
  output$eCluster <- renderPlot({
    ecluster <- effectorClusters()
    par(cex = input$etxtSize)
    plot(as.dendrogram(ecluster), horiz = input$ehoriz, main = paste0("PUTATIVE EFFECTORS - Distance: ", input$dist_effs, " Clustering: ", input$clust_effs))
  })
  # Download effector dendrogram
  output$downloadEdend <- downloadHandler(
    filename = function(){
      'effectorNewick.txt'
    },
    content = function(file) {
      ecluster <- effectorClusters()
      writeLines(paste(write.dendrogram(as.dendrogram(ecluster), file = "")), file)
    }
  )
  output$msaPlot <- renderMsaR({
    inmsa <- file.path('../output/03.presenceabsence/putative_effector_msas/', input$msaFile)
    msaR(inmsa)
  })
}

shinyApp(ui = ui, server = server)