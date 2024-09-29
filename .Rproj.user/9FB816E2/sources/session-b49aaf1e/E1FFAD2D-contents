# Load required packages
library(shiny)
library(Seurat)
library(bslib)
library(shinyjs)
library(DT)

# Load data
wt_integrate <- readRDS("mini_trajatlas.Rds")
wt_integrate@assays$originalexp@key <- "rna_"
descript <- read.csv("S2_Table_cell_description.csv")
anno <- readxl::read_xlsx("S3_Table_ssc_anno.xlsx")
objList <- readRDS("20240929_stemcell_highlight.Rds")
none_value <- rep("None", dim(descript)[2])
descript <- rbind(descript, none_value)
geneDatabase <- read.csv("S6_Table_OPCST_gene_combine_database.csv", skip = 20)


# Define the UI
ui <- page_navbar(
  title = "Overview of Differentiation Atlas",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel(
    "Annotation",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("name", "Select a slot", 
                    choices = c("level-1-anno", "level-2-anno", "level-3-anno", "level-4-anno", "level-5-anno", 
                                "level-6-anno", "level-7-anno", "Project", "Tissue origin", "Tissue location", "Tissue.Specific.", 
                                "Stage", "Gene.type", "Treatment", "Age", "Age.In.Detail.", "Machine", 
                                "Origin", "paper_label", "coarse_label")),
        uiOutput("ui"),
        actionButton("submit", "Submit", class = "btn-primary")
      ),
      layout_column_wrap(
        width = 1/2,
        card(
          card_header("Cell Information"),
          DTOutput("cellInfoTable")
        ),
        card(
          card_header("Overview"),
          imageOutput("image", height = "auto")
        )
      ),
      card(
        card_header("Visualization"),
        plotOutput("seuratPlot", width = "100%", height = "600px")
      )
    )
  ),
  nav_panel(
    "OPC Mapping",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("stemName", "Select stem cells' name", 
                    choices = c("Pdgfra+ BMSC", "Acta2+ BMSC", "Lepr+ BMSC", "Cxcl12+ BMSC", 
                                "Clec11a+ BMSC", "Nes+ BMSC", "Grem1+ BMSC", "pvSSC", "Cdh2+ BMSC", 
                                "Mcam+ BMSC", "Hypertrophic chondrocyte", "Chondrocyte", "Early perichondrial cell", 
                                "Hoxa11+ Mesenchyme", "ocSSC", "Pthlh+ chondrocyte", "Foxa2+ chondrocyte", 
                                "Cd168+ skeletal stem/progenitor cell", "Axin2+ Mes", "Hhip+ Mes", 
                                "Msx2+ Mes", "Gli1+ Mes", "Prrx1+ Mes", "Nfatc1+ Chondro", "Fgfr3+ Chondro", 
                                "Hes1+ Perichondrium", "Sstr2+ MSC", "Sox9+ chondrocyte")),
        actionButton("submitStem", "Submit", class = "btn-primary")
      ),
      layout_column_wrap(
        width = 1/2,
        card(
          card_header("Stem Cell Information"),
          DTOutput("stemCellInfoTable")
        ),
        card(
          card_header("OPC Mapping Overview"),
          imageOutput("imageOPC", height = "auto")
        )
      ),
      card(
        card_header("Visualization"),
        plotOutput("seuratPlotStem", width = "100%", height = "600px")
      )
    )
  ),
  nav_panel(
    "Gene Explorer",

      card(
        card_header("Gene Database"),
        DTOutput("geneTable")
      )
    ),
  nav_panel(
    "Citation",
    card(
      card_header("Please cite our work"),
      verbatimTextOutput("citation")
    )
  )
)

# Define the server
server <- function(input, output) {
  output$ui <- renderUI({
    if (is.null(input$name))
      return()
    selectInput("dynamic", "Select an idents",
                choices = unique(wt_integrate@meta.data[[input$name]]),
                selected = "option2"
    )
  })
  
  selected_row <- reactive({
    if (!is.null(input$dynamic)) {
      if (input$dynamic %in% descript$X) {
        descript[descript$X == input$dynamic, ]
      } else {
        descript[descript$X == "None", ]
      }
    }
  })
  
  selected_ssc <- reactive({
    anno[anno$`Stem cell name` == input$stemName,]
  })
  
  output$image <- renderImage({
    list(src = "www/mtree.png",
         contentType = "image/png",
         width = "100%")
  }, deleteFile = FALSE)
  
  output$imageOPC <- renderImage({
    list(src = "www/OPCMapping.png",
         contentType = "image/png",
         width = "100%")
  }, deleteFile = FALSE)
  
  observeEvent(input$submit, {
    row <- selected_row()
    
    output$cellInfoTable <- renderDT({
      data.frame(
        Category = c("Description", "Annotation level", "Predominant tissue origin", 
                     "Predominant tissue location", "Predominant tissue (histology)", 
                     "Predominant age", "Predominant age (detail)"),
        Value = c(row[["Description"]], row[["Anno_level"]], row[["Organ"]], 
                  row[["Tissue"]], row[["Tissue.Specific."]], row[["Age"]], 
                  row[["Age.In.Detail."]])
      )
    }, options = list(dom = 't', paging = FALSE, ordering = FALSE, searching = FALSE))
    
    highlighted_cell <- colnames(wt_integrate)[wt_integrate@meta.data[[input$name]] == input$dynamic]
    
    output$seuratPlot <- renderPlot({
      p1 <- DimPlot(wt_integrate, cells.highlight = highlighted_cell, raster = TRUE)
      if (input$name %in% c("level-5-anno", "level-6-anno", "level-7-anno")) {
        p2 <- DimPlot(wt_integrate, group.by = input$name, raster = TRUE) + NoLegend() 
      } else {
        p2 <- DimPlot(wt_integrate, group.by = input$name, raster = TRUE)
      }
      p1 + p2
    })
  })
  
  observeEvent(input$submitStem, {
    row2 <- selected_ssc()
    highlighted_cell <- objList[[input$stemName]]
    
    output$stemCellInfoTable <- renderDT({
      data.frame(
        Category = c("Stem cell name", "Marker", "Publication", "Cluster in Atlas", 
                     "Tissue location","Validation Methods", "Driver", "Confidence level"
        ),
        Value = c(
          row2[["Stem cell name"]],
          row2[["Marker"]],
          paste0("Published in ", row2[["Journal"]], ". Article: '", row2[["Publication"]], "', Year: ", row2[["Year"]]),
          row2[["Cluster_DifferentiationAtlas"]],
          paste0("Tissue source: ", row2[["Tissue source"]], ", Tissue location: ", row2[["Tissue location"]]),
          row2[["Validate Methods"]],
          row2[["Driver"]],
          paste0(row2[["Confidence level"]], " confidence that this stem cell matches the area in our atlas")
        )
      )
    }, options = list(dom = 't', paging = FALSE, ordering = FALSE, searching = FALSE))
    
    output$seuratPlotStem <- renderPlot({
      p1 <- DimPlot(wt_integrate, cells.highlight = highlighted_cell, raster = TRUE)
      p2 <- DimPlot(wt_integrate, group.by = "Tissue origin", raster = TRUE)
      p3 <- DimPlot(wt_integrate, group.by = "level-2-anno", raster = TRUE)
      p1 + p2 + p3
    })
  })
  
  output$geneTable <- renderDT({
    datatable(geneDatabase, 
              options = list(pageLength = 10, 
                             lengthMenu = c(10, 25, 50, 100),
                             scrollX = TRUE),
              filter = 'top',
              selection = 'single')
  })
  observeEvent(input$geneTable_rows_selected, {
    selected_row <- input$geneTable_rows_selected
    if (length(selected_row) > 0) {
      gene <- geneDatabase[selected_row, input$geneColumn]

    }
  })
  
  output$citation <- renderText({
    '@article {Han2024.05.28.596174,
    author = {Han, Litian and Ji, Yaoting and Yu, Yiqian and Ni, Yueqi and Zeng, Hao and Zhang, Xiaoxin and Liu, Huan and Zhang, Yufeng},
    title = {Trajectory-centric Framework TrajAtlas reveals multi-scale differentiation heterogeneity among cells, genes, and gene module in osteogenesis},
    elocation-id = {2024.05.28.596174},
    year = {2024},
    doi = {10.1101/2024.05.28.596174},
    publisher = {Cold Spring Harbor Laboratory},
    abstract = {Osteoblast differentiation is crucial for bone formation and maintaining skeletal integrity. Although it is now understood that this process exhibits significant heterogeneity across developmental stages and tissue microenvironments, the underlying mechanisms remain largely unexplored. In the present study, we introduce TrajAtlas, a comprehensive framework that addresses this gap in knowledge. TrajAtlas comprises four modules: a reference atlas (Differentiation Atlas), a differentiation model (Differentiation Model), a tool for differential pseudotime analysis (TrajDiff), and a method for pseudotemporal gene module detection (TRAVMap). By leveraging single-cell technologies, TrajAtlas offers a systematic approach to exploring the multi-scale heterogeneity among cells, genes, and gene modules within population-level trajectories across diverse tissues and age groups. We systematically investigate the impact of age and injury on osteogenesis, providing new insights into osteoporosis and bone regeneration. In conclusion, our comprehensive framework offers novel insights into osteogenesis and provides a valuable resource for understanding the complexities of bone formation.Author Summary Osteoblasts, the cells responsible for bone formation, can originate from various cellular sources. However, it\'s unclear how different progenitor cells differentiate into osteoblasts, and how this process is influenced by factors such as age and tissue location. This knowledge gap stems from the lack of comprehensive databases and tools to decipher the differentiation process. In this study, we introduce TrajAtlas, a comprehensive framework designed to bridge this gap. To explore the cellular origins of osteoblasts, we constructed an atlas centered on osteogenesis. To answer how progenitor cells differentiate to osteoblasts, we developed a model that reveals the dynamic regulatory landscape during this process. To elucidate the influence of age and tissue location on differentiation, we built a tool for differential analysis. Furthermore, to identify conserved patterns of differentiation, we developed an approach to detect pseudotemporal gene modules. We validated the effectiveness of this framework by applying it to more datasets, unveiling novel cell states associated with injury. Notably, this framework focuses on dynamic processes, with the potential for broader applications in studying cell differentiation and complementing cell-centric analyses.Competing Interest StatementThe authors have declared no competing interest.},
    URL = {https://www.biorxiv.org/content/early/2024/06/02/2024.05.28.596174},
    eprint = {https://www.biorxiv.org/content/early/2024/06/02/2024.05.28.596174.full.pdf},
    journal = {bioRxiv}
}'
  })
}

# Run the app
shinyApp(ui, server)