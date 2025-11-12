

source('Functions.R')

options(shiny.maxRequestSize = 100 * 1024^2)


sidebar <- dashboardSidebar(disable = TRUE, width = 165,
                            useShinyjs(),
                            sidebarMenu(
                              menuItem("Upload File", tabName = "load", selected = TRUE)
                            )
)

body <- dashboardBody(
  useShinyjs(),
  
  # Make label text bigger
  tags$style(HTML("
    label.control-label {
      font-size: 20px !important;
    }
  ")),
  
  
  tabItems(
    
    tabItem(tabName = "load",
            fluidRow(
              
              column(width=2,
                  fluidRow(
                    fluidRow(align = "center", 
                             tags$style(type = "text/css", "
                              .radio label { 
                                font-size: 18px;
                              }
                            "),
                         
                             img(src = 'astral.png', align = "center", width = 180, height = 180 ),
                            
                             hr(),
                             
                             radioButtons( 
                               inputId = "astral_id", 
                               label = "Select Astral", 
                               choices = list( 
                                 "OA10034" = 1, 
                                 "OA10222" = 2
                               ),
                               selected = 1,
                             ), 
                             
                             radioButtons( 
                               inputId = "plot_type", 
                               label = "Plot Type?", 
                               choices = list( 
                                 "Line" = 1, 
                                 "Bar" = 2
                               ),
                               selected = 1,
                             ), 
                             
                             br(),
                             
                             tags$h4("Add comment..."),
                             tags$style(HTML("
                                .selectize-dropdown-content {
                                  font-size: 14px !important;
                                }
                                .selectize-input {
                                  font-size: 14px !important;
                                }
                              ")),
                             selectInput(inputId = "comment_date", label = NULL, width = '80%', choices = c("")),
                             textInput(inputId = "comment_text", label = NULL, placeholder = "Type new comment here..."),
                             actionButton(inputId = "add_comment", label = "Add Comment", 
                                          style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; display:center"),
                             hr(),
                             
                             tags$h4("Load New HeLa tsv file..."),
                             shinyFilesButton('sfb_hela_tsv_file', label = 'Load TSV File', title = 'Please select HeLa tsv file', multiple = FALSE,
                                              style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; display:center"),
                             hr(),
                             textOutput("hela_tsv_file_name"), style = "color:blue; font-size:9px"),
                  )
                ),
                column(width=10,  
                  tabBox(id = 'table', width = 12, height = 850,
                         tabPanel("Table",
                                  fluidRow(
                                    box(title = "Quartile Ratio Over Time", width = 12, status = "primary", solidHeader = TRUE,
                                        DTOutput("mytable", height = "600px"))
                                  )),
                         
                         tabPanel("Quartiles",
                            tabBox(id = 'quartile_plots', width = 12, height = 650,
                              tabPanel("First Quartile",     
                                  fluidRow(
                                    box(title = "Sum of First Quartile", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("q1Plot", height = "600px", brush = 'q1_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                           DTOutput('q1_brush_table', height = "600px"))
                                  )),
                              tabPanel("Second Quartile",     
                                       fluidRow(
                                         box(title = "Sum of Second Quartile", width = 12, status = "primary", solidHeader = TRUE,
                                             plotOutput("q2Plot", height = "600px", brush = 'q2_brush'))),
                                       fluidRow(
                                         column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                                DTOutput('q2_brush_table', height = "600px"))
                                       )),
                              tabPanel("Third Quartile",     
                                       fluidRow(
                                         box(title = "Sum of Third Quartile", width = 12, status = "primary", solidHeader = TRUE,
                                             plotOutput("q3Plot", height = "600px", brush = 'q3_brush'))),
                                       fluidRow(
                                         column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                                DTOutput('q3_brush_table', height = "600px"))
                                       )),
                              tabPanel("Fourth Quartile",     
                                       fluidRow(
                                         box(title = "Sum of Fourth Quartile", width = 12, status = "primary", solidHeader = TRUE,
                                             plotOutput("q4Plot", height = "600px", brush = 'q4_brush'))),
                                       fluidRow(
                                         column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                            DTOutput('q4_brush_table', height = "600px"))
                                       ))
                              )),
                         
                         tabPanel("Protein",
                                  fluidRow(box(title = "Total Proteins", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("proteinPlot", height = "600px", brush = 'protein_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput('protein_brush_table', height = "600px"))
                                  )),
                         
                         tabPanel("Peptide",
                                  fluidRow(box(title = "Total Peptides", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("peptidePlot", height = "600px", brush = 'peptide_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput('peptide_brush_table', height = "600px"))
                                  )),
                         
                         tabPanel("Precursor",
                                  fluidRow(box(title = "Total Precursors", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("precursorPlot", height = "600px", brush = 'precursor_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput('precursor_brush_table', height = "600px"))
                                  )),
                         
                         tabPanel("PeakWidth",
                            tabBox(id = 'peakwidth_plots', width = 12, height = 700,
                              tabPanel("Ideal PeakWidth",
                                  fluidRow(box(title = "Ratio of Ideal Peak Width (8-30sec)", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("ideal_peakwidthPlot", height = "600px", brush = 'ideal_peakwidth_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput('ideal_peakwidth_brush_table', height = "600px"))
                                  )),
                         
                              tabPanel("Wide PeakWidth Graph",
                                  fluidRow(box(title = "Ratio of Ideal Peak Width (>30sec)", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("wide_peakwidthPlot", height = "600px", brush = 'wide_peakwidth_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput('wide_peakwidth_brush_table', height = "600px"))
                                  )),
                         
                              tabPanel("Narrow PeakWidth Graph",
                                  fluidRow(box(title = "Ratio of Narrow Peak Width (<8sec)", width = 12, status = "primary", solidHeader = TRUE,
                                               plotOutput("narrow_peakwidthPlot", height = "600px", brush = 'narrow_peakwidth_brush'))),
                                  fluidRow(
                                    column(width = 12, tags$b(tags$i('Rows corresponding to datapoints selected')), 
                                      DTOutput("narrow_peakwidth_brush_table", height = "600px"))
                                  ))
                              
                            ))
                         
                  )
                )
          )
        )
      )
    )


dashboardPage(
  dashboardHeader(title = "Duke Proteomics Astral Quality Control", titleWidth = 400),
  sidebar,
  body
)



