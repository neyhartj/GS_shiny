# The user interface script
library(shiny)

# Initialize the layout with a navbar
shinyUI(navbarPage("Simple Quantitative Genomics",
                   # First tab
                   tabPanel("Introduction",
                            titlePanel("Introduction"),
                            
                            mainPanel(
                              h2("Genomic Analysis of Complex Traits"),
                              br(),
                              h3("Why Study Complex Traits?"),
                              p("paragraph"),
                              br(),
                              h3("Associating Genomic Regions with Traits"),
                              br(),
                              h3("Genomic Prediction")
                            )
                   ),
                   
                   
                   # Second tab
                   tabPanel("Trait Analysis",
                            titlePanel("Trait Analysis"),
                            # Create a sidebar
                            sidebarLayout(
                              sidebarPanel(
                                
                                # First decide whether to use real or simulated data
                                helpText("Please choose whether you'd like to use simulated data, a historical barley dataset,
                                           or upload your own data."
                                ),
                                radioButtons("pheno_data_type1",
                                             label = "Choose a Data Type to Proceed",
                                             choices = c(
                                               "Real" = "real",
                                               "Simulated" = "simulated",
                                               "Historical" = "historical"
                                             ),
                                             selected = ""
                                ),
                                
                                # Conditional panel for historical data
                                conditionalPanel(
                                  condition = "input.pheno_data_type1 == 'historical'",
                                  
                                  # Help text to describe the data
                                  helpText("This historical data is malt extract from 768 breeding lines from
                                           the University of Minnesota and North Dakota State University barley
                                           breeding programs. The values have been centered at the mean."
                                  )
                                ),
                                
                                # Conditional panel for simulated data
                                conditionalPanel(
                                  condition = "input.pheno_data_type1 == 'simulated'",
                                  
                                  # Add some radio buttons for trait presets
                                  helpText("Here are some commonly measured traits in barley. Each has
                                           a preset as to the heritability and a constant genetic variance."),
                                  
                                  radioButtons("trait_presets",
                                               label = "Trait Presets or Custom",
                                               choices = c(
                                                 "Yield (h2 = 0.25)" = "yld",
                                                 "Height (h2 = 0.50)" = "ht",
                                                 "Diastatic Power (h2 = 0.75)" = "dp",
                                                 "Custom" = "custom"
                                               ),
                                               selected = "yld"
                                  )
                                ),
                                  
                                conditionalPanel(
                                  condition = "input.pheno_data_type1 == 'simulated' && input.trait_presets == 'custom'",

                                  # Genetic variance
                                  helpText("The genetic variance measures the variation at the genetic
                                           level. Use the slider to adjust genetic variance and observe its
                                           effect on the distribution."
                                  ),
                                  sliderInput("V_g",
                                              "Genetic Variance:",
                                              value = 5000,
                                              min = 0,
                                              max = 10000 
                                  ),
                                  
                                  # Heritability
                                  helpText("The heritability is the proportion of the observed variation
                                           that is due to variation at the genetic level. Use the slider
                                           to adjust the heritability."
                                           ),
                                  
                                  sliderInput("h2",
                                              "Heritability:",
                                              value = 0.5,
                                              min = 0.01,
                                              max = 1
                                  )
                                ),
                                  
                                # All other panels are conditional on using real data
                                conditionalPanel(
                                  condition = "input.pheno_data_type1 == 'real'",
                                  
                                  # Input options for trait data
                                  helpText("First upload a .csv of data on the plant traits (phenotypes)"
                                  ),
                                  fileInput(inputId = "pheno_file",
                                            label = "Upload Phenotype File",
                                            accept = c('.csv')
                                  )
                                )
                              ),
                              
                              # Main panel for the histogram
                              mainPanel(
                                tabsetPanel(type = "tabs",
                                            tabPanel("Plot", plotOutput("plot")),
                                            tabPanel("Summary", verbatimTextOutput("summary"))
                                )
                              )
                            )
                   ),
                   
                   # Third tab
                   tabPanel("Genetic Analysis",
                            titlePanel("Genetic Analysis"),
                              # Sidebar
                              sidebarLayout(
                                sidebarPanel(
                                  # Description
                                  helpText("Please choose whether you'd like to use a historical barley dataset
                                           or upload your own data."),
                                  
                                  radioButtons("geno_data_type1",
                                               label = "Choose a Data Type to Proceed",
                                               choices = c(
                                                 "Real" = "real",
                                                 "Historical" = "historical"
                                               ),
                                               selected = ""
                                  ),
                                  
                                  conditionalPanel(
                                    condition = "input.geno_data_type1 == 'real'",
                                    
                                    # Input options for trait data
                                    helpText("First upload a .csv of marker date (genotypes)."
                                    ),
                                    fileInput(inputId = "geno_file",
                                              label = "Upload Marker File",
                                              accept = c('.csv')
                                    )
                                  )
                                ),
                                  
                                # Main panel
                                mainPanel(
                                  tabsetPanel(type = "tabs",
                                              tabPanel("SNP Matrix", tableOutput("snp.mat")),
                                              tabPanel("Relationship Matrix", tableOutput("K.mat")),
                                              tabPanel("PCA Plot", plotOutput("K.PCA"))
                                              )
                                )
                              )
                   ),
                   
                   # Fourth tab
                   tabPanel("Marker-Trait Association",
                            titlePanel("Marker-Trait Association"),
                            # Side bar
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("gwas_data_type",
                                             label = "Choose a Data Type to Proceed",
                                             choices = c(
                                               "Real" = "real",
                                               "Simulated" = "simulated",
                                               "Historical" = "historical"
                                             ),
                                             selected = ""
                                ),
                                
                                # Conditional panel for simulated data
                                conditionalPanel(
                                  condition = "input.gwas_data_type == 'simulated'",
                                  
                                  # Add some radio buttons for trait presets
                                  helpText("Here are some commonly measured traits in barley. Each has
                                           a preset as to the heritability, number of QTL, and a constant genetic variance."),
                                  
                                  radioButtons("trait_presets2",
                                               label = "Trait Presets or Custom",
                                               choices = c(
                                                 "Yield (h2 = 0.25, #.QTL = 50)" = "yld",
                                                 "Height (h2 = 0.50, #.QTL = 25)" = "ht",
                                                 "Diastatic Power (h2 = 0.75, #.QTL = 10)" = "dp",
                                                 "Custom" = "custom"
                                               ),
                                               selected = "yld"
                                  )
                                  ),
                                
                                conditionalPanel(
                                  condition = "input.gwas_data_type == 'simulated' && input.trait_presets2 == 'custom'",
                                  
                                  # Number of QTL
                                  helpText("Detection of marker-trait associations or identifying putative QTL depends
                                           on a number of factors, including the number of QTL underlying a trait, the
                                           heritability of the trait, the number of markers assayed, and the number of 
                                           individuals observed. Use the following sliders to adjust these parameters:"
                                  ),
        
                                  sliderInput("gwas_n.qtl",
                                              "Number of QTL:",
                                              value = 5,
                                              min = 1,
                                              max = 50 
                                  ),
                                  
                                  sliderInput("gwas_h2",
                                              "Heritability:",
                                              value = 0.5,
                                              min = 0.01,
                                              max = 1
                                  ),
                                  
                                  sliderInput("gwas_n.markers",
                                              "Number of markers",
                                              value = 350,
                                              min = 7,
                                              max = 2000
                                  ),
                                  
                                  sliderInput("gwas_n.genos",
                                              "Population Size",
                                              value = 200,
                                              min = 10,
                                              max = 500
                                  ),
                                  
                                  helpText("Additionally, one can change the False Discovery Rate when analyzing the results
                                           of a GWAS. This parameter sets a threshold as to the number of significant detections,
                                           and also applies a Bonferroni correction to adjust for multiple testing."
                                  ),
                                  
                                  sliderInput("gwas_FDR",
                                              "False Discovery Rate",
                                              value = 0.05,
                                              min = 0.00001,
                                              max = 0.5
                                  )
                                  
                                ),
                                
                                
                                # Condition on historical data
                                conditionalPanel(
                                  condition = "input.gwas_data_type == 'historical'",
                                  
                                  # Help text
                                  helpText("The historical dataset uses 1074 SNP markers assayed on 768 breeding lines
                                             from the University of Minnesota and North Dakota State University barley breeding
                                             programs. The phenotype is malt extract. It may take some time to perform the GWAS.")
                                ),
                                
                                # All other panels are conditional on using real data
                                conditionalPanel(
                                  condition = "input.gwas_data_type == 'real'",
                                  
                                  # Input options for trait data
                                  helpText("Upload a .csv of data on the plant traits (phenotypes)"
                                  ),
                                  fileInput(inputId = "pheno_file2",
                                            label = "Upload Phenotype File",
                                            accept = c('.csv')
                                  ),
                                  
                                  # Input for marker data
                                  fileInput("geno_file2",
                                            label = "Upload Marker Genotype File",
                                            accept = c('.csv')
                                  )
                                )
                              ),
                              
                              # Main panel
                              mainPanel(
                                tabsetPanel(type = "tabs",
                                            tabPanel("Association Analysis", plotOutput("gwas"))
                                )
                              )
                            )
                            
                            
                   ),
                   
                   # Fifth tab
                   tabPanel("Genomic Prediction",
                            titlePanel("Genomic Prediction"),
                            sidebarLayout(
                              sidebarPanel(
                                # Radio buttons for the type of data
                                radioButtons("gs_data_type",
                                             label = "Choose a Data Type to Proceed",
                                             choices = c(
                                               "Real" = "real",
                                               "Simulated" = "simulated",
                                               "Historical" = "historical"
                                             ),
                                             selected = ""
                                ),
                                
                                # Conditional panel of cross-validation split
                                conditionalPanel(
                                  condition = "input.gs_data_type != ''",
                                  
                                  # Slider input
                                  sliderInput("cv_tp",
                                              "Cross-Validation Training Pop. Size:",
                                              value = 150,
                                              min = 10,
                                              max = 500
                                  )
                                )
                              ),
                              
                              # Main panel
                              mainPanel(
                                tabsetPanel(type = "tabs",
                                            tabPanel("Breeding Value Distribution", plotOutput("gebv")),
                                            tabPanel("Prediction Accuracy")
                                )
                              )
                            )
                   )
))
