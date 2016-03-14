# The server script
library(shiny)
library(rrBLUP, quietly = T)

# Load some CAP genotype data
load("CAP_data.RData")

M <- CAP.genos
K <- A.mat(M)
K.str <- svd(K)

# Define the manhattan plot function
# From Endelman et al 2011
manhattan <- function(input,fdr.level=0.05) {
  qvalue <- function(p) {
    smooth.df = 3
    
    if(min(p)<0 || max(p)>1) {
      print("ERROR: p-values not in valid range.")
      return(0)
    }
    
    lambda=seq(0,0.90,0.05)
    m <- length(p)
    
    pi0 <- rep(0,length(lambda))
    for(i in 1:length(lambda)) {pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])}
    
    spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
    pi0 <- predict(spi0,x=max(lambda))$y
    pi0 <- min(pi0,1)
    
    if(pi0 <= 0) {
      print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }
    
    #The estimated q-values calculated here
    u <- order(p)
    
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)
      
      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
      
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl
      
      return(tbl)
    }
    
    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)}
    
    return(qvalue)
  }
  
  #first column is marker name
  #second is chromosome
  #third is map position
  #fourth is score
  input <- input[order(input[,2],input[,3]),]
  
  chroms <- unique(input[,2])	
  n.chrom <- length(chroms)
  chrom.start <- rep(0,n.chrom)
  chrom.mid <- rep(0,n.chrom)
  
  if (n.chrom > 1) {
    for (i in 1:(n.chrom-1)) {chrom.start[i+1] <- chrom.start[i]+max(input[which(input[,2]==chroms[i]),3])+1}
  }
  x.max <- chrom.start[n.chrom]+max(input[which(input[,2]==chroms[n.chrom]),3])
  plot(0,0,type="n",xlim=c(0,x.max),ylim=c(0,max(input[,4])+1),ylab="-log(p)",xlab="Chromosome",xaxt="n")
  
  for (i in seq(1,n.chrom,by=2)) {
    ix <- which(input[,2]==chroms[i])
    chrom.mid[i] <- median(chrom.start[i]+input[ix,3])
    points(chrom.start[i]+input[ix,3],input[ix,4],col="dark blue",pch=16)
  }	
  
  if (n.chrom > 1){
    for (i in seq(2,n.chrom,by=2)) {
      ix <- which(input[,2]==chroms[i])
      chrom.mid[i] <- median(chrom.start[i]+input[ix,3])
      points(chrom.start[i]+input[ix,3],input[ix,4],col="cornflowerblue",pch=16)
    }	
  }
  
  q.ans <- qvalue(10^-input[,4])
  temp <- cbind(q.ans,input[,4])
  temp <- temp[order(temp[,1]),]	
  if (temp[1,1]<fdr.level) {
    temp2 <- tapply(temp[,2],temp[,1],mean)
    qvals <- as.numeric(rownames(temp2))
    x <- which.min(abs(qvals-fdr.level))
    first <- max(1,x-2)
    last <- min(x+2,length(qvals))
    if ((last-first)<4) {last <- first + 3}
    splin <- smooth.spline(x=qvals[first:last],y=temp2[first:last],df=3)
    lines(x=c(0,x.max),y=rep(predict(splin,x=fdr.level)$y,2),lty=2)
  }
  axis(side=1,at=chrom.mid,labels=chroms)
}


shinyServer( function(input, output) {
  
  traitData <- reactive({
    # Condition data based on datatype inpute
    ifelse(!any(c("real", "simulated", "historical") %in% input$pheno_data_type1),
      {return(NULL)},
      {
        if (input$pheno_data_type1 == "simulated") {
          # Set heritability based on presets or custom
          if (input$trait_presets == "yld") h2 = 0.25; V_g = 100
          if (input$trait_presets == "ht") h2 = 0.5; V_g = 100
          if (input$trait_presets == "dp") h2 = 0.75; V_g = 100
          if (input$trait_presets == "custom") h2 = input$h2; V_g = input$V_g
          n.genos = 200
          g <- rnorm(n = n.genos, mean = 0, sd = sqrt(V_g))
          Y <- g + rnorm(n.genos, mean = 0, sd = sqrt((1 - h2) / h2 * var(g)))
          
          return(Y)
        }
        
        if (input$pheno_data_type1 == "historical") {
          Y <- CAP.pheno[,2]
          Y <- Y - mean(Y, na.rm = T)
          return(Y)
        }
      
        if (input$pheno_data_type1 == "real") {
          # Only load the data if given a filepath
          if (is.null(input$pheno_file)) {
            return(NULL)
          }
          
          Y <- read.csv(input$pheno_file$datapath, row.names = 1, header = T)
          Y <- as.numeric(dataset[,1])
          return(Y)
        }
      }
    )
  })
  
  markerData <- reactive({
    # Condition on the selection of a radio button
    if (all(input$geno_data_type1 != "historical", input$geno_data_type1 != "real")) {
      return(NULL)
    } else {
      if (input$geno_data_type1 == "historical") {
        return(list(M = M, K = K, K.str = K.str, demog = demographics))
      }
      if (input$geno_data_type1 == "real") {
        if(is.null(input$geno_file)) {
          return(NULL)
        }
        
        M <- read.csv(input$geno_file$datapath, header = T, row.names = 1)
        K <- A.mat(M)
        K.str <- svd(K)
        
        return(list(M = M, K = K, K.str = K.str))
      }
    }
  })
  
  # Another marker data output
  markerData2 <- reactive({
    # Conditions
    ifelse(!any(c("real", "simulated", "historical") %in% input$gwas_data_type),
           return(NULL),
           {
             # Condition for real data
             if (input$gwas_data_type == "real") {
               if(is.null(input$geno_file2)) {
                 return(NULL)
               }
               
               M <- read.csv(input$geno_file2$datapath, header = T, row.names = 1)
               K <- A.mat(M)
               Y <- read.csv(input$pheno_file2$datapath, header = T, row.names = 1)
               
               return(list(geno = geno, K = K, pheno = pheno))
             }
             
             # Condition for historical data
             if (input$gwas_data_type == "historical") {
               M <- t(M)
               geno <- cbind(marker.data, data.frame(M, check.names = F))
               
               pheno <- data.frame(line = CAP.pheno[,1], Malt.extract = CAP.pheno[,2])
               pheno <- pheno[pheno$line %in% row.names(CAP.genos),]
               
               return(list(geno = geno, K = K, pheno = pheno))
             }
             
             # Condition on simulated data
             if (input$gwas_data_type == "simulated") {
               # Simulate marker data
               # Gather info
               m <- input$gwas_n.markers
               n <- input$gwas_n.genos
               chrom = sort(sample(1:7, size = m, replace = T))
               # Simulate
               M <- sapply(X = 1:n, FUN = function(x) ifelse(runif(m) < 0.5, -1, 1))
               colnames(M) <- paste("G", 1:n, sep = "")
               geno <- data.frame(
                 marker = paste("M", 1:m, sep = ""),
                 chrom = chrom,
                 pos = do.call("c", sapply(as.numeric(table(chrom)), function(n) 1:n)),
                 M
               )
               K <- A.mat(M)
               
               # Simulate QTL
               # Gather info
               ifelse(!any(c("yld", "ht", "dp") %in% input$trait_presets2),
                      {
                        n.QTL <- input$gwas_n.qtl
                        h2 <- input$gwas_h2
                      },
                      {
                        if (input$trait_presets2 == "yld") n.QTL = 50; h2 = 0.25
                        if (input$trait_presets2 == "ht") n.QTL = 25; h2 = 0.5
                        if (input$trait_presets2 == "dp") n.QTL = 10; h2 = 0.75
                      })
               
               QTL <- sample(m, n.QTL)
               u <- rep(0, m) # marker effects
               u[QTL] <- 1 # Assign QTL effects
               
               # Genotype values
               g <- as.vector(crossprod(M, u))
               Y <- g + rnorm(n, mean = 0, sd = sqrt((1-h2)/h2*var(g)))
               pheno <- data.frame(line = colnames(M), trait = Y)
               
               return(list(geno = geno, K = K, pheno = pheno))
             }
           })
  })
  
  # Marker data output for GS
  markerData3 <- reactive({
    # Conditions
    ifelse(!any(c("real", "simulated", "historical") %in% input$gs_data_type),
           {return(NULL)},
           {
             # Condition for real data
             if (input$gs_data_type == "real") {
               if(any(is.null(input$geno_file3), is.null(input$pheno_file3))) {
                 return(NULL)
               }
               
               M <- read.csv(input$geno_file3$datapath, header = T, row.names = 1)
               K <- A.mat(M)
               Y <- read.csv(input$pheno_file3$datapath, header = T, row.names = 1)
               
               return(list(geno = geno, K = K, pheno = pheno))
             }
             
             # Condition for historical data
             if (input$gs_data_type == "historical") {
               
               line.names <- row.names(M)
               n.pop <- length(line.names)
               
               # Size of the training pop
               n.TP <- input$cv_tp
               # Size of selection set
               n.VP <- n.pop - n.TP
               
               # Find the TP and VP lines
               TP.lines <- sample(line.names, size = n.TP)
               VP.lines <- setdiff(line.names, TP.lines)
               
               # Pull out data
               Y <- data.frame(row.names = CAP.pheno[,1], trait = CAP.pheno[,2])
               Y <- subset(Y, subset = row.names(Y) %in% row.names(M))
               
               # Create model parameters
               y.TP <- Y[TP.lines,]
               y.VP <- Y[VP.lines,]
               Z.TP <- M[TP.lines,]
               Z.VP <- M[VP.lines,]
               
               solve.out <- mixed.solve(y = y.TP, Z = Z.TP)
               mar.eff <- as.matrix(solve.out$u)
               # Calculate GEBVs
               GEBV <- Z.VP %*% mar.eff
               
               # Return a list
               return(list(mar.eff = mar.eff, GEBV = GEBV, y.VP = y.VP))
             }
             
             # Condition on simulated data
             if (input$gs_data_type == "simulated") {
               # Simulate marker data
               # Gather info
               m <- input$gs_n.markers
               n <- input$gs_n.genos
               
               # Simulate
               M <- t(sapply(X = 1:n, FUN = function(x) ifelse(runif(m) < 0.5, -1, 1)))
               row.names(M) <- paste("G", 1:n, sep = "")
               
               # Simulate QTL
               # Gather info
               ifelse(!any(c("yld", "ht", "dp") %in% input$trait_presets3),
                      {
                        n.QTL <- input$gs_n.qtl
                        h2 <- input$gs_h2
                      },
                      {
                        if (input$trait_presets3 == "yld") n.QTL = 50; h2 = 0.25
                        if (input$trait_presets3 == "ht") n.QTL = 25; h2 = 0.5
                        if (input$trait_presets3 == "dp") n.QTL = 10; h2 = 0.75
                      })
               
               QTL <- sample(m, n.QTL)
               u <- rep(0, m) # marker effects
               u[QTL] <- 1 # Assign QTL effects
               
               # Genotype values
               g <- as.vector(crossprod(t(M), u))
               Y <- g + rnorm(n, mean = 0, sd = sqrt((1-h2)/h2*var(g)))
               Y <- data.frame(row.names = row.names(M), trait = Y)
               
               # Start predictions
               line.names <- row.names(M)
               n.pop <- length(line.names)
               
               # Size of the training pop
               n.TP <- input$cv_tp
               # Size of selection set
               n.VP <- n.pop - n.TP
               
               # Find the TP and VP lines
               TP.lines <- sample(line.names, size = n.TP)
               VP.lines <- setdiff(line.names, TP.lines)
               
               # Create model parameters
               y.TP <- Y[TP.lines,]
               y.VP <- Y[VP.lines,]
               Z.TP <- M[TP.lines,]
               Z.VP <- M[VP.lines,]
               
               solve.out <- mixed.solve(y = y.TP, Z = Z.TP)
               mar.eff <- as.matrix(solve.out$u)
               # Calculate GEBVs
               GEBV <- Z.VP %*% mar.eff
               
               return(list(mar.eff = mar.eff, GEBV = GEBV, y.VP = y.VP))
             }
           })
    
  })


  # Plot the trait data
  output$plot <- renderPlot({
    if (!is.null(traitData())) {
    hist(traitData(),
         main = "Trait Distribution",
         xlab = "Phenotypic Value (Deviation from the Mean)",
         ylab = "Count",
         xlim = range(pretty(range(traitData(), na.rm = T))),
         breaks = 20)
    }
    })
      
  output$summary <- renderPrint({
    if (!is.null(traitData())) {
    summary(traitData())
    }
  })
  
  # Show the SNP matrix
  output$snp.mat <- renderTable({
    M <- markerData()$M
    M[c(1:10), c(1:10)]
  })
  
  # Show the relationship matrix
  output$K.mat <- renderTable({
    K <- markerData()$K
    K[c(1:10), c(1:10)]
  })
  
  # Plot population structure
  output$K.PCA <- renderPlot({
    # Plot PCA
    if (!is.null(markerData())) {
      K.str <- markerData()$K.str
      PC.x <- K.str$u[,1]
      PC.y <- K.str$u[,2]
      var.PC.x <- K.str$d[1] / sum(K.str$d)
      var.PC.y <- K.str$d[2] / sum(K.str$d)
      plot(K.str$u[,1], K.str$u[,2],
           xlab = paste("PC", "1", "(", round(var.PC.x, digits = 3)*100, "%)", sep = ""),
           ylab = paste("PC", "2", "(", round(var.PC.y, digits = 3)*100, "%)", sep = ""),
           main = "Genetic Relationship Principal Components Analysis",
           pch = 16,
           col = factor(markerData()$demog[,2])
      )
      legend(legend = levels(factor(markerData()$demog[,2])), "topleft", pch = 16, col = c("black", "red"))
    }
  })
  
  # Plot GWAS
  output$gwas <- renderPlot({
    if(!is.null(markerData2())) {
      pheno <- markerData2()$pheno
      geno <- markerData2()$geno
      K <- markerData2()$K
      tmp <- GWAS(pheno = pheno, geno = geno, K = K, plot = FALSE, P3D = T)
      manhattan(tmp, fdr.level = input$gwas_FDR)
    }
  })
  
  # Output GS results
  output$pred_plot <- renderPlot({
    
    if(!is.null(markerData3())) {
      GEBV <- markerData3()$GEBV
      y.VP <- markerData3()$y.VP
      plot(GEBV, y.VP, 
           main = "Genomic Prediction Accuracy",
           xlab = "Genomic Estimated Breeding Value",
           ylab = "Observed Phenotype")
      abline(lm(y.VP~GEBV))
    }
  })
  
  output$gebv <- renderPlot({
    
    if (!is.null(markerData3())) {
      GEBV <- markerData3()$GEBV
      hist(GEBV,
           main = "Distribution of Breeding Values"
      )
    }
  })
  
  output$pred_r <- renderText({
    
    if (!is.null(markerData3())) {
      GEBV <- markerData3()$GEBV
      y.VP <- markerData3()$y.VP
      cor(GEBV, y.VP, use = "complete.obs")
    }
  })
  
  
})













# model.functions <- list(
#   GBLUP = function(y, Z1, Z2) {
#     K <- A.mat(Z1)
#     solve.out <- mixed.solve(y = y, K = K)
#     GEBV <- return(solve.out$u)
#     return(GEBV)
#   },
#   RRBLUP = function(y, Z1, Z2) {
#     solve.out <- mixed.solve(y = y, Z = Z1)
#     GEBV <- Z2 %*% solve.out$u
#     return(GEBV)
#   }
# )
# 
# shinyServer(function(input, output) {
#   
#   selectedData <- reactive({
# 
#     # The inpute files are null originally
#     in.geno.train <- input$geno.train
#     in.geno.pred <- input$geno.pred
#     in.pheno.train <- input$pheno.train
# 
#     # Store the model choice
#     model <- input$model
# 
#     if (any(is.null(in.geno.train), is.null(in.pheno.train))) {
#     # if (any(is.null(in.geno.train), is.null(in.geno.pred), is.null(in.pheno.train))) {
#       return("No input to display")
#     } else {
# 
#       geno.train <- read.csv(in.geno.train$datapath, header = T, row.names = 1)
#       # read.csv(in.geno.pred$datapath)
#       pheno.train <- read.csv(in.pheno.train$datapath, header = T, row.names = 1)
# 
#       # Separate by trait
#       y <- pheno.train[,1]
#       Z1 <- geno.train
# 
#       GEBV <- model.functions[[model]](y = y, Z1 = Z1)
#       return(GEBV)
#     }
#   })
#   
#   # Output
#   output$summary <- renderPrint({
#     summary(selectedData)
#   })
#   
#   # Plot
#   # if (!is.null(selectedData)) {
#   #   output$plot <- renderPlot({
#   #     hist(selectedData)
#   #   })
#   # }
#   
# })
