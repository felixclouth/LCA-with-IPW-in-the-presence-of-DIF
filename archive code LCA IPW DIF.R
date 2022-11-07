###############################################################################################
# last edit: 04/11/2022

# accompanying code to the manuscript:
# "Three-step latent class analysis with inverse propensity weighting in the presence of differential item functioning"

# F.J. Clouth, MSc.
# Department of Methodology and Statistics, Tilburg University
# f.j.clouth@tilburguniversity.edu

library(foreign)
library(rio)
library(dplyr)
library(ggplot2)
library(ggsci)
library(survey)
library(tableone)
library(car)
library(mice)


setwd("adjust")
LG <- "adjust/LatentGOLD6.0/lg60.exe"
# latent Gold can be downloaded at https://www.statisticalinnovations.com/latent-gold-6-0/
# latent Gold syntax is generate and run from this code

###############################################################################################


# generate latent Glod syntax to simluate data that is perfectly in accordance with the population values (writeexemplary-option)
GenData <- function(syntaxName, infile, outfile, B, G, DE, D){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    options
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=0 nriterations=0;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=0;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 replicates=500 tolerance=1e-008;
                                                    missing  includeall;
                                                    output      
                                                    estimatedvalues parameters=first estimatedvalues profile reorderclasses iterationdetails
                                                    writeexemplary='", outfile,"';
                                                    variables
                                                    dependent Z nominal 2, y1 nominal 2, y2 nominal 2, y3 nominal 2, y4 nominal 2, y5 nominal 2, y6 nominal 2;
                                                    independent C1, C2;
                                                    latent Zlat nominal 2, Cluster nominal 3;
                                                    equations
                                                    Zlat <- 1 + C1 + C2;
                                                    Cluster <- 1 + C1 + C2 + Zlat;
                                                    Z <- (w~wei) 1 | Zlat; 
                                                    Y1 <- 1 | Cluster;
                                                    Y2 <- 1 | Cluster;
                                                    Y3 <- 1 | Cluster + ", DE,";
                                                    Y4 <- 1 | Cluster;
                                                    Y5 <- 1 | Cluster;
                                                    Y6 <- 1 | Cluster + ", DE,";
                                                    
                                                    w = {1 0 0 1};
                                                    
                                                    {0 1 ", B," .5 .5 1 1 1 1 1 ", G,"
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361 ", D,"
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361 ", D,"
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# varying parameters in our simulated data example: B=beta, G=gamma, DE=treatment or confounders have direct effects, D=strength of the direct effect
# note that not all scenarios are reported in the manuscript
B <- c(1, 2, 3)
G <- c(1, 2, 3)
D <- c(1, 2)
DE <- c("C1", "Z")


# run GenData function to generate the syntax files and run these syntax files in latent Gold
for(b in B) {
  for(g in G) {
    for(d in D) {
      for(de in DE) {
        des <- ifelse(de == "C1", "C1", "Z")
        GenData(syntaxName = paste0("GenData_B",b, "G",g, "D",d, "DE",de), infile = "example.sav", outfile = paste0("simB",b, "G",g, "D",d, "DE",des, ".dat"), B = b, G = g, D = d, DE = de)
        shell(paste0(LG, " ", "GenData_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        print(c(b, g, d, de))
      }
    }
  }
}


# estimate IPW and attach to simulated data sets
for(b in B) {
  for(g in G) {
    for(d in D) {
      for(de in DE) {
        des <- ifelse(de == "C1", "C1", "Z")
        sim <- read.table(paste0("simB",b, "G",g, "D",d, "DE",des, ".dat"), header = T)
        sim$Z[sim$Z==2] <- 0
        
        MyGLM <- glm(Z ~ C1 + C2, family = "binomial", data = sim, weights = frequency)
        PS_est <- MyGLM$fitted.values
        
        sim$PS <- PS_est
        sim$IPW <- ifelse(sim$Z == 1, 1/sim$PS, 1/(1-sim$PS))
        
        export(sim, paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"))
        
        print(c(b, g, d, de))
      }
    }
  }
}


# generate syntax files for the different estimation methods

# one-step regression adjustment without accounting for DIF
regadjmis <- function(syntaxName, infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title '1-step misspecified';
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first marginaleffects standarderrors profile bivariateresiduals estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    independent Z nominal, C1, C2;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + Z + C1 + C2;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0 0 0 0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# one-step regression adjustment with accounting for DIF
regadjcor <- function(syntaxName, infile, outfile, DE, D){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title '1-step correct';
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first marginaleffects standarderrors profile bivariateresiduals estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    independent Z nominal, C1, C2;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + Z + C1 + C2;
                                                    Y1 <- 1 | Cluster;
                                                    Y2 <- 1 | Cluster;
                                                    Y3 <- 1 | Cluster + ", DE,";
                                                    Y4 <- 1 | Cluster;
                                                    Y5 <- 1 | Cluster;
                                                    Y6 <- 1 | Cluster + ", DE,";
                                                    {0 0 0 0 0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361 ", D,"
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361 ", D,"
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# one-step IPW without accounting for DIF (Lanza et al., 2013)
Lanzamis <- function(syntaxName, infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'Lanza misspecified';
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first marginaleffects standarderrors profile bivariateresiduals estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    samplingweight IPW rescale;
                                                    independent Z nominal;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# one-step IPW with accounting for DIF (adjusted Lanza)
Lanzacor <- function(syntaxName, infile, outfile, DE, D, formula, form, par){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'Lanza correct';
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first marginaleffects standarderrors profile bivariateresiduals estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    samplingweight IPW rescale;
                                                    independent ", form,";
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + ", formula,";
                                                    Y1 <- 1 | Cluster;
                                                    Y2 <- 1 | Cluster;
                                                    Y3 <- 1 | Cluster + ", DE,";
                                                    Y4 <- 1 | Cluster;
                                                    Y5 <- 1 | Cluster;
                                                    Y6 <- 1 | Cluster + ", DE,";
                                                    {", par," 
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361 ", D,"
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361 ", D,"
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# step one of three-step IPW without accounting for DIF (Clouth et al., 2022)
step1mis <- function(syntaxName, infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'misspecified step 1';
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile bivariateresiduals estimatedvalues=model;
                                                    outfile  '", outfile,"' classification=posterior keep C1, C2, Z, IPW;
                                                    variables
                                                    caseweight frequency;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# step one of three-step IPW with accounting for DIF
step1cor <- function(syntaxName, infile, outfile, DE, D, form2){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'correct step 1';
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile bivariateresiduals estimatedvalues=model;
                                                    outfile  '", outfile,"' classification=posterior keep ", form2,", IPW;
                                                    variables
                                                    caseweight frequency;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    independent ", DE,";
                                                    latent
                                                    ClusterDIF nominal 3;
                                                    equations
                                                    ClusterDIF <- 1 + ", DE,";
                                                    Y1 <- 1 | ClusterDIF;
                                                    Y2 <- 1 | ClusterDIF;
                                                    Y3 <- 1 | ClusterDIF + ", DE,";
                                                    Y4 <- 1 | ClusterDIF;
                                                    Y5 <- 1 | ClusterDIF;
                                                    Y6 <- 1 | ClusterDIF + ", DE,";
                                                    {0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361 ", D,"
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361 ", D,"
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# step three of three-step IPW without accounting for DIF (Clouth et al., 2022)
step3mis1 <- function(syntaxName, infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'Step3 with misspecified step1';
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    step3 proportional ml; 
                                                    output      
                                                    parameters=first marginaleffects standarderrors=robust profile estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    samplingweight IPW rescale ipw;
                                                    independent Z nominal;
                                                    latent Cluster nominal posterior = ( Cluster#1 Cluster#2 Cluster#3 ) ;
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# step three of three-step IPW with accounting for DIF but with misspecified step one (not reported in the manuscript)
step3cor1 <- function(syntaxName, infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'Step3 with correct step1';
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    step3 proportional ml;
                                                    output      
                                                    parameters=first marginaleffects standarderrors=robust profile estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    samplingweight IPW rescale ipw;
                                                    independent Z nominal;
                                                    latent Cluster nominal posterior = ( ClusterDIF#1 ClusterDIF#2 ClusterDIF#3 ) ;
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# step three of three-step IPW with accounting for DIF
step3DIF <- function(syntaxName, infile, outfile, form, DE){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    model
                                                    title 'Step3 with correct step1 and DIF variable';
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=0 variances=0 latent=0 poisson=0;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    step3 proportional ml;
                                                    output      
                                                    parameters=first marginaleffects standarderrors=robust profile estimatedvalues=model append='", outfile,"';
                                                    variables
                                                    caseweight frequency;
                                                    samplingweight IPW rescale ipw;
                                                    independent ", form,";
                                                    latent Cluster nominal posterior = ( ClusterDIF#1 ClusterDIF#2 ClusterDIF#3 ) dif=", DE,";
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# generate syntax files and run in latent Gold
for(b in B) {
  for(g in G) {
    for(d in D) {
      for(de in DE) {
        des <- ifelse(de == "C1", "C1", "Z")
        form <- ifelse(de == "C1", "Z nominal, C1", "Z nominal")
        form2 <- ifelse(de == "C1", "Z, C2", "C1, C2")
        formula <- ifelse(de == "C1", "Z + C1", "Z")
        par <- ifelse(de == "C1", "0 0 0 0 0 0", "0 0 0 0")
        
        regadjmis(syntaxName = paste0("regadjmis_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("regadjmis_B",b, "G",g, "D",d, "DE",de, ".csv"))
        regadjcor(syntaxName = paste0("regadjcor_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("regadjcor_B",b, "G",g, "D",d, "DE",de, ".csv"), D = d, DE = de)
        
        Lanzamis(syntaxName = paste0("Lanzamis_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("Lanzamis_B",b, "G",g, "D",d, "DE",de, ".csv"))
        Lanzacor(syntaxName = paste0("Lanzacor_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("Lanzacor_B",b, "G",g, "D",d, "DE",de, ".csv"), D = d, DE = de, formula = formula, form = form, par = par)
        
        step1mis(syntaxName = paste0("step1mis_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("step1mis_B",b, "G",g, "D",d, "DE",de, ".sav"))
        step1cor(syntaxName = paste0("step1cor_B",b, "G",g, "D",d, "DE",de), infile = paste0("simB",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("step1cor_B",b, "G",g, "D",d, "DE",de, ".sav"), D = d, DE = de, form2 = form2)
        
        step3mis1(syntaxName = paste0("step3mis1_B",b, "G",g, "D",d, "DE",de), infile = paste0("step1mis_B",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("step3mis1_B",b, "G",g, "D",d, "DE",de, ".csv"))
        step3cor1(syntaxName = paste0("step3cor1_B",b, "G",g, "D",d, "DE",de), infile = paste0("step1cor_B",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("step3cor1_B",b, "G",g, "D",d, "DE",de, ".csv"))
        step3DIF(syntaxName = paste0("step3DIF_B",b, "G",g, "D",d, "DE",de), infile = paste0("step1cor_B",b, "G",g, "D",d, "DE",des, ".sav"), outfile = paste0("step3DIF_B",b, "G",g, "D",d, "DE",de, ".csv"), DE = de, form = form)
        
        
        shell(paste0(LG, " ", "regadjmis_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        shell(paste0(LG, " ", "regadjcor_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        
        shell(paste0(LG, " ", "Lanzamis_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        shell(paste0(LG, " ", "Lanzacor_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        
        shell(paste0(LG, " ", "step1mis_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        shell(paste0(LG, " ", "step1cor_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        
        shell(paste0(LG, " ", "step3mis1_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        shell(paste0(LG, " ", "step3cor1_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        shell(paste0(LG, " ", "step3DIF_B",b, "G",g, "D",d, "DE",de, ".lgs", " ", "/b"))
        
        print(c(b, g, d, de))
      }
    }
  }
}


# calculations for the final results presented in the manuscript
TrueATE <- c(0.0664870, 0.3019906, 0.4824338)

regadjmisATE <- array(NA, dim = c(3, 3, 2, 2))
regadjcorATE <- array(NA, dim = c(3, 3, 2, 2))

LanzamisATE <- array(NA, dim = c(3, 3, 2, 2))
LanzacorATE <- array(NA, dim = c(3, 3, 2, 2))

step3mis1ATE <- array(NA, dim = c(3, 3, 2, 2))
step3cor1ATE <- array(NA, dim = c(3, 3, 2, 2))
step3DIFATE <- array(NA, dim = c(3, 3, 2, 2))

regadjmisBIAS <- array(NA, dim = c(3, 3, 2, 2))
regadjcorBIAS <- array(NA, dim = c(3, 3, 2, 2))

LanzamisBIAS <- array(NA, dim = c(3, 3, 2, 2))
LanzacorBIAS <- array(NA, dim = c(3, 3, 2, 2))

step3mis1BIAS <- array(NA, dim = c(3, 3, 2, 2))
step3cor1BIAS <- array(NA, dim = c(3, 3, 2, 2))
step3DIFBIAS <- array(NA, dim = c(3, 3, 2, 2))


for(b in B) {
  for(g in G) {
    for(d in D) {
      for(de in DE) {
        des <- ifelse(de == "C1", 1, 2)
        
        regadjmisResults <- read.table(paste0("regadjmis_B",b, "G",g, "D",d, "DE",de, ".csv"))
        regadjmisMarginal <- as.character(regadjmisResults$V2[7])
        regadjmisMarginal <- strsplit(regadjmisMarginal, ",")
        regadjmisATE[b, g, d, des] <- matrix(as.numeric(regadjmisMarginal[[1]][11]), 1, 1, byrow = T)
        
        regadjcorResults <- read.table(paste0("regadjcor_B",b, "G",g, "D",d, "DE",de, ".csv"))
        regadjcorMarginal <- as.character(regadjcorResults$V2[7])
        regadjcorMarginal <- strsplit(regadjcorMarginal, ",")
        regadjcorATE[b, g, d, des] <- matrix(as.numeric(regadjcorMarginal[[1]][11]), 1, 1, byrow = T)
        
        
        LanzamisResults <- read.table(paste0("Lanzamis_B",b, "G",g, "D",d, "DE",de, ".csv"))
        LanzamisMarginal <- as.character(LanzamisResults$V2[7])
        LanzamisMarginal <- strsplit(LanzamisMarginal, ",")
        LanzamisATE[b, g, d, des] <- matrix(as.numeric(LanzamisMarginal[[1]][11]), 1, 1, byrow = T)
        
        LanzacorResults <- read.table(paste0("Lanzacor_B",b, "G",g, "D",d, "DE",de, ".csv"))
        LanzacorMarginal <- as.character(LanzacorResults$V2[7])
        LanzacorMarginal <- strsplit(LanzacorMarginal, ",")
        LanzacorATE[b, g, d, des] <- matrix(as.numeric(LanzacorMarginal[[1]][11]), 1, 1, byrow = T)
        
        
        step3mis1Results <- read.table(paste0("step3mis1_B",b, "G",g, "D",d, "DE",de, ".csv"))
        step3mis1Marginal <- as.character(step3mis1Results$V4[6])
        step3mis1Marginal <- strsplit(step3mis1Marginal, ",")
        step3mis1ATE[b, g, d, des] <- matrix(as.numeric(step3mis1Marginal[[1]][11]), 1, 1, byrow = T)
        
        step3cor1Results <- read.table(paste0("step3cor1_B",b, "G",g, "D",d, "DE",de, ".csv"))
        step3cor1Marginal <- as.character(step3cor1Results$V4[6])
        step3cor1Marginal <- strsplit(step3cor1Marginal, ",")
        step3cor1ATE[b, g, d, des] <- matrix(as.numeric(step3cor1Marginal[[1]][11]), 1, 1, byrow = T)
        
        step3DIFResults <- read.table(paste0("step3DIF_B",b, "G",g, "D",d, "DE",de, ".csv"))
        step3DIFMarginal <- as.character(step3DIFResults$V7[6])
        step3DIFMarginal <- strsplit(step3DIFMarginal, ",")
        step3DIFATE[b, g, d, des] <- matrix(as.numeric(step3DIFMarginal[[1]][11]), 1, 1, byrow = T)
        
        regadjmisBIAS[b, g, d, des] <- TrueATE[g] - regadjmisATE[b, g, d, des]
        regadjcorBIAS[b, g, d, des] <- TrueATE[g] - regadjcorATE[b, g, d, des]
        
        LanzamisBIAS[b, g, d, des] <- TrueATE[g] - LanzamisATE[b, g, d, des]
        LanzacorBIAS[b, g, d, des] <- TrueATE[g] - LanzacorATE[b, g, d, des]
        
        step3mis1BIAS[b, g, d, des] <- TrueATE[g] - step3mis1ATE[b, g, d, des]
        step3cor1BIAS[b, g, d, des] <- TrueATE[g] - step3cor1ATE[b, g, d, des]
        step3DIFBIAS[b, g, d, des] <- TrueATE[g] - step3DIFATE[b, g, d, des]
      }
    }
  }
}


regadjcorATE[1, , 1, 1]
regadjcorATE[1, , 1, 2]
regadjmisATE[1, , 1, 1]
regadjmisATE[1, , 1, 2]
LanzacorATE[1, , 1, 1]
LanzacorATE[1, , 1, 2]
LanzamisATE[1, , 1, 1]
LanzamisATE[1, , 1, 2]
step3DIFATE[1, , 1, 1]
step3DIFATE[1, , 1, 2]
step3mis1ATE[1, , 1, 1]
step3mis1ATE[1, , 1, 2]


method <- rep(c("regadjcor", "regadjmis", "Lanzacor", "LanzamisATE", "step3DIF", "step3cor1", "step3mis1"), each = 36)
betas <- rep(rep(c(1, 2, 3), each = 12), 7)
gammas <- rep(rep(c(1, 2, 3), each = 4), 21)
direct <- rep(rep(c(1, 2), each = 2), 63)
confounder <- rep(c(1, 2), 126)
bias <- rep(NA, 252)
bias1 <- rep(NA, 36)
bias2 <- rep(NA, 36)
bias3 <- rep(NA, 36)
bias4 <- rep(NA, 36)
bias5 <- rep(NA, 36)
bias6 <- rep(NA, 36)
bias7 <- rep(NA, 36)


Myfinalresults <- data.frame(method, betas, gammas, direct, confounder, bias, bias1, bias2, bias3, bias4, bias5, bias6, bias7)

Myfinalresults$bias1[1:2] <- regadjcorBIAS[1, 1, 1, ]
Myfinalresults$bias1[3:4] <- regadjcorBIAS[1, 1, 2, ]
Myfinalresults$bias1[5:6] <- regadjcorBIAS[1, 2, 1, ]
Myfinalresults$bias1[7:8] <- regadjcorBIAS[1, 2, 2, ]
Myfinalresults$bias1[9:10] <- regadjcorBIAS[1, 3, 1, ]
Myfinalresults$bias1[11:12] <- regadjcorBIAS[1, 3, 2, ]
Myfinalresults$bias1[13:14] <- regadjcorBIAS[2, 1, 1, ]
Myfinalresults$bias1[15:16] <- regadjcorBIAS[2, 1, 2, ]
Myfinalresults$bias1[17:18] <- regadjcorBIAS[2, 2, 1, ]
Myfinalresults$bias1[19:20] <- regadjcorBIAS[2, 2, 2, ]
Myfinalresults$bias1[21:22] <- regadjcorBIAS[2, 3, 1, ]
Myfinalresults$bias1[23:24] <- regadjcorBIAS[2, 3, 2, ]
Myfinalresults$bias1[25:26] <- regadjcorBIAS[3, 1, 1, ]
Myfinalresults$bias1[27:28] <- regadjcorBIAS[3, 1, 2, ]
Myfinalresults$bias1[29:30] <- regadjcorBIAS[3, 2, 1, ]
Myfinalresults$bias1[31:32] <- regadjcorBIAS[3, 2, 2, ]
Myfinalresults$bias1[33:34] <- regadjcorBIAS[3, 3, 1, ]
Myfinalresults$bias1[35:36] <- regadjcorBIAS[3, 3, 2, ]

Myfinalresults$bias2[1:2] <- regadjmisBIAS[1, 1, 1, ]
Myfinalresults$bias2[3:4] <- regadjmisBIAS[1, 1, 2, ]
Myfinalresults$bias2[5:6] <- regadjmisBIAS[1, 2, 1, ]
Myfinalresults$bias2[7:8] <- regadjmisBIAS[1, 2, 2, ]
Myfinalresults$bias2[9:10] <- regadjmisBIAS[1, 3, 1, ]
Myfinalresults$bias2[11:12] <- regadjmisBIAS[1, 3, 2, ]
Myfinalresults$bias2[13:14] <- regadjmisBIAS[2, 1, 1, ]
Myfinalresults$bias2[15:16] <- regadjmisBIAS[2, 1, 2, ]
Myfinalresults$bias2[17:18] <- regadjmisBIAS[2, 2, 1, ]
Myfinalresults$bias2[19:20] <- regadjmisBIAS[2, 2, 2, ]
Myfinalresults$bias2[21:22] <- regadjmisBIAS[2, 3, 1, ]
Myfinalresults$bias2[23:24] <- regadjmisBIAS[2, 3, 2, ]
Myfinalresults$bias2[25:26] <- regadjmisBIAS[3, 1, 1, ]
Myfinalresults$bias2[27:28] <- regadjmisBIAS[3, 1, 2, ]
Myfinalresults$bias2[29:30] <- regadjmisBIAS[3, 2, 1, ]
Myfinalresults$bias2[31:32] <- regadjmisBIAS[3, 2, 2, ]
Myfinalresults$bias2[33:34] <- regadjmisBIAS[3, 3, 1, ]
Myfinalresults$bias2[35:36] <- regadjmisBIAS[3, 3, 2, ]

Myfinalresults$bias3[1:2] <- LanzacorBIAS[1, 1, 1, ]
Myfinalresults$bias3[3:4] <- LanzacorBIAS[1, 1, 2, ]
Myfinalresults$bias3[5:6] <- LanzacorBIAS[1, 2, 1, ]
Myfinalresults$bias3[7:8] <- LanzacorBIAS[1, 2, 2, ]
Myfinalresults$bias3[9:10] <- LanzacorBIAS[1, 3, 1, ]
Myfinalresults$bias3[11:12] <- LanzacorBIAS[1, 3, 2, ]
Myfinalresults$bias3[13:14] <- LanzacorBIAS[2, 1, 1, ]
Myfinalresults$bias3[15:16] <- LanzacorBIAS[2, 1, 2, ]
Myfinalresults$bias3[17:18] <- LanzacorBIAS[2, 2, 1, ]
Myfinalresults$bias3[19:20] <- LanzacorBIAS[2, 2, 2, ]
Myfinalresults$bias3[21:22] <- LanzacorBIAS[2, 3, 1, ]
Myfinalresults$bias3[23:24] <- LanzacorBIAS[2, 3, 2, ]
Myfinalresults$bias3[25:26] <- LanzacorBIAS[3, 1, 1, ]
Myfinalresults$bias3[27:28] <- LanzacorBIAS[3, 1, 2, ]
Myfinalresults$bias3[29:30] <- LanzacorBIAS[3, 2, 1, ]
Myfinalresults$bias3[31:32] <- LanzacorBIAS[3, 2, 2, ]
Myfinalresults$bias3[33:34] <- LanzacorBIAS[3, 3, 1, ]
Myfinalresults$bias3[35:36] <- LanzacorBIAS[3, 3, 2, ]

Myfinalresults$bias4[1:2] <- LanzamisBIAS[1, 1, 1, ]
Myfinalresults$bias4[3:4] <- LanzamisBIAS[1, 1, 2, ]
Myfinalresults$bias4[5:6] <- LanzamisBIAS[1, 2, 1, ]
Myfinalresults$bias4[7:8] <- LanzamisBIAS[1, 2, 2, ]
Myfinalresults$bias4[9:10] <- LanzamisBIAS[1, 3, 1, ]
Myfinalresults$bias4[11:12] <- LanzamisBIAS[1, 3, 2, ]
Myfinalresults$bias4[13:14] <- LanzamisBIAS[2, 1, 1, ]
Myfinalresults$bias4[15:16] <- LanzamisBIAS[2, 1, 2, ]
Myfinalresults$bias4[17:18] <- LanzamisBIAS[2, 2, 1, ]
Myfinalresults$bias4[19:20] <- LanzamisBIAS[2, 2, 2, ]
Myfinalresults$bias4[21:22] <- LanzamisBIAS[2, 3, 1, ]
Myfinalresults$bias4[23:24] <- LanzamisBIAS[2, 3, 2, ]
Myfinalresults$bias4[25:26] <- LanzamisBIAS[3, 1, 1, ]
Myfinalresults$bias4[27:28] <- LanzamisBIAS[3, 1, 2, ]
Myfinalresults$bias4[29:30] <- LanzamisBIAS[3, 2, 1, ]
Myfinalresults$bias4[31:32] <- LanzamisBIAS[3, 2, 2, ]
Myfinalresults$bias4[33:34] <- LanzamisBIAS[3, 3, 1, ]
Myfinalresults$bias4[35:36] <- LanzamisBIAS[3, 3, 2, ]

Myfinalresults$bias5[1:2] <- step3DIFBIAS[1, 1, 1, ]
Myfinalresults$bias5[3:4] <- step3DIFBIAS[1, 1, 2, ]
Myfinalresults$bias5[5:6] <- step3DIFBIAS[1, 2, 1, ]
Myfinalresults$bias5[7:8] <- step3DIFBIAS[1, 2, 2, ]
Myfinalresults$bias5[9:10] <- step3DIFBIAS[1, 3, 1, ]
Myfinalresults$bias5[11:12] <- step3DIFBIAS[1, 3, 2, ]
Myfinalresults$bias5[13:14] <- step3DIFBIAS[2, 1, 1, ]
Myfinalresults$bias5[15:16] <- step3DIFBIAS[2, 1, 2, ]
Myfinalresults$bias5[17:18] <- step3DIFBIAS[2, 2, 1, ]
Myfinalresults$bias5[19:20] <- step3DIFBIAS[2, 2, 2, ]
Myfinalresults$bias5[21:22] <- step3DIFBIAS[2, 3, 1, ]
Myfinalresults$bias5[23:24] <- step3DIFBIAS[2, 3, 2, ]
Myfinalresults$bias5[25:26] <- step3DIFBIAS[3, 1, 1, ]
Myfinalresults$bias5[27:28] <- step3DIFBIAS[3, 1, 2, ]
Myfinalresults$bias5[29:30] <- step3DIFBIAS[3, 2, 1, ]
Myfinalresults$bias5[31:32] <- step3DIFBIAS[3, 2, 2, ]
Myfinalresults$bias5[33:34] <- step3DIFBIAS[3, 3, 1, ]
Myfinalresults$bias5[35:36] <- step3DIFBIAS[3, 3, 2, ]

Myfinalresults$bias6[1:2] <- step3cor1BIAS[1, 1, 1, ]
Myfinalresults$bias6[3:4] <- step3cor1BIAS[1, 1, 2, ]
Myfinalresults$bias6[5:6] <- step3cor1BIAS[1, 2, 1, ]
Myfinalresults$bias6[7:8] <- step3cor1BIAS[1, 2, 2, ]
Myfinalresults$bias6[9:10] <- step3cor1BIAS[1, 3, 1, ]
Myfinalresults$bias6[11:12] <- step3cor1BIAS[1, 3, 2, ]
Myfinalresults$bias6[13:14] <- step3cor1BIAS[2, 1, 1, ]
Myfinalresults$bias6[15:16] <- step3cor1BIAS[2, 1, 2, ]
Myfinalresults$bias6[17:18] <- step3cor1BIAS[2, 2, 1, ]
Myfinalresults$bias6[19:20] <- step3cor1BIAS[2, 2, 2, ]
Myfinalresults$bias6[21:22] <- step3cor1BIAS[2, 3, 1, ]
Myfinalresults$bias6[23:24] <- step3cor1BIAS[2, 3, 2, ]
Myfinalresults$bias6[25:26] <- step3cor1BIAS[3, 1, 1, ]
Myfinalresults$bias6[27:28] <- step3cor1BIAS[3, 1, 2, ]
Myfinalresults$bias6[29:30] <- step3cor1BIAS[3, 2, 1, ]
Myfinalresults$bias6[31:32] <- step3cor1BIAS[3, 2, 2, ]
Myfinalresults$bias6[33:34] <- step3cor1BIAS[3, 3, 1, ]
Myfinalresults$bias6[35:36] <- step3cor1BIAS[3, 3, 2, ]

Myfinalresults$bias7[1:2] <- step3mis1BIAS[1, 1, 1, ]
Myfinalresults$bias7[3:4] <- step3mis1BIAS[1, 1, 2, ]
Myfinalresults$bias7[5:6] <- step3mis1BIAS[1, 2, 1, ]
Myfinalresults$bias7[7:8] <- step3mis1BIAS[1, 2, 2, ]
Myfinalresults$bias7[9:10] <- step3mis1BIAS[1, 3, 1, ]
Myfinalresults$bias7[11:12] <- step3mis1BIAS[1, 3, 2, ]
Myfinalresults$bias7[13:14] <- step3mis1BIAS[2, 1, 1, ]
Myfinalresults$bias7[15:16] <- step3mis1BIAS[2, 1, 2, ]
Myfinalresults$bias7[17:18] <- step3mis1BIAS[2, 2, 1, ]
Myfinalresults$bias7[19:20] <- step3mis1BIAS[2, 2, 2, ]
Myfinalresults$bias7[21:22] <- step3mis1BIAS[2, 3, 1, ]
Myfinalresults$bias7[23:24] <- step3mis1BIAS[2, 3, 2, ]
Myfinalresults$bias7[25:26] <- step3mis1BIAS[3, 1, 1, ]
Myfinalresults$bias7[27:28] <- step3mis1BIAS[3, 1, 2, ]
Myfinalresults$bias7[29:30] <- step3mis1BIAS[3, 2, 1, ]
Myfinalresults$bias7[31:32] <- step3mis1BIAS[3, 2, 2, ]
Myfinalresults$bias7[33:34] <- step3mis1BIAS[3, 3, 1, ]
Myfinalresults$bias7[35:36] <- step3mis1BIAS[3, 3, 2, ]


Myfinalresults$bias[1:36] <- Myfinalresults$bias1[1:36]
Myfinalresults$bias[37:72] <- Myfinalresults$bias2[1:36]
Myfinalresults$bias[73:108] <- Myfinalresults$bias3[1:36]
Myfinalresults$bias[109:144] <- Myfinalresults$bias4[1:36]
Myfinalresults$bias[145:180] <- Myfinalresults$bias5[1:36]
Myfinalresults$bias[181:216] <- Myfinalresults$bias6[1:36]
Myfinalresults$bias[217:252] <- Myfinalresults$bias7[1:36]

Myfinalresults <- Myfinalresults %>%
  select(method, betas, gammas, direct, confounder, bias)


export(Myfinalresults, "Myfinalresults.sav")

write.csv(Myfinalresults, "Myfinalresults.csv", row.names = T)



###############################################################################################

# replication of the results reported in Lanza et al. (2013)

LanzaSelection <- read.csv("LanzaSelection.csv")

# data cleaning
LanzaSelection <- LanzaSelection %>%
  rename(Age = R0000600,
         Household_income = R0217900,
         Maternal_education = R0006500,
         Number_Sib = R0009100,
         Education_aspiration_1 = R0023500)

LanzaSelection$Household_income[LanzaSelection$Household_income < 0] <- NA
LanzaSelection$Maternal_education[LanzaSelection$Maternal_education < 0] <- NA
LanzaSelection$Number_Sib[LanzaSelection$Number_Sib < 0] <- NA
LanzaSelection$Education_aspiration_1[LanzaSelection$Education_aspiration_1 < 0] <- NA

LanzaSelection <- LanzaSelection %>%
  filter(R0015700 == 12 & R0150000 == 0)

LanzaSelection <- LanzaSelection %>%
  rename(Gender = R0214800,
         Race_ethnicity = R0172700,
         Metropolitan_status = R0215100)

LanzaSelection$Race_ethnicity[LanzaSelection$Race_ethnicity < 0] <- NA

LanzaSelection$Language_home <- NA
LanzaSelection$Language_home[LanzaSelection$R0001100 == 0] <- 1
LanzaSelection$Language_home[LanzaSelection$R0001100 == 1 & LanzaSelection$R0001200 == 1] <- 2
LanzaSelection$Language_home[LanzaSelection$R0001100 == 1 & LanzaSelection$R0001200 > 1] <- 3

LanzaSelection$Parent_figure <- car::recode(LanzaSelection$R0001900, "15 = 1; 51 = 1; 53 = 1; 54 = 1; 55 = 1; 91 = 1; 11 = 2; 12 = 2; 21 = 2; 22 = 2; 31 = 2; 33 = 2; 41 = 2; 44 = 2; 80 = NA")

LanzaSelection <- LanzaSelection %>%
  rename(Education_aspiration_2 = R0171800,
         College_Prep = R0019600)
LanzaSelection$Education_aspiration_2[LanzaSelection$Education_aspiration_2 < 0] <- NA
LanzaSelection$College_Prep[LanzaSelection$College_Prep < 0] <- NA

LanzaSelection$college_enroll <- NA
LanzaSelection$college_enroll[LanzaSelection$R0228100 == 0] <- 0
LanzaSelection$college_enroll[LanzaSelection$R0228100 == 1 & LanzaSelection$R0228500 == 0] <- 0
LanzaSelection$college_enroll[LanzaSelection$R0228500 == 1 & LanzaSelection$R0228600 == 13] <- 1
LanzaSelection$college_enroll[LanzaSelection$R0228500 == 1 & LanzaSelection$R0228600 == 14] <- 1

LanzaSelection <- LanzaSelection %>%
  rename(cig_daily_94 = R5052700,
         mar_use_94 = R5053300,
         coc_use_94 = R5053700,
         cra_use_94 = R5054100,
         alcohol_use_past_month_94 = R4979200,
         alcohol_ge_6_past_month_94 = R4979300)

LanzaSelection$Alcohol <- NA

LanzaSelection$alcohol_use_past_month_94[LanzaSelection$R4979100 == 0] <- 0
LanzaSelection$alcohol_ge_6_past_month_94[LanzaSelection$R4979100 == 0] <- 0

LanzaSelection$Alcohol[LanzaSelection$alcohol_use_past_month_94 == 0] <- 1
LanzaSelection$Alcohol[LanzaSelection$alcohol_use_past_month_94 == 1 & LanzaSelection$alcohol_ge_6_past_month_94 == 0] <- 2
LanzaSelection$Alcohol[LanzaSelection$alcohol_use_past_month_94 == 1 & LanzaSelection$alcohol_ge_6_past_month_94 > 0] <- 3

LanzaSelection$cig_daily_94[LanzaSelection$R5052400 == 0] <- 3
LanzaSelection$cig_daily_94[LanzaSelection$R5052400 == 1 & LanzaSelection$R5052500 == 2] <- 3

LanzaSelection$Cigarette <- NA
LanzaSelection$Cigarette[LanzaSelection$cig_daily_94 == 1] <- 1
LanzaSelection$Cigarette[LanzaSelection$cig_daily_94 == 2] <- 2
LanzaSelection$Cigarette[LanzaSelection$cig_daily_94 == 3] <- 3

LanzaSelection$CigaretteDummy <- NA
LanzaSelection$CigaretteDummy[LanzaSelection$Cigarette == 1] <- 1
LanzaSelection$CigaretteDummy[LanzaSelection$Cigarette == 2] <- 1
LanzaSelection$CigaretteDummy[LanzaSelection$Cigarette == 3] <- 0

LanzaSelection$mar_use_94[LanzaSelection$R5053100 == 0] <- 4
LanzaSelection$Marijuana <- NA
LanzaSelection$Marijuana[LanzaSelection$mar_use_94 == 0] <- 1
LanzaSelection$Marijuana[LanzaSelection$mar_use_94 > 0] <- 0


LanzaSelection$coc_use_94[LanzaSelection$R5053500 == 0] <- 4
LanzaSelection$Cocaine <- NA
LanzaSelection$Cocaine[LanzaSelection$coc_use_94 == 0] <- 1
LanzaSelection$Cocaine[LanzaSelection$coc_use_94 > 0] <- 0

LanzaSelection$cra_use_94[LanzaSelection$R5053900 == 0] <- 4
LanzaSelection$Crack <- NA
LanzaSelection$Crack[LanzaSelection$cra_use_94 == 0] <- 1
LanzaSelection$Crack[LanzaSelection$cra_use_94 > 0] <- 0

LanzaSelection$CrackCocaine <- NA
LanzaSelection$CrackCocaine[LanzaSelection$Cocaine == 0 & LanzaSelection$Crack == 0] <- 0
LanzaSelection$CrackCocaine[LanzaSelection$Cocaine == 1] <- 1
LanzaSelection$CrackCocaine[LanzaSelection$Crack == 1] <- 1

LanzaSelection$CSweight79 <- LanzaSelection$R0216101
LanzaSelection$weight79 <- LanzaSelection$R0216100
LanzaSelection$CSweight80 <- LanzaSelection$R0405201
LanzaSelection$weight80 <- LanzaSelection$R0405200
LanzaSelection$CSweight94 <- LanzaSelection$R5080401
LanzaSelection$weight94 <- LanzaSelection$R5080400

LanzaSelection$ID <- c(1:1092)

LanzaSelection <- LanzaSelection %>%
  select(ID, Age, Gender, Race_ethnicity, Household_income, Number_Sib, Language_home, Maternal_education,  
         Education_aspiration_1, Education_aspiration_2, Parent_figure, Metropolitan_status,
         College_Prep, college_enroll, Alcohol, Cigarette, CigaretteDummy, Marijuana, Cocaine, Crack, 
         CrackCocaine, CSweight79, weight79, CSweight80, weight80, CSweight94, weight94)

export(LanzaSelection, "LanzaSelectionClean.sav")

LanzaSelectionClean <- read.spss("LanzaSelectionClean.sav", to.data.frame = T)

# multiple imputation of missing values
LanzaSelectionClean <- LanzaSelectionClean %>%
  mutate(Gender = as.factor(Gender)) %>%
  mutate(Race_ethnicity = as.factor(Race_ethnicity)) %>%
  mutate(Number_Sib = as.factor(Number_Sib)) %>%
  mutate(Language_home = as.factor(Language_home)) %>%
  mutate(Education_aspiration_2 = as.factor(Education_aspiration_2)) %>%
  mutate(Parent_figure = as.factor(Parent_figure)) %>%
  mutate(Metropolitan_status = as.factor(Metropolitan_status)) %>%
  mutate(College_Prep = as.factor(College_Prep)) %>%
  mutate(college_enroll = as.factor(college_enroll)) %>%
  mutate(Alcohol = as.factor(Alcohol)) %>%
  mutate(Cigarette = as.factor(Cigarette)) %>%
  mutate(CigaretteDummy = as.factor(CigaretteDummy)) %>%
  mutate(Marijuana = as.factor(Marijuana)) %>%
  mutate(Cocaine = as.factor(Cocaine)) %>%
  mutate(Crack = as.factor(Crack)) %>%
  mutate(CrackCocaine = as.factor(CrackCocaine))

LanzaSelectionClean <- LanzaSelectionClean %>%
  select(Age, Gender, Race_ethnicity, Household_income, Number_Sib, Language_home, Maternal_education,  
         Education_aspiration_1, Education_aspiration_2, Parent_figure, Metropolitan_status,
         College_Prep, college_enroll, Alcohol, Cigarette, CigaretteDummy, Marijuana, Cocaine, Crack, 
         CrackCocaine)

init = mice(LanzaSelectionClean, maxit=0) 
meth = init$method
predM = init$predictorMatrix

predM[, 1] <- 0
predM[, 22] <- 0
predM[, 23] <- 0
predM[, 24] <- 0
predM[, 25] <- 0
predM[, 26] <- 0
predM[, 27] <- 0

meth[c("Alcohol")]="" 
meth[c("Cigarette")]="" 
meth[c("CigaretteDummy")]="" 
meth[c("Marijuana")]="" 
meth[c("Cocaine")]="" 
meth[c("Crack")]="" 
meth[c("CrackCocaine")]="" 
meth[c("college_enroll")]="" 

LanzaImputed <- mice(LanzaSelectionClean, method=meth, predictorMatrix=predM, m=5, seed = 4753262)

ImputedData <- complete(LanzaImputed, action = "long")

ImputedData <- ImputedData %>%
  rename(id = .id,
         imp = .imp)

ImputedData_1 <- ImputedData %>%
  filter(imp == 1)
ImputedData_2 <- ImputedData %>%
  filter(imp == 2)
ImputedData_3 <- ImputedData %>%
  filter(imp == 3)
ImputedData_4 <- ImputedData %>%
  filter(imp == 4)
ImputedData_5 <- ImputedData %>%
  filter(imp == 5)

# estimating IPW
MyGLM_1 <- glm(college_enroll ~ Age + Gender + Race_ethnicity + Household_income + Number_Sib + Language_home + 
                 Maternal_education + Education_aspiration_1 + Education_aspiration_2 + Parent_figure +
                 Metropolitan_status + College_Prep, family = "binomial", data = ImputedData_1)
PS_est_1 <- MyGLM_1$fitted.values
ImputedData_1$PS <- PS_est_1
ImputedData_1$IPW <- ifelse(ImputedData_1$college == 1, 1/ImputedData_1$PS, 1/(1-ImputedData_1$PS))

MyGLM_2 <- glm(college_enroll ~ Age + Gender + Race_ethnicity + Household_income + Number_Sib + Language_home + 
                 Maternal_education + Education_aspiration_1 + Education_aspiration_2 + Parent_figure +
                 Metropolitan_status + College_Prep, family = "binomial", data = ImputedData_2)
PS_est_2 <- MyGLM_2$fitted.values
ImputedData_2$PS <- PS_est_2
ImputedData_2$IPW <- ifelse(ImputedData_2$college == 1, 1/ImputedData_2$PS, 1/(1-ImputedData_2$PS))

MyGLM_3 <- glm(college_enroll ~ Age + Gender + Race_ethnicity + Household_income + Number_Sib + Language_home + 
                 Maternal_education + Education_aspiration_1 + Education_aspiration_2 + Parent_figure +
                 Metropolitan_status + College_Prep, family = "binomial", data = ImputedData_3)
PS_est_3 <- MyGLM_3$fitted.values
ImputedData_3$PS <- PS_est_3
ImputedData_3$IPW <- ifelse(ImputedData_3$college == 1, 1/ImputedData_3$PS, 1/(1-ImputedData_3$PS))

MyGLM_4 <- glm(college_enroll ~ Age + Gender + Race_ethnicity + Household_income + Number_Sib + Language_home + 
                 Maternal_education + Education_aspiration_1 + Education_aspiration_2 + Parent_figure +
                 Metropolitan_status + College_Prep, family = "binomial", data = ImputedData_4)
PS_est_4 <- MyGLM_4$fitted.values
ImputedData_4$PS <- PS_est_4
ImputedData_4$IPW <- ifelse(ImputedData_4$college == 1, 1/ImputedData_4$PS, 1/(1-ImputedData_4$PS))

MyGLM_5 <- glm(college_enroll ~ Age + Gender + Race_ethnicity + Household_income + Number_Sib + Language_home + 
                 Maternal_education + Education_aspiration_1 + Education_aspiration_2 + Parent_figure +
                 Metropolitan_status + College_Prep, family = "binomial", data = ImputedData_5)
PS_est_5 <- MyGLM_5$fitted.values
ImputedData_5$PS <- PS_est_5
ImputedData_5$IPW <- ifelse(ImputedData_5$college == 1, 1/ImputedData_5$PS, 1/(1-ImputedData_5$PS))

ImputedData_1$Alcohol <- LanzaSelectionClean$Alcohol
ImputedData_1$CigaretteDummy <- LanzaSelectionClean$CigaretteDummy
ImputedData_1$Marijuana <- LanzaSelectionClean$Marijuana
ImputedData_1$CrackCocaine <- LanzaSelectionClean$CrackCocaine

ImputedData_2$Alcohol <- LanzaSelectionClean$Alcohol
ImputedData_2$CigaretteDummy <- LanzaSelectionClean$CigaretteDummy
ImputedData_2$Marijuana <- LanzaSelectionClean$Marijuana
ImputedData_2$CrackCocaine <- LanzaSelectionClean$CrackCocaine

ImputedData_3$Alcohol <- LanzaSelectionClean$Alcohol
ImputedData_3$CigaretteDummy <- LanzaSelectionClean$CigaretteDummy
ImputedData_3$Marijuana <- LanzaSelectionClean$Marijuana
ImputedData_3$CrackCocaine <- LanzaSelectionClean$CrackCocaine

ImputedData_4$Alcohol <- LanzaSelectionClean$Alcohol
ImputedData_4$CigaretteDummy <- LanzaSelectionClean$CigaretteDummy
ImputedData_4$Marijuana <- LanzaSelectionClean$Marijuana
ImputedData_4$CrackCocaine <- LanzaSelectionClean$CrackCocaine

ImputedData_5$Alcohol <- LanzaSelectionClean$Alcohol
ImputedData_5$CigaretteDummy <- LanzaSelectionClean$CigaretteDummy
ImputedData_5$Marijuana <- LanzaSelectionClean$Marijuana
ImputedData_5$CrackCocaine <- LanzaSelectionClean$CrackCocaine

export(ImputedData_1, "LanzaImputed_imp1.sav")
export(ImputedData_2, "LanzaImputed_imp2.sav")
export(ImputedData_3, "LanzaImputed_imp3.sav")
export(ImputedData_4, "LanzaImputed_imp4.sav")
export(ImputedData_5, "LanzaImputed_imp5.sav")

ImputedData_v2 <- rbind(ImputedData_1, ImputedData_2)
ImputedData_v2 <- rbind(ImputedData_v2, ImputedData_3)
ImputedData_v2 <- rbind(ImputedData_v2, ImputedData_4)
ImputedData_v2 <- rbind(ImputedData_v2, ImputedData_5)

export(ImputedData_v2, "LanzaImputed_v2.sav")

ImputedData_1$PSt <- ifelse(ImputedData_1$college_enroll == 1, ImputedData_1$PS, 0)
ImputedData_1$PSnt <- ifelse(ImputedData_1$college_enroll == 0, ImputedData_1$PS, 0)

ImputedData_1$PSt[ImputedData_1$PSt ==0] <- NA
ImputedData_1$PSnt[ImputedData_1$PSnt ==0] <- NA

# Figure 3. Overlap of propoesnity scores
png("PSoverlap2.png", width = 3100, height = 2100, res = 300)
ggplot(ImputedData_1, aes(x=PS) ) +
  geom_histogram( aes(x = PSt, y = ..density..), fill="#155F83FF" ) +
  geom_label( aes(x=.2, y=3, label="treated"), color="#155F83FF") + 
  geom_histogram( aes(x = PSnt, y = -..density..), fill= "#C16622FF") +
  geom_label( aes(x=.21, y=-4, label="untreated"), color="#C16622FF") +
  xlab("propensity score") +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15))
dev.off()

# Table 3. Balance of treated vs. untreated
xvars <- c("Age", "Gender", "Race_ethnicity", "Household_income", "Number_Sib", "Language_home", 
           "Maternal_education", "Education_aspiration_1", "Education_aspiration_2", 
           "Parent_figure", "Metropolitan_status", "College_Prep")

NWDesign <- svydesign(ids = ~1, data = ImputedData_1)
table1.nw <- svyCreateTableOne(vars = xvars, strata = "college_enroll", data = NWDesign, test = T)
t1 <- print(table1.nw, smd = T)

IPWDesign <- svydesign(ids = ~1, data = ImputedData_1, weights = ImputedData_1$IPW)
table1.ipw <- svyCreateTableOne(vars = xvars, strata = "college_enroll", data = IPWDesign, test = T)
t2 <- print(table1.ipw, smd = T)
