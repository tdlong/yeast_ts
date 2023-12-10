args <- commandArgs(TRUE)
pool <- as.character(args[1])
folder <- as.character(args[2])
SNPtable <- as.character(args[3])
foundernames <- as.character(args[4])
NSNPs <- as.numeric(args[5])

library(limSolve)
source("scripts/haplotyper.limSolve.simple.code.R")
runscan(pool, folder, SNPtable, foundernames, NSNPs)

