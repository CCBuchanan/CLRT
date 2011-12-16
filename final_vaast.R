#VAAST CLR method
#This utilizes the VAAST composite likelihood method employed by Yandall et al (2011).
#The idea is to evaluate if there is an increased burden of rare variants in affected individuals
#in a case control study.  
#Data: See write-up for example of data format.

#The main method is called 'clrt'
#Required inputs: data file which will become a data table
####################################################
#calculates the null and alternate likelihoods
L.null <- function(af_all,total_var,n,xi){
  return((af_all^xi)*(1-af_all)^((2*total_var*n)-xi))
}

L.alt <- function(af_u, xu,n_u, total_var,af_a,xa, n_a){
  return((af_u^xu)*(1-af_u)^((2*total_var*n_u)-xu)*(af_a^xa)*(1-af_a)^((2*total_var*n_a)-xa))
}
####################################################
#sums the number of minor alleles in each bin for AFFECTED individuals
suma <- function(exon_bins){
  xa <- rep(0,2)
  total_row <- length(exon_bins$ID)
  for (i in 3:length(exon_bins)){
    count <-0
    for (j in 2:total_row){
      
      if (exon_bins[j,2]>1){
      #this means the row is an affected individual
      count <- count+exon_bins[j,i]
      }
    }
    xa <- c(xa,count)
  }
  return(xa)
}

###################################################
#sums the number of minor alleles in each bin for ALL individuals
sumi <- function(exon_bins){
  xi <- rep(0,2)
  #show(exon_bins[2,2])
  for (i in 3:length(exon_bins)){
    col_total <- sum(exon_bins[, i]) - exon_bins[1,i]
    #show(col_total)
    xi <- c(xi,col_total)
  }
  return(xi)
}
###################################################
#returns the number of unaffected individuals in the data according to the status column
#by default, 1 is a control and 2 is a case
n_unaff <- function(exon_bins){
  count2 <- 0
  total_row <- length(exon_bins$ID)
  for (j in 2:total_row){
    if (exon_bins[j,2]>1){
      #this means the row is an affected individual
      count2 <- count2+1
    }
  }
  return(count2)
}
####################################################
#returns the total number of variant sites for each bin
total_variant_sites <- function(exon_bins){
  xv <- rep(0,2)
  total_row <- length(exon_bins$ID)
  for (i in 3:length(exon_bins)){
    count <-0
    for (j in 2:total_row){
      
      if (exon_bins[j,i]>0){
        count <- count+1
      }
    }
    xv[i]<-count
  }
  return(xv)
}
#######################################
#returns the number of variant sites by reviewing only cases
aff_variants_sites <- function(exon_bins){
  xv <- rep(0,2)
  total_row <- length(exon_bins$ID)
  for (i in 3:length(exon_bins)){
    count <-0
    for (j in 2:total_row){
      if (exon_bins[j,2]>1){
        if (exon_bins[j,i]>0){
            count <- count+1
        }
      }
    }
    xv[i]<-count
  }
  return(xv)
}
#####################################################
#main function
#composite likelihood ratio test similar to CLT used in VAAST method
clrt <- function(name){
#data should be in csv with header row
exon_bins <- read.csv(name)
#opens dataset
View(exon_bins)

#total number of bins being tested (-2 represents the status and ID column)
bin <- length(exon_bins)-2

#bonferroni corrected p-value which corrects for multiple testing
bonferroni <- 0.05/bin

#Number of affected and unaffected individuals in the file 
n_u <- n_unaff(exon_bins)
n_tot <- length(exon_bins$ID)-1
n_a <- (n_tot-n_u)

#total number of variant SITES in each bin
#NOTE: sites refer to an actual physical variant location, at each location (because there are 2 chromosomes)
#you can zero, one, or two  minor alleles
total_variants <- total_variant_sites(exon_bins)
aff_variants <- aff_variants_sites(exon_bins)
unaff_variants <- total_variants - aff_variants

#x is the number of copies of the minor allele for each group of individuals 
#(can be greater than number of variant sites)
xi <- sumi(exon_bins)
xa <- suma(exon_bins)
xu <- xi-xa
#could test if xu is negative here and throw a format error
#show(xu)

#m is the number of collapsing categories, ie. greater than one if using functional info
m <- 1
#l is the total number of variant sites in a given m
l <- exon_bins[1,]

#allele frequency
#could use this as an input a matrix with 3 rows, header is same as big file
#could do a test to be sure the header rows match
af_all <- total_variants/n_tot
show(aff_variants)
show(n_a)
af_a <- aff_variants/n_a
af_u <- unaff_variants/n_u


#1 to bin will do this calculation for the entire file
lambda <- rep(0,2)
p_val <- rep(0,2)
sig <- mat.or.vec(1,0)
sig_chi <- mat.or.vec(1,0)

#i iterates over each column (which are bins)
for (i in 3:length(exon_bins)){
  likeli_null <- L.null(af_all[i],total_variants[i],n_tot,xi[i])
  likeli_alt <- L.alt(af_u[i],xu[i],n_u,total_variants[i],af_a[i],xa[i],n_a)
  lambda[i] <- log2(likeli_null/likeli_alt)
  if (is.na(lambda[i])==T){
    lambda[i] <- 0
  }
  #calculates p-value from chi-square distribution
  p_val[i] <- dchisq(lambda[i],df=bin)
  #creates second table of only bins which are below the bonferronni corrected p-value
  if (p_val[i] < bonferroni && p_val[i] > 0 ){
    Significant_bin <- c(sig, colnames(exon_bins)[i])
    P_value <- c(sig_chi,p_val[i])
  }
}

whole_list <- cbind(colnames(exon_bins),lambda, p_val, xa, xu, xi,af_a, af_u, af_all)
final_ans <- cbind(Significant_bin,P_value)
show(whole_list)
return(final_ans)
}  

#Testing:
#I performed a great deal of testing but because the file format is so rigid, the best way to test is to look at the output
#from the first table and then the data to be sure the values in the table match those from the dataset.  This is easy to do in excel
#when the dataset is small.


