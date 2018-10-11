rm(list = ls())


#Following code is a modification of the one written for the analysis of acetyl and alcohol groups in the polysaccharide (Original: mass_data.R).
#The aim is to identify genes signifiant in regression analysis with respect to the indicidual polysaccharides.

library(dplyr)
library(tidyr)
library(reshape2)
library(MASS)
library(ggplot2)
library(dplyr)
library(plotly)
library(lme4)

#Load the data file in csv format
pan_genes <- read.csv("pangenome_gene_presence_absence.csv")
#Transform the data
pan_genes <- t(pan_genes)
#Assign gene names as column names
colnames(pan_genes) <- pan_genes[1,]
#Remove the row containing gene names
pan_genes <- pan_genes[-1,]
#Convert pan_genes matrix into a data frame
pan_genes <- data.frame(pan_genes)

#Change rownames to just the sample names
#Import data for Israel samples
isr_seqs <- read.csv("WGS_Isr_corrected.csv")
#Select first four columns of the data frame
isr_seqs <- isr_seqs[,1:4]
#Add string "Pool_DMW_" to samples names so that it matches with pangenome output
isr_seqs$corrected.well <- paste("Pool_DMW_",isr_seqs$corrected.well,sep="")
#Assign value of 0 to miscalls
isr_seqs$miscall[which(isr_seqs$miscall == 1)] = 0
#Assign value 1 to samples for which results from PneumoCAT matched
isr_seqs$miscall[which(is.na(isr_seqs$miscall)==TRUE)] = 1

#Import data for Hungary samples
hun_seqs <- read.csv("WGS_hun_corrected.csv")
#Select first four columns of the data frame
hun_seqs <- hun_seqs[,1:4]
#Assign column names same as isr_seqs
colnames(hun_seqs) <- colnames(isr_seqs)
#Add string "Pool_hung_" to samples names so that it matches with pangenome output
hun_seqs$corrected.well <- paste("Pool_hung_",hun_seqs$corrected.well,sep="")

#Combine data for Israel and Hungary samples
seq_info <- rbind(isr_seqs, hun_seqs)

#In pan_genes data frame, make a new column "corrected.wells" that will help in matching this column with the seq_info data set
pan_genes$corrected.well <- rownames(pan_genes)
#In pan_genes, assign numerical row names
rownames(pan_genes) <- 1:nrow(pan_genes)

#Combine seq_info and pan_genes datasets
pan_genes_st <- merge(x = seq_info, y = pan_genes, by = "corrected.well", all.x = TRUE)

#Remove rows with NA
pan_genes_st <- na.omit(pan_genes_st)
#Remove samples for which miscall values are 0 i.e. serotype from pneumocat didn't match with the reported one
pan_genes_st <- pan_genes_st[-which(pan_genes_st$miscall==0),]
#Remove the NT samples. These are also the ones for which pneumocat failed
pan_genes_st <- pan_genes_st[-which(pan_genes_st$pneumocat.serotype == "Failed"),]
#Remove the samples for which assigned serotype is NT
pan_genes_st <- pan_genes_st[-which(pan_genes_st$assigned.serotype == "NT"),]
#Convert assigned serotype column to a character datatype
pan_genes_st$assigned.serotype <- as.character(pan_genes_st$assigned.serotype)
#Fix the serotype names for "19A"
pan_genes_st$assigned.serotype[which(pan_genes_st$pneumocat.serotype == "19A")] <- "19A"
#Fix serotype for Pool_hung_H11
pan_genes_st$assigned.serotype[which(pan_genes_st$corrected.well == "Pool_hung_H11")] <- "38"
#Fix serotype for Pool_hung_F8
pan_genes_st$assigned.serotype[which(pan_genes_st$corrected.well == "Pool_hung_F8")] <- "24B"
#Fix serotype for Pool_hung_F8
pan_genes_st$assigned.serotype[which(pan_genes_st$corrected.well == "Pool_hung_F8")] <- "24B"

#Add ps composition to pan_genes dataset based on the serotype definitions from pneumocat
#Import dataset for polysaccharide composition
ps_composition <- read.csv("ps_composition.csv")
#Combine pan_genes_st dataset with ps_composition
pan_genes_st_ps <- merge(x = pan_genes_st, y = ps_composition, by.x = "assigned.serotype", by.y = "Serotype", all.x = TRUE)
pan_genes_st_ps <- na.omit(pan_genes_st_ps)

fdat2 <- pan_genes_st_ps[,5:ncol(pan_genes_st_ps)]
#Make ps composition in fdat2 a separate dataset
ps_fdat2 <- fdat2[,(ncol(fdat2)-23):ncol(fdat2)]
#Assign remaining dataset to fdat2
fdat2 <- fdat2[,1:(ncol(fdat2)-24)]
#Logistic Regression for Acetyl
reg.ace.nb1 <- vector("list", ncol(fdat2)) 
aic.ace.nb1 <- 1: ncol(fdat2)
pred.ace.nb1<-matrix(nrow = nrow(fdat2), ncol = ncol(fdat2) )
for (i in 1:ncol(fdat2)){
  if(nlevels(as.factor(fdat2[,i])) > 1){
    reg.ace.nb1[[i]]<-glm(as.factor(ps_fdat2[,1]) ~ as.factor(fdat2[,i]), family = binomial)
    aic.ace.nb1[[i]]<-reg.ace.nb1[[i]]$aic
    pred.ace.nb1[,i]<-predict(reg.ace.nb1[[i]])
  }
  else{
    #AIC values cannot be calculated since there is only one level
    reg.ace.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    aic.ace.nb1[[i]]<-1000
    #Predicted value will be same as existing value
    pred.ace.nb1[,i]<-fdat2[1,i]
  }
}
#Compute weights for each position
wgt_acetyl <- exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1)))/sum(exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
wgt1_acetyl <- wgt_acetyl
wgt1_acetyl[wgt1_acetyl == wgt1_acetyl[1]] <- 0
pos1 <- colnames(fdat2)
#Make a dataframe containing positions and weights
plot_data_ace <- data.frame(pos1,wgt1_acetyl) 
names(plot_data_ace) <- c("Position","Weight")
plot_data_ace$Product <- ps_gene$Product[ps_gene$Position %in% plot_data_ace$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(plot_data_ace,"plot_data_acetyl.csv")

#Plot weights with respect to position
plot_data_ace$Position <- factor(plot_data_ace$Position, as.character(plot_data_ace$Position))
ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight), group=1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
plot_acetyl <- ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight^0.005)^10, group=1, text=plot_data_ace$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(plot_acetyl, tooltip = c("x","text"))



#===================================================================================================================#

#Make a dataframe with accession_number, serotype, and sequence type (ST)
#Serotype and ST information can be accessed using the accession number
st_accession <- mass_S1[1:2,]
st_accession[,1:3] <- NULL
st_accession <- as.data.frame(t(st_accession))
colnames(st_accession) <- c("Serotype","ST")
st_accession <- st_accession[-1,]
length(unique(rownames(st_accession))) #616 accession numbers
length(unique(st_accession$Serotype)) #31 serotypes
length(unique(st_accession$ST)) #140 STs

#Add pathway details to the data frame
pathway <- read.csv("pathway.csv")

#Function to convert factors into character data type
fac2char <- function(df_name){
  df_name[] <- lapply(df_name, as.character)
  df_name <- as.data.frame(df_name)
  return(df_name)
}
st_accession <- fac2char(st_accession)
pathway <- fac2char(pathway)
st_accession$acetyl <- pathway$acetyl[match(st_accession$Serotype, pathway$Serotype)]
#Serotype 15B and 15C both have acetyl so provide value 1
st_accession$acetyl[st_accession$Serotype == "15B/C"] <- "1"

st_accession$alcohol <- pathway$alcohol[match(st_accession$Serotype, pathway$Serotype)]
#Serotype 15B and 15C both have alcohol so provide value 1
st_accession$alcohol[st_accession$Serotype == "15B/C"] <- "1"

#Write a csv file containing the accession number, serotype, ST, and binary for acetyl and alcohol
write.csv(st_accession, "st_accession.csv")

#Make copy to mass_S1 into a new variable mass_data
#The third column X.2 contains position number that is unique. 
#This will be unique id, and gene symbols and names are deleted (Refer to ps_gene data frame).
mass_data <- mass_S1
mass_data[,1:2] <- NULL
mass_data <- mass_data[-c(1:2),]
mass_data[,2] <- NULL
colnames(mass_data)[1] <- "Position" 
#Assign position numbers as row names and remove the position column
rownames(mass_data) <- mass_data$Position 
mass_data$Position <- NULL

#Transpose data
tr_mass_data <- as.data.frame(t(mass_data))

#The data type for the columns in the data frame are factors containing allelic profile as levels.
#Convert factors to character data type (Refer fac2char function)
tr_mass_data <- fac2char(tr_mass_data)

#Replace cells not containg "X" with 1
tr_mass_data[tr_mass_data!="X"] <- "1"
#Replace cells containing "X" with 0
tr_mass_data[tr_mass_data=="X"] <- "0"

#Columns to be removed
to_remove <- c("300883","305006","306448","307152","308550","311225","312203","313437","314830","316014","318770","319640","320246","321361","323780")
#Add acetyl and alcohol columns to the main data frame
#fdat1 <- cbind(tr_mass_data, st_accession[,3:4])
#Remove above mentioned accession numbers from fdat1
fdat1 <- tr_mass_data[,!(names(tr_mass_data) %in% to_remove)]

#Write the fdat1 data frame into a csv file.
#This file contains gene product presence/absence information in the binary format
#The positions that need to be removed have already been removed (expand the details)
write.csv(fdat1,"mass_data.csv")

#Serotypes for which acetyl column is NA
na_acetyl <- unique(st_accession$Serotype[which(is.na(st_accession$acetyl)==TRUE)])
#Serotypes for which alcohol column is NA
na_alcohol <- unique(st_accession$Serotype[which(is.na(st_accession$alcohol)==TRUE)])

#Remove the serotypes with missing values.
#Serotypes with missing values are same for both acetyl and alcohol in the current dataset.
#If dataset changes, check values for na_acetyl and na_alcohol before proceeding
#Serotypes excluded in this analysis: "38"  "23A" "23B" "7C"  "NT" 
#Extract row names from st_accession dataframe that contain missing values (NA)
na_row_acetyl <- rownames(st_accession)[which(is.na(st_accession$acetyl)==TRUE)]
#Remove rows from fdat1 and st_accession dataframe
fdat2 <- fdat1[!rownames(fdat1) %in% na_row_acetyl,]
st_access1 <- st_accession[!rownames(st_accession) %in% na_row_acetyl,]

#Function to convert characters into factor data type
#char2fac <- function(df_name){
#  bin_levels <- unique(unlist( lapply( df_name , levels )))
  #Set binary "0" and "1" levels for each column even though only one
  #level of either 0 or 1 might be present throughout the column
#  df_name <- data.frame(lapply(df_name , factor , levels = bin_levels))
#  return(df_name)
#}

#Convert fdat2 and st_access1 to factor
#fdat2 <- char2fac1(fdat2)
#st_access1 <- char2fac1(st_access1)

#The following if_factor variable counts the number of levels for each column in a dataframe
if_factor <- sapply(fdat2, function(x) nlevels(as.factor(x)))
#To test the number of columns with more than one levels
length(which(if_factor>1))
#To test the number of columns with only one level
length(which(if_factor==1))

#Logistic Regression for Acetyl
reg.ace.nb1 <- vector("list", ncol(fdat2)) 
aic.ace.nb1 <- 1: ncol(fdat2)
pred.ace.nb1<-matrix(nrow = nrow(fdat2), ncol = ncol(fdat2) )
for (i in 1:ncol(fdat2)){
  if(nlevels(as.factor(fdat2[,i])) > 1){
    reg.ace.nb1[[i]]<-glm(as.factor(st_access1$acetyl) ~ as.factor(fdat2[,i]), family = binomial)
    aic.ace.nb1[[i]]<-reg.ace.nb1[[i]]$aic
    pred.ace.nb1[,i]<-predict(reg.ace.nb1[[i]])
  }
  else{
    #AIC values cannot be calculated since there is only one level
    reg.ace.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    aic.ace.nb1[[i]]<-1000
    #Predicted value will be same as existing value
    pred.ace.nb1[,i]<-fdat2[1,i]
  }
}
#Compute weights for each position
wgt_acetyl <- exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1)))/sum(exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
wgt1_acetyl <- wgt_acetyl
wgt1_acetyl[wgt1_acetyl == wgt1_acetyl[1]] <- 0
pos1 <- colnames(fdat2)
#Make a dataframe containing positions and weights
plot_data_ace <- data.frame(pos1,wgt1_acetyl) 
names(plot_data_ace) <- c("Position","Weight")
plot_data_ace$Product <- ps_gene$Product[ps_gene$Position %in% plot_data_ace$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(plot_data_ace,"plot_data_acetyl.csv")

#Plot weights with respect to position
plot_data_ace$Position <- factor(plot_data_ace$Position, as.character(plot_data_ace$Position))
ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight), group=1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
plot_acetyl <- ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight^0.005)^10, group=1, text=plot_data_ace$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(plot_acetyl, tooltip = c("x","text"))

#One of the positions has disproportainately high value of weight
#plot_data$Position[which(plot_data$Weight==max(plot_data$Weight))] #1023758

#Remove the weight for that position
#plot_data1 <- plot_data[-which(plot_data$Position==1023758),]

#Replot
#plot_data1$Position <- factor(plot_data1$Position, as.character(plot_data1$Position))
#ggplot(data=plot_data1, aes(x=plot_data1$Position, y=(plot_data1$Weight^0.005)^10, group=1)) +
#  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
#  scale_x_discrete(breaks = plot_data1$Position[which(plot_data1$Weight!=0)]) +
#  labs(x = "Position", y = "Weight (adjusted)")

#Logistic Regression for alcohol
reg.alc.nb1 <- vector("list", ncol(fdat2)) 
aic.alc.nb1 <- 1: ncol(fdat2)
pred.alc.nb1<-matrix(nrow = nrow(fdat2), ncol = ncol(fdat2) )
for (i in 1:ncol(fdat2)){
  if(nlevels(as.factor(fdat2[,i])) > 1){
    reg.alc.nb1[[i]]<-glm(as.factor(st_access1$alcohol) ~ as.factor(fdat2[,i]), family = binomial)
    aic.alc.nb1[[i]]<-reg.alc.nb1[[i]]$aic
    pred.alc.nb1[,i]<-predict(reg.alc.nb1[[i]])
  }
  else{
    #AIC values cannot be calculated since there is only one level
    reg.alc.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    aic.alc.nb1[[i]]<-1000
    #Predicted value will be same as existing value
    pred.alc.nb1[,i]<-fdat2[1,i]
  }
}
#Compute weights for each position
wgt_alcohol <- exp(-0.5*(aic.alc.nb1-min(aic.alc.nb1)))/sum(exp(-0.5*(aic.alc.nb1-min(aic.alc.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
wgt1_alcohol <- wgt_alcohol
wgt1_alcohol[wgt1_alcohol == wgt1_alcohol[1]] <- 0
#Make a dataframe containing positions and weights
plot_data_alc <- data.frame(pos1,wgt1_alcohol) #Refer to plot_data_ace for pos1
names(plot_data_alc) <- c("Position","Weight")
plot_data_alc$Product <- ps_gene$Product[ps_gene$Position %in% plot_data_alc$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(plot_data_alc,"plot_data_alcohol.csv")

#Plot weights with respect to position
plot_data_alc$Position <- factor(plot_data_alc$Position, as.character(plot_data_alc$Position))
ggplot(data=plot_data_alc, aes(x=plot_data_alc$Position, y=(plot_data_alc$Weight), group=1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_alc$Position[which(plot_data_alc$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
plot_alcohol <- ggplot(data=plot_data_alc, aes(x=plot_data_alc$Position, y=(plot_data_alc$Weight^0.005)^10, group=1, text=plot_data_alc$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_alc$Position[which(plot_data_alc$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(plot_alcohol, tooltip = c("x","text"))

#Mixed model with random effects
fdat3 <- fdat2

#Model for Acetyl
lmereg.ace.nb1 <- vector("list", ncol(fdat3)) 
lmeaic.ace.nb1 <- 1:(ncol(fdat3))
lmepred.ace.nb1<-matrix(nrow = nrow(fdat3), ncol = ncol(fdat3))

#fdat3[nrow(fdat3)+1,] <- sapply(fdat3, function(x) sum(as.numeric(x)) ) 

#fdat3 <- dplyr::filter(fdat3, fdat3[nrow(fdat3)+1,] < nrow(fdat3))

for (i in 1:ncol(fdat3)){
  if(nlevels(as.factor(fdat3[,i])) > 1){
    tryCatch({lmereg.ace.nb1[[i]]<-glmer(as.factor(st_access1$acetyl) ~ as.factor(fdat3[,i]) + 
                                           (1|as.factor(rownames(st_access1))), family = binomial)}, error = function(err){})
    if(is.null(lmereg.ace.nb1[[i]]) == FALSE){
      lmeaic.ace.nb1[[i]]<-extractAIC(lmereg.ace.nb1[[i]])[2]
      #lmepred.alc.nb1[,i]<-predict(lmereg.alc.nb1[[i]])
    }
    else {
      lmeaic.ace.nb1[[i]]<-3000
      #lmepred.alc.nb1[,i]<-"NA"
    }
  }
  else{
    #AIC values cannot be calculated since there is only one level
    lmereg.ace.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    lmeaic.ace.nb1[i]<-1000
    #Predicted value will be same as existing value
    #lmepred.alc.nb1[,i]<-fdat3[1,i]
  }
}

#Compute weights for each position
lmewgt_acetyl <- exp(-0.5*(lmeaic.ace.nb1-min(lmeaic.ace.nb1)))/sum(exp(-0.5*(lmeaic.ace.nb1-min(lmeaic.ace.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
lmewgt1_acetyl <- lmewgt_acetyl
lmewgt1_acetyl[lmewgt1_acetyl == lmewgt1_acetyl[1]] <- 0

#Make a dataframe containing positions and weights
lmeplot_data_ace <- data.frame(pos1,lmewgt1_acetyl) #Refer to plot_data_ace for pos1
names(lmeplot_data_ace) <- c("Position","Weight")
lmeplot_data_ace$Product <- ps_gene$Product[ps_gene$Position %in% lmeplot_data_ace$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(lmeplot_data_ace,"lme_plot_data_acetyl.csv")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
lmeplot_acetyl <- ggplot(data=lmeplot_data_ace, aes(x=lmeplot_data_ace$Position, y=(lmeplot_data_ace$Weight^0.005)^10, group=1, text=lmeplot_data_ace$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = lmeplot_data_ace$Position[which(lmeplot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(lmeplot_acetyl, tooltip = c("x","text"))


#Model for Alcohol
#fdat2$isolate <- as.factor(rownames(fdat2))
lmereg.alc.nb1 <- vector("list", ncol(fdat3)) 
lmeaic.alc.nb1 <- 1:(ncol(fdat3))
lmepred.alc.nb1<-matrix(nrow = nrow(fdat3), ncol = ncol(fdat3))

#fdat3[nrow(fdat3)+1,] <- sapply(fdat3, function(x) sum(as.numeric(x)) ) 

#fdat3 <- dplyr::filter(fdat3, fdat3[nrow(fdat3)+1,] < nrow(fdat3))

for (i in 1:ncol(fdat3)){
  if(nlevels(as.factor(fdat3[,i])) > 1){
    #Code to loop through the error
    tryCatch({lmereg.alc.nb1[[i]]<-glmer(as.factor(st_access1$alcohol) ~ as.factor(fdat3[,i]) + 
                                           (1|as.factor(rownames(st_access1))), family = binomial)}, error = function(err){})
    if(is.null(lmereg.alc.nb1[[i]]) == FALSE){
      lmeaic.alc.nb1[[i]]<-extractAIC(lmereg.alc.nb1[[i]])[2]
      #lmepred.alc.nb1[,i]<-predict(lmereg.alc.nb1[[i]])
    }
    else {
      lmeaic.alc.nb1[[i]]<-3000
      #lmepred.alc.nb1[,i]<-"NA"
    }
  }
  else{
    #AIC values cannot be calculated since there is only one level
    lmereg.alc.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    lmeaic.alc.nb1[i]<-1000
    #Predicted value will be same as existing value
    #lmepred.alc.nb1[,i]<-fdat3[1,i]
  }
}

#Compute weights for each position
lmewgt_alcohol <- exp(-0.5*(lmeaic.alc.nb1-min(lmeaic.alc.nb1)))/sum(exp(-0.5*(lmeaic.alc.nb1-min(lmeaic.alc.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
lmewgt1_alcohol <- lmewgt_alcohol
lmewgt1_alcohol[lmewgt1_alcohol == lmewgt1_alcohol[1]] <- 0

#Make a dataframe containing positions and weights
lmeplot_data_alc <- data.frame(pos1,lmewgt1_alcohol) #Refer to plot_data_ace for pos1
names(lmeplot_data_alc) <- c("Position","Weight")
lmeplot_data_alc$Product <- ps_gene$Product[ps_gene$Position %in% lmeplot_data_alc$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(lmeplot_data_alc,"lme_plot_data_alcohol.csv")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
lmeplot_alcohol <- ggplot(data=lmeplot_data_alc, aes(x=lmeplot_data_alc$Position, y=(lmeplot_data_alc$Weight^0.005)^10, group=1, text=lmeplot_data_alc$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = lmeplot_data_alc$Position[which(lmeplot_data_alc$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(lmeplot_alcohol, tooltip = c("x","text"))



