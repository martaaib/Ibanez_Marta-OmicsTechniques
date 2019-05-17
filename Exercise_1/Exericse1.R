# Exercise 1
#Data loading
samples_data <- read.delim("~/Desktop/Universitat/2_ANY/3_TRIMESTRE/Omics Techniques/Ibanez_Marta-OmicsTechniques/samples_data.txt", header=FALSE, comment.char="#")
class(samples_data)
#The data is a matrix and it has 8 different columns, which I assume are the different samples.
sampleNames
#Columwise
summary(samples_data)

#Plot data frame

boxplot(samples_data, las = 2)
#As far as I can see the values are pretty similar, the valyes are between 2 and 8