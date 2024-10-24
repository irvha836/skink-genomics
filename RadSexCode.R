# Installation from GitHub requires devtools
#install.packages("devtools")
#devtools::install_github("SexGenomicsToolkit/sgtr")

#code for plotting heat map
library(sgtr)

input_file= "distribution.tsv"

radsex_distrib(input_file=input_file, groups =c("M", "F"))



#code for plotting bootstrap distributions
inputfile= "500radsexbootsrap.txt"

hist(500radsexbootstrap.txt")
# Load the data
bootstrap_data <- read.table(500radsexbootstrap.txt, header = FALSE)

# Extract the column with bootstrap values
bootstrap_values <- bootstrap_data$V1

# Plot the histogram
hist(bootstrap_values, 
     main = "Histogram of Bootstrap Results", 
     xlab = "Bootstrap Values", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black")

# Optional: Add a line for the mean
abline(v = mean(bootstrap_values), col = "red", lwd = 2, lty = 2)


# Load the data
bootstrap_data <- read.table("500radsexbootstrap.txt", header = FALSE)

# Assign column names for clarity
colnames(bootstrap_data) <- c("Female_Sum", "Male_Sum")

# Plot histogram for Female_Sum
hist(bootstrap_data$Female_Sum, 
     main = "Histogram of Bootstrap Female Sum", 
     xlab = "Sum of Markers (Females)", 
     col = "lightblue", 
     border = "black")
     
     abline(v = 30, col = "red", lwd = 2, lty = 2)

# Plot histogram for Male_Sum
hist(bootstrap_data$Male_Sum, 
     main = "Histogram of Bootstrap Male Sum", 
     xlab = "Sum of Markers (Males)", 
     col = "lightcoral", 
     border = "black")
     abline(v = 30, col = "red", lwd = 2, lty = 2)



