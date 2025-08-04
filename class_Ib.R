#create subfolders
dir.create("raw_data") #storing raw data
dir.create("clean_data") #storing clean data
dir.create("script") #saving scripts
dir.create("results") #saving analysis output
dir.create("plots") #visualization
#import data
data <- read.csv(file.choose())
#view data
View(data)
#structure of dataset
str(data)
#convert character into factor
data$diagnosis_fac <- as.factor(data$diagnosis)
str(data)
#convert factor(categorical data) into numeric factor
data$diagnosis_fac <- ifelse(data$diagnosis_fac == "Cancer", 1, 0)
class(data$diagnosis_num)
#gender to factor
data$gender_fac <- as.factor(data$gender)
str(data)
#convert factor(categorical data) into numeric factor
data$gender_num <- ifelse(data$gender_fac == "Female", 1, 0)
class(data$gender_num)
#convert from numeric factor to factor 
data$gender_num <- as.factor(data$gender_num)
class(data$gender_num)
#convert numeric to integer
data$bmi
class(data$bmi)
as.integer(data$bmi)
print <- as.integer(data$bmi)
data$bmi_fac <- as.integer(data$bmi)
str(data)
#convert factor to numeric factor
my_data= smoker_fac[,-11]
data$smoker_fac <- 
str(data)
data <- data[,-11]
#smoking status as a binary factor
smoking_status <- c(0,1)
  labels= c("No", "Yes")


