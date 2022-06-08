# Epidemiology-Data-wrangling-and-Data-linkage
Epidemiology Covid data: R demo for data wrangling and data linkage

###Author: Ling Wang, Ph.D.
The University of Chicago Pritzker School of Medicine

-   [Personal Web](https://www.lingwangmicrobiome.com/)
-   [GitHub](https://github.com/lwang18/RNAseq-demo-using-Midway2/)
-   [LinkedIn](https://www.linkedin.com/in/lingwangncsu/)
-   [Google Scholar ](https://scholar.google.com/citations?user=dBCVIigAAAAJ&hl=en&authuser=1)


## 1: Read file into R. 
```{r eval = FALSE}
Vaccine <- read.csv("Vaccine.csv", header=T)

```


#############################################################
## 2: Convert to uppercase. 

```{r eval = FALSE}
Vaccine[,-1] <- as.data.frame(sapply(Vaccine[,-1], toupper))


```

#############################################################
## 3: Remove all punctuation, spaces, hyphens, etc.

```{r eval = FALSE}
Vaccine[,c(3,4)] <- as.data.frame(sapply(Vaccine[,c(3,4)], function(x) gsub("[[:punct:]]", "", x)))

Vaccine[,c(3,4)] <- as.data.frame(sapply(Vaccine[,c(3,4)], function(x) gsub(" ", "", x)))

#I also want to check if there are NAs that need to be removed.
which(is.na(Vaccine)) # row 53 has a NA
Vaccine <- na.omit(Vaccine) # remove row 53


```


#############################################################
## 4: Data cleaning: unique ID. 

```{r eval = FALSE}
# First examine the object class
class(Vaccine$VaccineDate) ## Results: "character"
Vaccine$VaccineDate<-as.Date(Vaccine$VaccineDate, "%m/%d/%Y") # convert to objects of class “Date” 
class(Vaccine$VaccineDate) ## Results: "Date"


library(dplyr)
Vaccine_ID <- as.data.frame(Vaccine %>% 
                            group_by(FirstName,LastName) %>% 
                            filter(Create_Date==min(Create_Date)))

```

#############################################################
## 5: Order the data by `ID` and `VaccineDate`. Create a new variable “dose” which indicates the order of each dose for each person.
```{r eval = FALSE}
Vaccine_dose <- as.data.frame(Vaccine_ID %>% 
                              group_by(FirstName,LastName) %>%
                              mutate(dose = rank(VaccineDate)) %>%
                              arrange(ID,VaccineDate) )

dim(Vaccine_dose) ## Results: [1] 31  7


```


#############################################################
## 6: Exclude people with more than 2 doses. 

```{r eval = FALSE}
Vaccine_keep2 <- as.data.frame(Vaccine_dose %>% 
                               group_by(ID) %>% 
                               filter(n() <= 2))

```


#############################################################
## 7: Create a new variable `Dose_DIFF` which represents the different days between the first dose and second dose. Bonus if you also show how to solve this using your own function.

```{r eval = FALSE}
Dose_DIFF <- as.data.frame(Vaccine_keep2 %>% 
                           group_by(FirstName,LastName) %>%
                           mutate(Dose_DIFF
                                  =max(VaccineDate)-min(VaccineDate)))

```

#############################################################
## 8: If the second dose date is received less than 22 days compared to the first dose date. Please delete all records for this person. 

```{r eval = FALSE}
Vaccine_22days <- as.data.frame(Dose_DIFF %>% 
                                group_by(ID) %>% 
                                filter(!any(Dose_DIFF < 22)) )

```

#############################################################
## 9: Exclude those people if their `DateOfBirth` is the same as the `Create_Date`.

```{r eval = FALSE}
# Convert class of objects to Date
Vaccine_22days$DateOfBirth <- as.Date(Vaccine_22days$DateOfBirth, "%m/%d/%Y")
Vaccine_22days$Create_Date <- as.Date(Vaccine_22days$Create_Date, "%m/%d/%Y")

Vaccine_DateOfBirth <- 
  as.data.frame(Vaccine_22days %>% 
                group_by(ID) %>% 
                filter(!any(DateOfBirth == Create_Date)))

 
# Not run
# Vaccine_DateOfBirth2 <- 
# as.data.frame(Vaccine_DateOfBirth %>% 
#               group_by(ID) %>% 
#               filter(!any(DateOfBirth > VaccineDate)))

```

#############################################################
## 10: Exclude those people if their “DateOfBirth” is recorded as future date. Assuming today is 9/15/2021, please exclude those people if their `DateOfBirth` is larger than 9/15/2021.

```{r eval = FALSE}
Vaccine_future <- as.data.frame(Vaccine_DateOfBirth %>% 
                                group_by(ID) %>% 
                                filter(!any(DateOfBirth > "2021-09-15" )))

```


#############################################################
## 11: Transform the data from long format into wide format such that the information for a single person is in one row. Drop the `Dose` and the `Dose_DIFF` columns and create two new variables `First_dose_date` and `Second_dose_date`

```{r eval = FALSE}
Vaccine_wide <- as.data.frame(Vaccine_future %>% 
                              group_by(ID) %>%
                              mutate(First_dose_date=min(VaccineDate),
                                     Second_dose_date=max(VaccineDate)) %>% 
                              select(-Dose_DIFF, -VaccineDate, -dose) %>%
                              distinct(FirstName,LastName,DateOfBirth, .keep_all= TRUE) )

```


#############################################################
## 12: Merge data 
```{r eval = FALSE}
Covid <- read.csv("COVID.csv", header=T)

# Fix uppercases,punctuation, spaces, etc
Covid[,-1] <- as.data.frame(sapply(Covid[,-1], toupper))
Covid[,c(3,4)] <- as.data.frame(sapply(Covid[,c(3,4)], function(x) gsub("[[:punct:]]", "", x)))
Covid[,c(3,4)] <- as.data.frame(sapply(Covid[,c(3,4)], function(x) gsub(" ", "", x)))

# Fix data class
Covid$DOB <- as.Date(Covid$DOB, "%m/%d/%Y")

# Rename columns
Covid <- Covid %>%
  rename(FirstName = First_Name,
         LastName = Last_Name,
         DateOfBirth = DOB) # To be consistent with Vaccine_wide dataset

# Merge datasets
MergeData <- merge(Vaccine_wide, Covid, 
                   by.x = c("FirstName","LastName","DateOfBirth"), 
                   by.y = c("FirstName","LastName","DateOfBirth"), all=F) 
# If we want to keep those mismatches, change all=T 


dim(MergeData) ## Results: [1] 6 9
## Looks like #s of rows are reduced from 7 to 6 (1 more row for GIROLAMO FRACASTORO as he's likely been tested positive twice, and ABIGAL SALK & TONY POWERS are dropped).

# Lastly, output data in csv file
write.csv(MergeData,"LingWang_MergeData.csv", row.names = FALSE)

```


#############################################################
#############################################################
## Data Linkage

```{r eval = FALSE}
# The goal is to identify those with the same "LastName" & "DateOfBirth" since first name can be confused by nicknames. So unlike in Q12 where "FirstName" is used in by.x and by.y to merge, here I exclude "FirstName". 
MergeData2 <- merge(Vaccine_wide, Covid, 
                    by.x = c("LastName","DateOfBirth"), 
                    by.y = c("LastName","DateOfBirth"), 
                    all=F, sort = F) 


# Now instead of 6 rows (or 6 people) that have the exact same "FirstName","LastName" and "DateOfBirth", I get 8 rows (or 8 people) that have exact "LastName" and "DateOfBirth". This means that there are 2 people that probably don't have the same "FirstName" between Vaccine_wide & Covid datasets. 
```

```{r eval = FALSE}
# Then I want to pull out the data of these two specific people.
WhoHaveNicknames <- 
  as.data.frame(MergeData2 %>% 
                  group_by(ID) %>% 
                  filter(!any(FirstName.x == FirstName.y)))

```

```{r eval = FALSE}
# Option 1
# I can  visually inspect their first names (ABIGAL vs ABBY, TONY vs TONNY) in columns FirstName.x and FirstName.y to decide whether they are the same first name. 


# Option 2
# I can add a column "match" to visually indicate whether the first names are matched between Vaccine_wide & Covid datasets for each person. 
MergeData2$match <- lengths(Map(intersect, 
                                MergeData2$FirstName.x,
                                MergeData2$FirstName.y)) > 0


```

```{r eval = FALSE}
# Option 3
# Use a package 'stringsim' which returns a vector with similarities using values between 0 and 1 where 1 corresponds to perfect similarity (distance 0) and 0 to complete dissimilarity. 

library(stringdist)

# Calculate the similarity using the default method of optimal string alignment (osa)
# Playing around with different names
stringsim("ca", "abc") ## Results: 0
stringsim("ca", "ca") ## Results: 1
stringsim("ca", "caa") ## Results: 0.6666667
stringsim("ca", "cva") ## Results: 0.6666667
stringsim("Tony", "Tonny") ## Results: 0.8
stringsim("ABIGAL", "ABBY") ## Results: 0.3333333
stringsim("Joann", "Joanna") ## Results: 0.8333333
stringsim("Joey", "Joseph") ## Results: 0.5

# Calculate the similarity using the Jaro-Winkler method (from my research it seems to be the most popular method in this field). 
stringsim('Tony','Tonny',method='jw') ## Results: 0.9333333
stringsim('ABIGAL','ABBY',method='jw') ## Results: 0.6111111
stringsim('Joann','Joanna',method='jw') ## Results: 0.9444444
stringsim('Joey','Joseph',method='jw') ## Results: 0.75
# Looks to me a threshold of 0.6 seems to be a good start. 
# If the similarity value is over 0.6 for a row, this means that this person's first names between two datasets are "similar enough" and its row data can be kept (likely the same person). 


# Next I want to add these results to new columns
MergeData2 <- MergeData2 %>% 
    mutate(Similarity = 
           stringsim(FirstName.x, FirstName.y, method ='jw'))

MergeData2$Keep <- ifelse(MergeData2$Similarity>0.6 , "YES", "NO")


# Lastly: rearrange output
MergeData2 <- MergeData2 %>%
  rename(FirstName = FirstName.x)  %>%
  select(ID, FirstName, LastName, DateOfBirth,
         Create_Date, First_dose_date, 
         Second_dose_date,COVID_ID, COVID_Date, Keep) %>%
  filter(any(Keep == "YES" )) %>%
  select(-Keep)

```

##############################################################
## More attempts:
## From researching the field data linkage, I also tried 'reclin'

```{r eval = FALSE}
library(reclin)
p <- pair_blocking(Vaccine_wide, Covid, "LastName", large = FALSE)
print(p)

## Results:
## Simple blocking
## Blocking variable(s): LastName
## First data set:  7 records
## Second data set: 16 records
## Total number of pairs: 8 pairs

## Showing all pairs:
##   x y
## 1 1 1 
## 2 2 3
## 3 3 4
## 4 4 5
## 5 5 6
## 6 6 7 
## 7 6 8 
## 8 7 9
# Record 1 from Vaccine_wide is compared to record 1 from Covid.
# Record 6 from Vaccine_wide is compared to records 7 and 8 from Covid.

p <- compare_pairs(p, by = c("FirstName", "LastName", "DateOfBirth"))
print(p)

## Results:
##   x y FirstName LastName DateOfBirth
## 1 1 1      TRUE     TRUE        TRUE
## 2 2 3      TRUE     TRUE        TRUE
## 3 3 4      TRUE     TRUE        TRUE
## 4 4 5      TRUE     TRUE        TRUE
## 5 5 6     FALSE     TRUE        TRUE 
## 6 6 7      TRUE     TRUE        TRUE
## 7 6 8      TRUE     TRUE        TRUE
## 8 7 9     FALSE     TRUE        TRUE 

p <- compare_pairs(p, by = c("FirstName", "LastName", "DateOfBirth"),
                   default_comparator = jaro_winkler(0.9), overwrite = TRUE)
print(p)

## Results:
##     x y FirstName LastName DateOfBirth
##   1 1 1 1.0000000        1           1
##   2 2 3 1.0000000        1           1
##   3 3 4 1.0000000        1           1
##   4 4 5 1.0000000        1           1
##   5 5 6 0.6111111        1           1
##   6 6 7 1.0000000        1           1
##   7 6 8 1.0000000        1           1
##   8 7 9 0.9333333        1           1

```

##############################################################
## More attempts:
## From researching data linkage, I also tried 'fastLink'

```{r eval = FALSE}
library(fastLink) # Fast Probabilistic Record Linkage

dfA = Vaccine_wide
dfB = Covid

matches.out <- fastLink(
  dfA = dfA, dfB = dfB, 
  varnames = c("FirstName", "LastName", "DateOfBirth"),
  stringdist.match = c("FirstName", "LastName", "DateOfBirth"),
  partial.match = c("FirstName", "LastName")  )

matched_dfs <- getMatches(
  dfA = dfA, dfB = dfB, 
  fl.out = matches.out, threshold.match = 0.85)

matched_dfs

summary(matches.out)

## Results:
##                95%  85%  75%   Exact
## 1 Match Count    7    7    7       6
## 2  Match Rate 100% 100% 100% 85.714%
## 3         FDR   0%   0%   0%        
## 4         FNR   0%   0%   0% 

```

## I think data linkage is very interesting and definitely something worth studying post covid era as it is not surprising to receive data with names misspelled or nicknames (or maybe even pseudo names if people are worried about releasing private medical data).
