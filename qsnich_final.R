library(readxl)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)
library(dplyr)
library(stringr)
library(ggplot2)
library(rstatix)
library(ggridges)


# gps coordinates and ELISA results
elisa =
    "00-RawData/QSNICHRetrospectiveDataGeocoding.csv" %>%
    read_csv %>%
    rename_all(function(x) gsub(' +','', x)) %>%
    rename_all(function(x) gsub('\r|\n|\\)','', x)) %>%
    rename_all(function(x) gsub('\\(','_', x)) %>%
    # correct longitude of one person who was likely mis-typed
    mutate(HomeLongitude = ifelse(HomeLongitude > 50
        , HomeLongitude
        , gsub('^1', '10', as.character(HomeLongitude)) %>% as.numeric
    )) %>%
    mutate_at(vars(matches('[Dd]ate$')), as.Date, format = "%d/%m/%Y") %>%
    # make unavailable titer data NA
    mutate_at(vars(matches('Ig[GM]')), as.numeric)

# dates
f1 = "00-RawData/Requested information QSNICH and KPPH to Henrik_25JULY23.xlsx"
f2 = "00-RawData/Requested information QSNICH to Henrik_04AUG23.xlsx"

# QSNICH and elisa data
qs1 = f1 %>%
    read_excel(sheet = grep('QSNICH', excel_sheets(f1), value = T)) %>%
    rename_all(function(x) gsub(' +','', x)) %>%
    rename_all(function(x) gsub('-','', x)) %>%
    mutate_if(is.POSIXct, as.Date)

qs2 =
    f2 %>%
    read_excel() %>%
    rename_all(function(x) gsub(' +','', x)) %>%
    rename_all(function(x) gsub('-','', x)) %>%
    mutate_if(is.POSIXct, as.Date)

elisa %>%
    select(matches('Ig[GM]')) %>%
    lapply(class)


# function for data cleaning
cleanData <- function(df, filterList) {

    # converting all dates from Buddhist Era to Common Era (e.g. SubjectNo = SV20-0412)
    df <- mutate(df, Dateill1st = case_when(
        Dateill1st > as.Date("2023-01-01") ~ Dateill1st - 365*543, 
        TRUE    ~ Dateill1st
        ))
    
    df <- mutate(df, AdmitDate = case_when(
        AdmitDate > as.Date("2023-01-01") ~ AdmitDate - 365*543, 
        TRUE    ~ AdmitDate
        ))

    # addition of Days, Age and PCR columns
    df <- df %>%
        mutate(
            Days = as.numeric(dateofcollection - Dateill1st, unit = "days") %>% as.integer
            , Age = Ageyr + Agemo/12
            , PCR = ifelse(DENVPCR != "NEG", "PCR+", "PCR-")
        )

    # simplification of clinical presentation data into DF, DHF and NA
    df$ClinicalDiagnosis[df$ClinicalDiagnosis %in% c("Dengue Fever", "Viral Infection DF")] <- "DF"
    df$ClinicalDiagnosis[substr(df$ClinicalDiagnosis, 1, 3) %in% c("DHF", "DSS")] <- "DHF"
    df$ClinicalDiagnosis[substr(df$ClinicalDiagnosis, 1, 21) %in% c("Dengue Shock Syndrome", "Dengue shock syndrome", "Dengue Hemorrhagic Fe")] <- "DHF"
    df$ClinicalDiagnosis[!(df$ClinicalDiagnosis %in% c("DF", "DHF"))] <- NA

    dfFinal <- df

    # filtering of data
    for (filterChoice in filterList) {
        if (filterChoice %in% c(1,2)) {
        
            # filtering of dates
            # removes dates where the time period is too great for best acute Ig measurements 
            # and rows where dateofcollection is earlier than Dayill1st
            dfFinal <- dfFinal %>%
                filter(dateofcollection >= Dateill1st
                    , abs(dateofcollection - Dateill1st) < 10
                    , abs(AdmitDate - Dateill1st) < 10
                    , Days > 0
                )
        } 

        # filtering of PCR=NA
        if (filterChoice %in% c(1,3)) {
            dfFinal <- dfFinal %>%
                filter (PCR %in% c("PCR+", "PCR-"))
        }

        # filtering of ages < 1.5
        if (filterChoice %in% c(1,4)) {
            dfFinal <- dfFinal %>%
                filter (Age >= 1.5)
        }
    }
    # orders by dateill1st
    dfFinal <- dfFinal[order(dfFinal$Dateill1st),]
        
    return(dfFinal)
}


# function removing rows with the same AFRIMSID
cleanElisa <- function(dfElisa) {
    duplicateRows = (dfElisa[duplicated(dfElisa$AFRIMSID),])

    IDtoremove = list()

    # looping through every duplicated AFRIMSID
    for (id in unique(duplicateRows$AFRIMSID)) {
        # creates a temporary dataframe with rows with duplicate AFRIMSID
        dfDuplicateID = dfElisa %>%
            select(AFRIMSID, CaseNo) %>%
            filter(AFRIMSID == id) %>%
            mutate(editedCaseNo = str_sub(CaseNo, 2, 5)) %>%
            arrange(desc(CaseNo))
        # storing every row with a duplicate AFRIMSID except for the row added last
        for (eachDuplicate in 1:length(dfDuplicateID)) {
            if (eachDuplicate != 1) {
                IDtoremove = c(IDtoremove, dfDuplicateID[eachDuplicate,"CaseNo"])
            }
        }
    }

    # removing duplicate rows
    dfElisaFinal = dfElisa %>%
        filter(!(CaseNo %in% IDtoremove))
    
    # renaming column with clinical presentations and adding columns with IgG/IgM log ratios
    dfElisaFinal = dfElisaFinal %>%
        rename(FinalDiagnosis = `FinalDiagnosis_Notdengue,DF,DHF,DHF1,DHF2,DHF3,DHF4,DSS`) %>%
        mutate(logRatio = log2(AcuteIgG_1stblooddraw + 1) - log2(AcuteIgM_1stblooddraw + 1)
        , logRatioConv = log2(ConvalescentIgG_2ndblooddraw + 1) - log2(ConvalescentIgM_2ndblooddraw + 1))

    # simplification of clinical presentation data into DF, DHF and NA
    dfElisaFinal$FinalDiagnosis[dfElisaFinal$FinalDiagnosis %in% c("DF", "Influenza with DF", "Df")] <- "DF"
    dfElisaFinal$FinalDiagnosis[substr(dfElisaFinal$FinalDiagnosis, 1, 3) %in% c("DHF", "DSS")] <- "DHF"
    dfElisaFinal$FinalDiagnosis[dfElisaFinal$FinalDiagnosis == "Dengue Hemorrhagic Fever without shock"] <- "DHF"
    dfElisaFinal$FinalDiagnosis[!(dfElisaFinal$FinalDiagnosis %in% c("DF", "DHF"))] <- NA

    return(dfElisaFinal)
}


# function merging dataframe with elisa dataframe
mergeElisa <- function(df, dfElisa, filterChoice) {
    dfFinal <- df %>%
        inner_join(dfElisa
        , by = c('SubjectNo' = 'AFRIMSID')
    )

    # removal of conflicting data entries
    if (filterChoice == 1) {
        dfFinal <- dfFinal %>%
            filter( #serointerpretation == Interpretation | !is.na(serointerpretation) | !is.na(Interpretation)
                !is.na(Sex),
                !is.na(Gender),
                !is.na(logRatio),
                substr(Gender, 1, 1) == Sex, 
                AdmitDate == AdmissionDate,
                , ClinicalDiagnosis == FinalDiagnosis
                , PCR == "PCR+")
    }
    return (dfFinal)
}


# merging of both datasets
qs <- rbind(qs1, qs2)

# reformatting of merged dataset without filtering
qsFormatted <- cleanData(qs,0)
qs1Formatted <- cleanData(qs1,0)
qs2Formatted <- cleanData(qs2,0)

# cleaning of merged and unmerged datasets
qsCleaned <- cleanData(qs,1)
qs1Cleaned <- cleanData(qs1,1)
qs2Cleaned <- cleanData(qs2,1)

cleanDatesOnly <- cleanData(qs,2)
cleanPcrOnly <- cleanData(qs,3)
cleanAgeOnly <- cleanData(qs,4)

cleanDatesPcr <- cleanData(qs,list(2,3))

# cleaning of elisa data
datElisa = cleanElisa(elisa)

# merging of datasets with elisa dataset and cleaning 
dat <- mergeElisa(qsCleaned, datElisa, 1)
#dat <- mergeElisa(qsCleaned, datElisa, 2)
dat1 <- mergeElisa(qs1Cleaned, datElisa, 1)
dat2 <- mergeElisa(qs2Cleaned, datElisa, 1)

# merging of datasets with elisa dataset without cleaning
datUncleaned <- mergeElisa(qsCleaned, datElisa, 0)

# merging of datasets with elisa dataset without cleaning elisa or qs
datUncleanedAll <- mergeElisa(qsFormatted, datElisa, 0)

# validation of cleaning of dataset
View(dat)


# Descriptions of data

# earliest and latest study years
print(min(qsFormatted$Studyyear))
print(max(qsFormatted$Studyyear))

print(min(datOld$Studyyear))
print(max(datOld$Studyyear))

print(min(dat$Studyyear))
print(max(dat$Studyyear))

# earliest and latest dates (roughly - for months and years)
print(min(na.omit(qsFormatted$AdmitDate)))
print(min(na.omit(qsFormatted$Dateill1st)))
print(min(na.omit(qsFormatted$dateofcollection)))
print(max(na.omit(qsFormatted$AdmitDate)))
print(max(na.omit(qsFormatted$Dateill1st)))
print(max(na.omit(qsFormatted$dateofcollection)))

print(min(na.omit(dat$AdmitDate)))
print(min(na.omit(dat$Dateill1st)))
print(min(na.omit(dat$dateofcollection)))
print(max(na.omit(dat$AdmitDate)))
print(max(na.omit(dat$Dateill1st)))
print(max(na.omit(dat$dateofcollection)))

print(min(elisa$AdmissionDate))
print(max(elisa$AdmissionDate))

# number of unique individuals
print(length(unique(qsFormatted$SubjectNo)))
print(length(unique(dat$SubjectNo)))

# number of individuals with acute and convalescent data
print(nrow(dat))
print(nrow(dat[!(is.na(dat$ConvalescentIgG_2ndblooddraw == "PCR-")),]))

# number and percentage of males and females
genderDemo <- function(dat) {
    output = c()
    maleTotal = nrow(dat[dat$Sex == "M",])
    femaleTotal = nrow(dat[dat$Sex == "F",])
    output = c(output, "Male", maleTotal, maleTotal * 100 / nrow(dat), "Female", femaleTotal, femaleTotal * 100 / nrow(dat))
    datConv <- dat[!(is.na(dat$ConvalescentIgM_2ndblooddraw)),]
    maleTotalConv = nrow(datConv[datConv$Sex == "M",])
    femaleTotalConv = nrow(datConv[datConv$Sex == "F",])
    output = c(output, "Male Conv", maleTotalConv, maleTotalConv * 100 / nrow(datConv), "Female Conv", femaleTotalConv, femaleTotalConv * 100 / nrow(datConv))
}
genderDemoResults = genderDemo(dat)
print(genderDemoResults[1:3])
print(genderDemoResults[4:6])
print(genderDemoResults[7:9])
print(genderDemoResults[10:12])

# number and percentage of patients with DF or DHF
severityDemo <- function(dat) {
    output = c()
    dfTotal = nrow(dat[dat$ClinicalDiagnosis == "DF",])
    dhfTotal = nrow(dat[dat$ClinicalDiagnosis == "DHF",])
    output = c(output, "DF", dfTotal, dfTotal * 100 / nrow(dat), "DHF", dhfTotal, dhfTotal * 100 / nrow(dat))
    datConv <- dat[!(is.na(dat$ConvalescentIgM_2ndblooddraw)),]
    dfTotalConv = nrow(datConv[datConv$ClinicalDiagnosis == "DF",])
    dhfTotalConv = nrow(datConv[datConv$ClinicalDiagnosis == "DHF",])
    output = c(output, "DF Conv", dfTotalConv, dfTotalConv * 100 / nrow(datConv), "DHF Conv", dhfTotalConv, dhfTotalConv * 100 / nrow(datConv))
    
}
severityDemoResults = severityDemo(dat)
print(severityDemoResults[1:3])
print(severityDemoResults[4:6])
print(severityDemoResults[7:9])
print(severityDemoResults[10:12])

# mean and sd of patient ages
mean(dat$Age)
sd(dat$Age)
datConv <- dat[!(is.na(dat$ConvalescentIgM_2ndblooddraw)),]
mean(datConv$Age)
sd(datConv$Age)

# sample sizes for days 6-9
print(nrow(dat[dat$Days == 6,]))
print(nrow(dat[dat$Days == 7,]))
print(nrow(dat[dat$Days == 8,]))
print(nrow(dat[dat$Days == 9,]))



# Descriptive plots of data


# histograms with disease severity by age with 1-year age bins
dat %>%
    ggplot(aes(x = Age, fill = ClinicalDiagnosis), group = ClinicalDiagnosis) +    # fill=AgeGroup
    #geom_histogram()
    geom_histogram(bins = max(dat$Ageyr))

# box plot of Age x Days, grouped by PCR and disease severity
dat %>%
    ggplot(aes(x = Days, y = Age, group = paste(Days,ClinicalDiagnosis), fill = ClinicalDiagnosis)) +
    geom_boxplot()

# acute/conv IgM/IgG x Days, comparing differences in disease severity over time
ggarrange(
    dat %>%
        ggplot(aes(x = Days, y = AcuteIgM_1stblooddraw, group = paste(Days,ClinicalDiagnosis), fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute IgM titres") +
        geom_boxplot(outlier.size = 0.5) +
        guides(fill=guide_legend(title="Disease severity")) +
        theme(legend.position = "none") 
    , dat %>%
        ggplot(aes(x = ClinicalDiagnosis, y = ConvalescentIgM_2ndblooddraw, fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        ylab("Convalescent IgM titres") +
        xlab("Disease severity") +
        geom_boxplot(outlier.size = 0.5) +
        theme(legend.position = "none")
    , dat %>%
        ggplot(aes(x = Days, y = AcuteIgG_1stblooddraw, group = paste(Days,ClinicalDiagnosis), fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute IgG titres") +
        geom_boxplot(outlier.size = 0.5) +
        theme(legend.position = "none")
    , dat %>%
        ggplot(aes(x = ClinicalDiagnosis, y = ConvalescentIgG_2ndblooddraw, fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        ylab("Convalescent IgG titres") +
        xlab("Disease severity") +
        geom_boxplot(outlier.size = 0.5) +
        theme(legend.position = "none")
    , dat %>%
        ggplot(aes(x = Days, y = logRatio, group = paste(Days,ClinicalDiagnosis), fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute log(IgG:IgM)") +
        geom_boxplot(outlier.size = 0.5) +
        theme(legend.position = "none")
    , dat %>%
        ggplot(aes(x = ClinicalDiagnosis, y = logRatioConv, fill = ClinicalDiagnosis)) +  
        scale_fill_brewer(palette="Set2") +
        ylab("Convalescent log(IgG:IgM)") +
        xlab("Disease severity") +
        geom_boxplot(outlier.size = 0.5) +
        theme(legend.position = "none")
    , ncol = 2
    , nrow = 3
    , widths = c(1,0.35)
    , common.legend = TRUE
    , legend = "right"
)

# box plots of Ig titres x Days, grouped into three age categories
ggarrange(
    dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = AcuteIgM_1stblooddraw, x = Days, group = paste(Days, as.integer(AgeGroup)), fill = AgeGroup)) +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute IgM titres") +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_discrete(name = "Age groups (years)", labels = c("1.5-6", "6-12", ">12")) +
        #guides(fill=guide_legend(title="Age groups")) +
        theme(legend.position = "none")
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = ConvalescentIgM_2ndblooddraw, x = AgeGroup, fill = AgeGroup)) +
        ylab("Convalescent IgM titres") +
        scale_x_discrete(name = "Age groups (years)", labels = c("1.5-6", "6-12", ">12")) +
        geom_boxplot(outlier.size=0.5) +
        theme(legend.position = "none")
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = AcuteIgG_1stblooddraw, x = Days, group = paste(Days, as.integer(AgeGroup)), fill = AgeGroup)) +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute IgG titres") +
        geom_boxplot(outlier.size=0.5) +
        theme(legend.position = "none")
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = ConvalescentIgG_2ndblooddraw, x = AgeGroup, fill = AgeGroup)) +
        ylab("Convalescent IgG titres") +
        scale_x_discrete(name = "Age groups (years)", labels = c("1.5-6", "6-12", ">12")) +
        geom_boxplot(outlier.size=0.5) +
        theme(legend.position = "none")
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = logRatio, x = Days, group = paste(Days, as.integer(AgeGroup)), fill = AgeGroup)) +
        scale_x_continuous(breaks = seq(0, 10, by = 1)) +
        ylab("Acute log(IgG:IgM)") +
        geom_boxplot(outlier.size=0.5) +
        theme(legend.position = "none")
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(y = logRatioConv, x = AgeGroup, fill = AgeGroup)) +
        ylab("Convalescent log(IgG:IgM)") +
        scale_x_discrete(name = "Age groups (years)", labels = c("1.5-6", "6-12", ">12")) +
        geom_boxplot(outlier.size=0.5) +
        theme(legend.position = "none")
    , ncol = 2
    , nrow = 3
    , widths = c(0.9,0.25)
    , common.legend = TRUE
    , legend = "right"
)



    #   Infer latent class of data points using EM
    #   (EM code modified from: https://rpubs.com/H_Zhu/246450)
    #   .................................

source("C:/Users/kaavy/VSCodeProjects/dengueProject/dengueAgAbKinetics/Scripts/functions/em_Kaavya.R")

# Models fitted with acute antibody sample data

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM) data
fit_all3 =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5)

# fitting the model with only log(IgG) data
fit_IgG_only =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputIgM = FALSE, inputRatio = FALSE)

# fitting the model with log(IgG) and log(IgM), with only IgG kinetics visualised in plots
fit_Igs_IgG =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputRatio = FALSE, view = 2)

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM), with only IgG kinetics visualised in plots
fit_all3_IgG = 
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, view = 2)

# fitting the model with only log(IgM)
fit_IgM_only =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputIgG = FALSE, inputRatio = FALSE)

# fitting the model with log(IgG) and log(IgM), with only IgM kinetics visualised in plots
fit_Igs_IgM =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputRatio = FALSE, view = 3)

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM), with only IgM kinetics visualised in plots
fit_all3_IgM = 
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, view = 3)

# fitting the model with only log(IgG:IgM), with only IgM kinetics visualised in plots
fit_ratio_only =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputIgM = FALSE, inputIgG = FALSE)

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM), with only IgG:IgM kinetics visualised in plots
fit_all3_ratio = 
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, view = 4)

# fitting the model with log(IgG) and log(IgM)
fit_Igs =
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, inputRatio = FALSE, view = 5)

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM), with only IgG and IgM kinetics visualised in plots
fit_all3_Igs = 
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init) %>%
    fitResponse(nPoly=5, view = 5)


# improvement of likelihood over the iterations
plot(fit$Q, ylab = 'LL', xlab = 'Iteration', type = 'l')

# fit acts as a placeholder for any fitted model of interest
fit = fit_all3

# colours used to plot antibody kinetics
colourList = c()
colourList = c(colourList,"#00B934","#e9b141") # for IgG
colourList = c(colourList, "#5D99FF","#F561E2") # for IgM
colourList = c(colourList, "#F8746B","#00BEC3") # for IgG:IgM

# fitted average lines of the kinetics
ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
    scale_color_manual(values=colourList) +
    scale_fill_manual(values=colourList) +
    geom_hline(yintercept = 0, linetype = 3)+
    # initial splines
    geom_line(data = fit$init, linetype = 2, linewidth = 0.6)+
    # final fitted splines
    geom_ribbon(data = fit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
    geom_line(data = fit$final, linewidth = 0.5)+
    ylab("log(Ig) titres")+
    scale_x_continuous(breaks = seq(0, 10, by = 1))


# bimodal histogram of P(Primary)
dat %>%
    mutate(pPrimary = fit_all3$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")


# plots of the model, isolated predictions of antibody kinetics and bimodal histogram of P(Primary)
ggarrange(
    ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=colourList) +
        scale_fill_manual(values=colourList) +
        geom_hline(yintercept = 0, linetype = 3)+
        # final fitted splines
        geom_ribbon(data = fit_all3$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_all3$final, linewidth = 0.5)+
        ylab("log(Ig) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , NULL
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#00B934","#e9b141")) +
        scale_fill_manual(values=c("#00B934","#e9b141")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_all3_IgG$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_all3_IgG$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_all3_IgG$final, linewidth = 0.5)+
        ylab("log(IgG) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#5D99FF","#F561E2")) +
        scale_fill_manual(values=c("#5D99FF","#F561E2")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_all3_IgM$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_all3_IgM$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_all3_IgM$final, linewidth = 0.5)+
        ylab("log(IgM) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))+
        theme(legend.position = "none")
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#F8746B","#00BEC3")) +
        scale_fill_manual(values=c("#F8746B","#00BEC3")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_all3_ratio$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_all3_ratio$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_all3_ratio$final, linewidth = 0.5)+
        ylab("log(IgG:IgM) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))+
        theme(legend.position = "none")
    , dat %>%
    mutate(pPrimary = fit_all3$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")
    , ncol = 2
    , nrow = 3
    , common.legend = TRUE
    , legend = "right"
)


# plots of the models fitted with subsets of the available antibody response variables
ggarrange(
    ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#00B934","#e9b141")) +
        scale_fill_manual(values=c("#00B934","#e9b141")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_IgG_only$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_IgG_only$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_IgG_only$final, linewidth = 0.5)+
        ylab("log(IgG) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , dat %>%
    mutate(pPrimary = fit_IgG_only$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#5D99FF","#F561E2")) +
        scale_fill_manual(values=c("#5D99FF","#F561E2")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_IgM_only$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_IgM_only$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_IgM_only$final, linewidth = 0.5)+
        ylab("log(IgM) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , dat %>%
    mutate(pPrimary = fit_IgM_only$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#F8746B","#00BEC3")) +
        scale_fill_manual(values=c("#F8746B","#00BEC3")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_ratio_only$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_ratio_only$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_ratio_only$final, linewidth = 0.5)+
        ylab("log(IgG:IgM) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , dat %>%
    mutate(pPrimary = fit_ratio_only$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")
    , ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
        scale_color_manual(values=c("#00B934","#e9b141","#5D99FF","#F561E2")) +
        scale_fill_manual(values=c("#00B934","#e9b141","#5D99FF","#F561E2")) +
        geom_hline(yintercept = 0, linetype = 3)+
        # initial splines
        geom_line(data = fit_Igs$init, linetype = 2, linewidth = 0.6)+
        # final fitted splines
        geom_ribbon(data = fit_Igs$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
        geom_line(data = fit_Igs$final, linewidth = 0.5)+
        ylab("log(Ig) titres")+
        scale_x_continuous(breaks = seq(0, 10, by = 1))
    , dat %>%
    mutate(pPrimary = fit_Igs$pPrimary) %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("Density")+
    xlab("P(Primary)")
    , ncol = 2
    , nrow = 4
    , widths = c(1,0.5)
)


# plots of p(Primary) by disease severity and age
ggarrange(
    dat %>%
        ggplot(aes(x = ClinicalDiagnosis, y = fit_all3$pPrimary, fill = ClinicalDiagnosis)) +
        scale_fill_brewer(palette="Set2") +
        geom_violin(adjust=1.5) +
        geom_boxplot(width = 0.1, fill = NA) + #+#adjust=1.5, 
        ylab("P(Primary)")+
        xlab("Disease Severity")+
        guides(fill=guide_legend(title="Disease Severity"))
    , dat %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(x = AgeGroup, y = fit_all3$pPrimary, fill = AgeGroup)) +
        geom_violin(adjust=1.5, ) +
        geom_boxplot(width = 0.1) + #+#adjust=1.5, 
        ylab("P(Primary)")+
        scale_x_discrete(name = "Age (years)", labels = c("1.5-6", "6-12", ">12")) +
        scale_fill_discrete(name = "Age (years)", labels = c("1.5-6", "6-12", ">12")) 
    , ncol = 1
    , nrow = 2
)


# calculation of mean absolute deviation per person from pPrimary = 0.5 by day
dat %>% 
    mutate(pPrimaryMAD_allmodel = fit_all3$P1,
    pPrimaryMAD_IgGmodel = fit_IgG_only$P1,
    pPrimaryMAD_IgMmodel = fit_IgM_only$P1,
    pPrimaryMAD_ratiomodel = fit_ratio_only$P1,
    pPrimaryMAD_Igsmodel = fit_Igs$P1) %>%
    group_by(Days) %>%
    summarise(mad_all = sum(abs(0.5 - pPrimaryMAD_allmodel)) / n(),
    mad_IgG = sum(abs(0.5 - pPrimaryMAD_IgGmodel)) / n(),
    mad_IgM = sum(abs(0.5 - pPrimaryMAD_IgMmodel)) / n(),
    mad_ratio = sum(abs(0.5 - pPrimaryMAD_ratiomodel)) / n(),
    mad_Igs = sum(abs(0.5 - pPrimaryMAD_Igsmodel)) / n()
    ) %>%
    ggplot() +
    expand_limits(y=0) +
    geom_line(aes(x = Days, y = mad_all, color = "all"), linewidth = 0.6) +  #colourListMAD[1]
    geom_line(aes(x = Days, y = mad_IgG, color = "IgG"), linewidth = 0.6) + 
    geom_line(aes(x = Days, y = mad_IgM, color = "IgM"), linewidth = 0.6) + 
    geom_line(aes(x = Days, y = mad_ratio, color = "ratio"), linewidth = 0.6) + 
    geom_line(aes(x = Days, y = mad_Igs, color = "Igs"), linewidth = 0.6) + 
    scale_x_continuous(breaks = seq(0, 9, by = 1)) +
    ylab("MAD") +
    scale_color_manual(name='Models (labelled using \nfitted antibody variables)',
                     labels = c("IgG, IgM and IgG:IgM", "IgG", "IgM", "IgG and IgM", "IgG:IgM"),
                     values=c("all" = "#000000", "IgG" = "#00B934", "IgM" = "#5D99FF", "ratio" = "#F8746B", "Igs" = "#FF9922"))#+



# sensitivity analysis of final model for IgG by fitting polynomials 3-7
fit_all3_sensitivity = 
    dat %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1)
    ) %>%
    select(SubjectNo, Age, Days, logIgM, logIgG, logRatio, pPrimary_init)

# creation of a base sensitivity plot
sensitivityPlot = ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
geom_hline(yintercept = 0, linetype = 3)+
scale_x_continuous(breaks = seq(0, 14, by = 1))

# list of colours for each degree of polynomial
colorList = c("#000000", "#FF0000","#007700", "#0000FF","777700")

# predictions of IgG kinetics with varying degrees of polynomial
IgG_sensitivity = sensitivityPlot + ylab("log(IgG) titres") + ylim(-3,10) + scale_color_manual(values=c("#00B934","#e9b141")) + scale_fill_manual(values=c("#00B934","#e9b141")) 
IgG_sensitivity2 = IgG_sensitivity

for (eachPoly in c(3:5)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 2)
    IgG_sensitivity = IgG_sensitivity + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

for (eachPoly in c(3:7)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 2)
    IgG_sensitivity2 = IgG_sensitivity2 + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

# predictions of IgM kinetics with varying degrees of polynomial
IgM_sensitivity = sensitivityPlot + ylim(-3,10) +  ylab("log(IgM) titres") +scale_color_manual(values=c("#5D99FF","#F561E2")) + scale_fill_manual(values=c("#5D99FF","#F561E2")) 
IgM_sensitivity2 = IgM_sensitivity

for (eachPoly in c(3:5)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 3)
    IgM_sensitivity = IgM_sensitivity + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

for (eachPoly in c(3:7)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 3)
    IgM_sensitivity2 = IgM_sensitivity2 + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

# predictions of IgG:IgM kinetics with varying degrees of polynomial
ratio_sensitivity = sensitivityPlot + ylim(-8,10) +  ylab("log(IgG:IgM) titres") + scale_color_manual(values=c("#F8746B","#00BEC3")) + scale_fill_manual(values=c("#F8746B","#00BEC3")) 
ratio_sensitivity2 = ratio_sensitivity

for (eachPoly in c(3:5)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 4)
    ratio_sensitivity = ratio_sensitivity + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

for (eachPoly in c(3:7)) {
    eachFit = fitResponse(fit_all3_sensitivity, nPoly=eachPoly, view = 4)
    ratio_sensitivity2 = ratio_sensitivity2 + geom_ribbon(data = eachFit$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.05)+
    geom_line(data = eachFit$final, linewidth = 0.5, color = colorList[eachPoly-2])
}

# joining of all sensitivity plots
ggarrange(
    IgG_sensitivity + theme(legend.position = "none"),
    IgG_sensitivity2 + theme(legend.position = "none"),
    IgM_sensitivity + theme(legend.position = "none"), 
    IgM_sensitivity2 + theme(legend.position = "none"),
    ratio_sensitivity + theme(legend.position = "none"),
    ratio_sensitivity2 + theme(legend.position = "none")
    , ncol = 2
    , nrow = 3
)


# Models fitted with both acute and convalescent antibody sample data

# fitting the model with log(IgG), log(IgM) and log(IgG:IgM)
fitConv = 
    dat %>%
    filter(
        !(is.na(ConvalescentIgM_2ndblooddraw)),
        !(is.na(ConvalescentIgG_2ndblooddraw)),
        !(is.na(logRatioConv))
    ) %>%
    mutate(
        pPrimary_init = as.integer(Age < 8),
        logIgM = log2(AcuteIgM_1stblooddraw + 1),
        logIgG = log2(AcuteIgG_1stblooddraw + 1),
        logIgMConv = log2(ConvalescentIgM_2ndblooddraw + 1),
        logIgGConv = log2(ConvalescentIgG_2ndblooddraw + 1),
        DaysConv = 14,
    ) %>%
    select(SubjectNo, Days, logIgM, logIgG, logRatio, pPrimary_init, DaysConv, logIgMConv, logIgGConv, logRatioConv) %>% #Age, 
    fitResponseConv(nPoly=5)


# improvement of likelihood over the iterations
plot(fitConv$Q, ylab = 'LL', xlab = 'Iteration', type = 'l')


colourList = c()
colourList = c(colourList,"#00B934","#e9b141") # for IgG
colourList = c(colourList, "#5D99FF","#F561E2") # for IgM
colourList = c(colourList, "#F8746B","#00BEC3") # for IgG:IgM

# fitted average lines of the kinetics, along with initial model's prediction
ggplot(mapping = aes(x = Days, y = logIg, color = Response, fill = Response))+
    scale_color_manual(values=colourList) +
    scale_fill_manual(values=colourList) +
    geom_hline(yintercept = 0, linetype = 3)+
    # initial splines
    #geom_line(data = fitConv$init, linetype = 2, linewidth = 0.6)+
    # final fitted splines
    geom_ribbon(data = fitConv$final, aes(ymin = lowerVal, ymax = upperVal), color = NA, alpha = 0.1)+
    geom_line(data = fitConv$final, linewidth = 0.5)+
    geom_line(data = fit_all3$final, linetype = 4)+
    ylab("log(Ig) titres")+
    scale_x_continuous(breaks = seq(0, 14, by = 1))


# filtering of dataframe for use with convalescent antibody data 
datConv = dat %>%
    filter(
        !(is.na(ConvalescentIgM_2ndblooddraw)),
        !(is.na(ConvalescentIgG_2ndblooddraw)),
        !(is.na(logRatioConv))
        ) %>%
    mutate(pPrimary = fitConv$pPrimary[(seq(1, length(fitConv$pPrimary), 2))],
    pPrimaryMAD = fitConv$P1)

# sample size of patients with samples collected on days 8 and 9
print(nrow(dat[(dat$Days == 8),]))
print(nrow(dat[(dat$Days == 9),]))
print(nrow(datConv[(datConv$Days == 8),]))
print(nrow(datConv[(datConv$Days == 9),]))

# bimodal histogram of P(Primary) for model with convalescent data
datConv %>%
    ggplot(aes(x = pPrimary))+
    geom_histogram(aes(y = after_stat(density)), color = NA, fill = '#dddddd', bins = 100)+
    ylab("P(Primary)")


# plots of P(Primary) x Disease Severity and Age
ggarrange(
    datConv %>%
        ggplot(aes(x = ClinicalDiagnosis, y = pPrimary, fill = ClinicalDiagnosis)) +
        scale_fill_brewer(palette="Set2") +
        geom_violin(adjust=1.5) +
        geom_boxplot(width = 0.1, fill = NA) + #+#adjust=1.5, 
        ylab("P(Primary)")+
        xlab("Disease Severity")+
        guides(fill=guide_legend(title="Disease Severity"))
    , datConv %>%
        mutate(
            AgeGroup = cut(Age, breaks = c(0,6,12,Inf), include.lowest = T)
            ) %>%
        ggplot(aes(x = AgeGroup, y = pPrimary, fill = AgeGroup)) +
        geom_violin(adjust=1.5 ) +
        geom_boxplot(width = 0.1) + #+#adjust=1.5, 
        ylab("P(Primary)")+
        scale_x_discrete(name = "Age (years)", labels = c("1.5-6", "6-12", ">12")) +
        scale_fill_discrete(name = "Age (years)", labels = c("1.5-6", "6-12", ">12")) 
    , ncol = 1
    , nrow = 2
)


# calculation of mean absolute difference using convalescent model
datConv %>% 
    group_by(Days) %>%
    summarise(mad = sum(abs(0.5 - pPrimaryMAD)) / n()) %>%
    add_row(Days = 14, mad = sum(abs(0.5 - datConv$pPrimaryMAD)) / nrow(datConv)) %>%
    ggplot(aes(x = Days, y = mad)) +
    expand_limits(y=0) +
    geom_line(linewidth = 0.6) +
    scale_x_continuous(breaks = seq(0, 14, by = 1)) +
    ylab("MAD")
