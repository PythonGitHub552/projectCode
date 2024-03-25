avgResponseDf = function(fit, Sd, responseList, ci = c(0.025,0.975)){
    q = qnorm(ci)
    seq_along(fit) %>%
    lapply(function(iResponse){
        out = tibble(Days = 0:9)
        out %>%
            mutate(
                logIg = predict(fit[[iResponse]], newdata = out)
                , lowerVal = logIg + Sd[iResponse] * q[1]
                , upperVal = logIg + Sd[iResponse] * q[2]
                , Response = responseList[iResponse]
            )

    }) %>%
    do.call(what = rbind)    
}


avgResponseDfConv = function(fit, Sd, ci = c(0.025,0.975)){
    q = qnorm(ci)
    seq_along(fit) %>%
    lapply(function(iResponse){
        out = tibble(Days = c(0:9,14))
        out %>%
            mutate(
                logIg = predict(fit[[iResponse]], newdata = out)
                , lowerVal = logIg + Sd[iResponse] * q[1]
                , upperVal = logIg + Sd[iResponse] * q[2]
                , Response = c('IgM Primary','IgM Secondary','IgG Primary','IgG Secondary', 'IgG:IgM Primary', 'IgG:IgM Secondary')[iResponse]
            )

    }) %>%
    do.call(what = rbind)    
}


dnormCalc = function(logRatio, predicted, sd){
    dnormDf <- data.frame(logRatio = logRatio, predicted = predicted, sd = sd)
    #View(dnormDf)
    P <- apply(dnormDf, 1, function(x) {dnorm(x[1], x[2], x[3])})
    #apply(x, margin, func)
    #dTest$logRatio, predictedTest[[1]], Sd[1]
    P
}


fitResponse = function(dTrain, nPoly, inputIgM = TRUE, inputIgG = TRUE, inputRatio = TRUE, view = 1){
    # view indicates which variables should be visualised in plots
    #   1 = all variables being fitted
    #   2 = only IgG
    #   3 = only IgM
    #   4 = only IgG:IgM
    #   5 = IgG and IgM

    # use age group as initial value for response type of individuals
    pPrimary = dTrain$pPrimary_init

    # modified sum only considers finite values
    sum.finite <- function(x) {
        sum(x[is.finite(x)])
    }

    responseList = c()
    if (inputIgM & view %in% c(1,3,5)) {
        responseList = c(responseList, 'IgM Primary','IgM Secondary')
    }
    if (inputIgG & view %in% c(1,2,5)) {
        responseList = c(responseList, 'IgG Primary','IgG Secondary')
    }
    if (inputRatio & view %in% c(1,4)) {
        responseList = c(responseList, 'IgG:IgM Primary', 'IgG:IgM Secondary')
    }

    Q = 0; k = 1
    while(TRUE){

        if(k > 1){
            # ... E step
            pPrimary = L1/(L1+L2)            
        }        

            # ... M step


        # fit splines for each group   
        fit = c()  

        if (inputIgM) {
            fitIgM =  
                list(pPrimary, 1-pPrimary) %>%
                lapply(function(w){
                    lm(logIgM ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain)            
                })
            if (view %in% c(1,3,5)) {
                fit = c(fit, fitIgM) 
            }
        }

        if (inputIgG) {
            fitIgG = 
                list(pPrimary, 1-pPrimary) %>%
                lapply(function(w){
                    lm(logIgG ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain) 
                })
            if (view %in% c(1,2,5)) {
                fit = c(fit, fitIgG) 
            }
        }
        
        if (inputRatio) {
            fitRatio = 
                list(pPrimary, 1-pPrimary) %>%
                lapply(function(w){
                    lm(logRatio ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain) 
                })
            if (view %in% c(1,4)) {
                fit = c(fit, fitRatio) 
            }
        }

        # make predictions from each of the splines and calculate SD of each group
        Sd = NULL

        if (inputIgM) {
            predictedIgM = lapply(fitIgM, function(Fit){ predict(Fit, newdata = dTrain)}) 
            SdIgM = sqrt(sum.finite(pPrimary     * (dTrain$logIgM - predictedIgM[[1]])^2) / sum.finite(pPrimary))
            SdIgM[2] = sqrt(sum.finite((1-pPrimary) * (dTrain$logIgM - predictedIgM[[2]])^2) / sum.finite(1 - pPrimary))
            Sd = SdIgM
        }

        if (inputIgG) {
            predictedIgG = lapply(fitIgG, function(Fit){ predict(Fit, newdata = dTrain)}) 
            SdIgG = sqrt(sum.finite(pPrimary     * (dTrain$logIgG - predictedIgG[[1]])^2) / sum.finite(pPrimary))
            SdIgG[2] = sqrt(sum.finite((1-pPrimary) * (dTrain$logIgG - predictedIgG[[2]])^2) / sum.finite(1 - pPrimary))
            if (length(Sd) == 0) {
                Sd = SdIgG
            } else {
                Sd[3] = SdIgG[1]
                Sd[4] = SdIgG[2]
            }
        }

        if (inputRatio) {
            predictedRatio = lapply(fitRatio, function(Fit){ predict(Fit, newdata = dTrain)}) 
            SdRatio = sqrt(sum.finite(pPrimary     * (dTrain$logRatio - predictedRatio[[1]])^2) / sum.finite(pPrimary))
            SdRatio[2] = sqrt(sum.finite((1-pPrimary) * (dTrain$logRatio - predictedRatio[[2]])^2) / sum.finite(1 - pPrimary))
            if (length(Sd) == 0) {
                Sd = SdRatio
            } else {
                Sd[5] = SdRatio[1]
                Sd[6] = SdRatio[2]
            }
        }  

        # calculate proportions in each group
        prop1 = sum(pPrimary)/length(pPrimary)
        prop2 = 1 - prop1

        # keep initial states for comparison
        if(k==1){
            #out = list(init = avgResponseDf(fit[1:4], Sd[1:4]))
            out = list(init = avgResponseDf(fit, Sd, responseList))
        }

        # calculate log-likelihood of primary (L1) and secondary (L2) infections
        k = k + 1
        L1 = prop1
        L2 = prop2

        if (inputIgM) {
            L1 = L1 * dnorm(dTrain$logIgM, predictedIgM[[1]], SdIgM[1])
            L2 = L2 * dnorm(dTrain$logIgM, predictedIgM[[2]], SdIgM[2])
        }

        if (inputIgG) {
            L1 = L1 * dnorm(dTrain$logIgG, predictedIgG[[1]], SdIgG[1])
            L2 = L2 * dnorm(dTrain$logIgG, predictedIgG[[2]], SdIgG[2])
        }

        if (inputRatio) {
            L1 = L1 * dnorm(dTrain$logRatio, predictedRatio[[1]], SdRatio[1])
            L2 = L2 * dnorm(dTrain$logRatio, predictedRatio[[2]], SdRatio[2])
        }

        Q[k] <- sum(log(L1+L2))


        # stop if LL difference is lower than threshold
        if(abs(Q[k]-Q[k-1]) < 1e-6){ break }
    }
    
    out$final = avgResponseDf(fit, Sd, responseList)

    # calculations for MAD plot
    P1 = 1

    if (inputIgM) {
        P1 = P1 * dnorm(dTrain$logIgM, predictedIgM[[1]], SdIgM[1])
    }

    if (inputIgG) {
        P1 = P1 * dnorm(dTrain$logIgG, predictedIgG[[1]], SdIgG[1])
    }

    if (inputRatio) {
        P1 = P1 * dnorm(dTrain$logRatio, predictedRatio[[1]], SdRatio[1])
    }

    c(out, list(
        pPrimary = pPrimary, Q = Q[-1], fit = fit, nPoly = nPoly, P1 = P1
    ))
}



fitResponseConv = function(dTrain, nPoly){

    # use age group as initial value for response type of individuals
    # doubling the length of pPrimary to account for convalescent data
    pPrimary = rep(dTrain$pPrimary_init, each = 2)

    # modified sum only considers finite values
    sum.finite <- function(x) {
        sum(x[is.finite(x)])
    }

    # transforming dTrain log IgM, IgG and IgG:IgM into long format
    dTrain_logIgM_long <- dTrain %>%
        select(-c(Days, DaysConv, logIgG, logIgGConv, logRatio, logRatioConv)) %>%
        pivot_longer(cols = c(logIgM, logIgMConv), names_to = "Timing", values_to = "logIgM") %>%
        mutate(Timing = recode(Timing, 'logIgM' = 'Acute', 'logIgMConv' = 'Convalescent'))

    dTrain_logIgG_long <- dTrain %>%
        select(-c(Days, DaysConv, logIgM, logIgMConv, logRatio, logRatioConv, pPrimary_init)) %>%
        pivot_longer(cols = c(logIgG, logIgGConv), names_to = "Timing", values_to = "logIgG") %>%
        mutate(Timing = recode(Timing, 'logIgG' = 'Acute', 'logIgGConv' = 'Convalescent'))

    dTrain_logRatio_long <- dTrain %>%
        select(-c(Days, DaysConv, logIgM, logIgMConv, logIgG, logIgGConv, pPrimary_init)) %>%
        pivot_longer(cols = c(logRatio, logRatioConv), names_to = "Timing", values_to = "logRatio") %>%
        mutate(Timing = recode(Timing, 'logRatio' = 'Acute', 'logRatioConv' = 'Convalescent'))

    # transforming dTrain acute and convalescent day data into long format
    dTrain_days_long <- dTrain %>%
        select(SubjectNo, Days, DaysConv) %>%
        pivot_longer(cols = c(Days, DaysConv), names_to = "Timing", values_to = "Days") %>%
        mutate(Timing = recode(Timing, 'Days' = 'Acute', 'DaysConv' = 'Convalescent'))

    # merging of both long format dataframes in a way that ensure the correct log Ig data and Days values share the same row based on Timing
    dTrain_long = merge(dTrain_logIgM_long, dTrain_logIgG_long, by = c("SubjectNo","Timing"))
    dTrain_long = merge(dTrain_long, dTrain_logRatio_long, by = c("SubjectNo","Timing"))
    dTrain_long = merge(dTrain_long, dTrain_days_long, by = c("SubjectNo","Timing"))
 
    
    Q = 0; k = 1; count = 0
    while(TRUE){

        if(k > 1){
            # ... E step
            pPrimary = rep(L1/(L1+L2), each = 2)           
        }        

            # ... M step 
        
        # fit splines for each group 
        fitIgM = 
            list(pPrimary, 1-pPrimary) %>%
            lapply(function(w){
                lm(logIgM ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain_long) 
            })

        fitIgG = 
            list(pPrimary, 1-pPrimary) %>%
            lapply(function(w){
                lm(logIgG ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain_long) 
            })

        fitlogRatio =  
            list(pPrimary, 1-pPrimary) %>%
            lapply(function(w){
                lm(logRatio ~ poly(Days, nPoly, raw = TRUE), weights = w, data = dTrain_long)   
            })
        
        
        fit = c(fitIgM, fitIgG, fitlogRatio)
       
    
        # make predictions from each of the splines
        predictedIgM = lapply(fit[1:2], function(Fit){ predict(Fit, newdata = dTrain_long)}) 
        predictedIgG = lapply(fit[3:4], function(Fit){ predict(Fit, newdata = dTrain_long)}) 
        predictedlogRatio = lapply(fit[5:6], function(Fit){ predict(Fit, newdata = dTrain_long)}) 

        # calculate SD of each group
        Sd = sqrt(sum.finite(pPrimary     * (dTrain_long$logIgM - predictedIgM[[1]])^2) / sum.finite(pPrimary))
        Sd[2] = sqrt(sum.finite((1-pPrimary) * (dTrain_long$logIgM - predictedIgM[[2]])^2) / sum.finite(1 - pPrimary))
        Sd[3] = sqrt(sum.finite(pPrimary     * (dTrain_long$logIgG - predictedIgG[[1]])^2) / sum.finite(pPrimary))
        Sd[4] = sqrt(sum.finite((1-pPrimary) * (dTrain_long$logIgG - predictedIgG[[2]])^2) / sum.finite(1 - pPrimary))
        Sd[5] = sqrt(sum.finite(pPrimary     * (dTrain_long$logRatio - predictedlogRatio[[1]])^2) / sum.finite(pPrimary))
        Sd[6] = sqrt(sum.finite((1-pPrimary) * (dTrain_long$logRatio - predictedlogRatio[[2]])^2) / sum.finite(1 - pPrimary))



        # calculate proportions in each group
        prop1 = sum(pPrimary)/length(pPrimary)
        prop2 = 1 - prop1

        # keep initial states for comparison
        if(k==1){
            out = list(init = avgResponseDfConv(fit, Sd))
        }


        # separating data into acute and convalescent to calculate log-likelihood
        dTrainAcute = dTrain_long[dTrain_long$Timing == "Acute",]
        dTrainConv = dTrain_long[dTrain_long$Timing == "Convalescent",]

        # separates predictedIgM, predictedIgG and predictedlogRatio into predicted primary acute, primary convalescent, secondary acute and secondary convalescent
        predictedIgM_longer = c(list(predictedIgM[[1]][(seq(1, length(predictedIgM[[1]]), 2))]), 
                                list(predictedIgM[[1]][(seq(2, length(predictedIgM[[1]]), 2))]),
                                list(predictedIgM[[2]][(seq(1, length(predictedIgM[[2]]), 2))]),
                                list(predictedIgM[[2]][(seq(2, length(predictedIgM[[2]]), 2))]))

        predictedIgG_longer = c(list(predictedIgG[[1]][(seq(1, length(predictedIgG[[1]]), 2))]), 
                                list(predictedIgG[[1]][(seq(2, length(predictedIgG[[1]]), 2))]),
                                list(predictedIgG[[2]][(seq(1, length(predictedIgG[[2]]), 2))]),
                                list(predictedIgG[[2]][(seq(2, length(predictedIgG[[2]]), 2))]))

        
        predictedlogRatio_longer = c(list(predictedlogRatio[[1]][(seq(1, length(predictedlogRatio[[1]]), 2))]), 
                                list(predictedlogRatio[[1]][(seq(2, length(predictedlogRatio[[1]]), 2))]),
                                list(predictedlogRatio[[2]][(seq(1, length(predictedlogRatio[[2]]), 2))]),
                                list(predictedlogRatio[[2]][(seq(2, length(predictedlogRatio[[2]]), 2))]))

        # calculate log-likelihood of primary (L1) and secondary (L2) infections
        k = k + 1
        L1 = prop1 * dnorm(dTrainAcute$logIgM, predictedIgM_longer[[1]], Sd[1]) * dnorm(dTrainConv$logIgM, predictedIgM_longer[[2]], Sd[1]) * 
            dnorm(dTrainAcute$logIgG, predictedIgG_longer[[1]], Sd[3]) * dnorm(dTrainConv$logIgG, predictedIgG_longer[[2]], Sd[3]) *
            dnorm(dTrainAcute$logRatio, predictedlogRatio_longer[[1]], Sd[5]) * dnorm(dTrainConv$logRatio, predictedlogRatio_longer[[2]], Sd[5])
        L2 = prop2 * dnorm(dTrainAcute$logIgM, predictedIgM_longer[[3]], Sd[2]) * dnorm(dTrainConv$logIgM, predictedIgM_longer[[4]], Sd[2]) * 
            dnorm(dTrainAcute$logIgG, predictedIgG_longer[[3]], Sd[4]) * dnorm(dTrainConv$logIgG, predictedIgG_longer[[4]], Sd[4]) *
            dnorm(dTrainAcute$logRatio, predictedlogRatio_longer[[3]], Sd[6]) * dnorm(dTrainConv$logRatio, predictedlogRatio_longer[[4]], Sd[6])

        Q[k] <- sum(log(L1+L2))

        # stop if LL difference is lower than threshold
        if(abs(Q[k]-Q[k-1]) < 1e-6){ break }
    }
    
    out$final = avgResponseDfConv(fit,Sd)

    # calculations for MAD plot
    P1 = dnorm(dTrainAcute$logIgM, predictedIgM_longer[[1]], Sd[1]) * dnorm(dTrainConv$logIgM, predictedIgM_longer[[2]], Sd[1]) * 
            dnorm(dTrainAcute$logIgG, predictedIgG_longer[[1]], Sd[3]) * dnorm(dTrainConv$logIgG, predictedIgG_longer[[2]], Sd[3]) *
            dnorm(dTrainAcute$logRatio, predictedlogRatio_longer[[1]], Sd[5]) * dnorm(dTrainConv$logRatio, predictedlogRatio_longer[[2]], Sd[5])

    c(out, list(
        pPrimary = pPrimary, Q = Q[-1], fit = fit, nPoly = nPoly, P1 = P1
    ))
}