## This R-file contains functions required for running the interference modeling script ("IM.Rmd")




### Write function that extracts required spectral information
extract_spectral_features_in_parallel <- function(rawfilefolder_filepath = rawfilefolder_filepath,
                                                  msms,
                                                  scan_number = "Scan.number",
                                                  rawfile = "Raw.file",
                                                  charge = "Charge",
                                                  reporter_ions,
                                                  mass_error_tolerance = 0.005){
  ## This function requires: 
  ## 1) The file path to where Thermo .raw-files are stored + respective raw file tables provided by rawStallion; specified by parameter "rawfilefolder_filepath". 
  ## 2) A PSM-table data frame (for example  MQ search result "msms.txt"); specified by the parameter "msms". Required column names are to be specified via parameters "scan_number", "rawfile" and "charge".
  ## 3) A reporter ion object as implemented in the R Bioconductor package "MSnBase"; specified by parameter "reporter_ions".
  ## note: This function takes some time to finish (~30 min if there exists a separate core for each raw file to be processed. Else even longer)
  ## note: Choose mass error tolerance according to the distribution of mass errors in your data. Per default, using 0.005 Th.
  ## This function outputs: A data frame containing PSM-specific metrics like the noise value MS1 isolation windows, the precursor purity fraction PPF (similar to MaxQuant's PIF), reporter ion intensities, and more.
  
  ## Loop in parallel over the different raw files in the experiment
  df_res <- foreach (r = unique(msms[,rawfile]), .combine = rbind, .packages = c("tidyverse","foreach","MSnbase")) %dopar% {
    
    ## Create unique identifier for each PSM
    key_msms <- paste0(msms[,rawfile], "_" , msms[,scan_number])
  
    ## Extract list of all rawfiles + noise tables found in rawfilefolder_filepath
    rawfiles_filepath <- list.files(rawfilefolder_filepath, full.names = TRUE)
    rawfiles_filepath
  
    ## Define helper function that extracts mass, intensity and noise for a given scan number. 
    extract_single_spectrum <- function(df, scannr, window=NULL){
      df_spec <- df[df$ScanNumber == scannr,]
      res <- data.frame(mz = df_spec$Mass, i = df_spec$Intensity , n = df_spec$Noise)
      if(!is.null(window)){
        lwr <- window[1]
        upr <- window[2]
        res <- res[res$mz > lwr & res$mz < upr,]
      } 
      return(res)
    }
    
    ## Find corresponding filepath of rawfile_r and noise table of rawfile_r
    filepath_index_r <- grep(x= rawfiles_filepath, pattern = paste0(r, "[-]index[.]tsv"), value = TRUE)
    filepath_noisevalues_r <- grep(x= rawfiles_filepath, pattern = paste0(r, "[-]noise[.]tsv"), value = TRUE)
    
    ## Read in index file for rawfile r (which includes metadata like scan type etc). Also save some variables
    df_index_r <- read.delim(file=filepath_index_r, sep="\t", header=TRUE)
    df_index_r$CV <- df_index_r$CV %>% as.character(.) %>% replace_na(., replace = "")
    scannr_r <- df_index_r$ScanNumber
    scannr_ms1_r <- scannr_r[df_index_r$MsOrder == 1]
    cv_ms1_r <- df_index_r$CV[df_index_r$MsOrder == 1]
    
    ## Read in noise table for rawfile r (provided by rawStallion)
    df_noise_r <- read.delim(file=filepath_noisevalues_r, sep="\t", header=TRUE)
    
    ## Find all PSM IDs (unique PSM identifier) belonging to rawfile r
    key_msms_r <- key_msms[msms[,rawfile]==r]
    
    
    ## Within rawfile r, loop over each unique PSM i to ...
    ## A) extract reporter ion intensities.
    ## B) extract noise values etc. and all peaks in the isolation window range of adjacent MS1-scans. 
    ## C) calculate PPF (precursor purity fraction) and TIW (total intensity in the isolation window).
    df_res_r <- foreach (i = key_msms_r, .combine = rbind, .packages = c("MSnbase")) %do% {
      
      # Initiate output dataframe w/o reporter ion intensities for PSM i in rawfile r 
      df_res_r_i <- data.frame(key_msms=i,
                               rawfileCharge = numeric(1),
                               noiseValue = numeric(1),
                               PPF = numeric(1),
                               TIW = numeric(1),
                               PIC = numeric(1),
                               precursorMz = numeric(1),
                               minIntensity_MS2 = numeric(1),
                               parent_MS1 = "preceding",
                               compensationVoltage=character(1))
      
      
      # Initiate reporter intensity vector for PSM i in rawfile r
      reporter_intensities <- matrix(numeric(0), nrow=1, ncol=length(reporter_ions))
      rownames(reporter_intensities) <- i
      colnames(reporter_intensities) <- paste0("reporters_", reporter_ions@reporterNames)
      
      
      # Save data of PSM i from msms.txt and extract scan number and precursor charge 
      df_msms_i <- msms[key_msms==i,,drop=FALSE]
      scannr_i <- df_msms_i[,scan_number]
      charge_i <- df_msms_i[,charge]
      
      # Save MS2 spectrum i from index file
      df_index_r_i <- df_index_r[df_index_r$ScanNumber == scannr_i, ]
      
      # Save rawfileCharge and precursorMZ and compensation voltage
      df_res_r_i[,"rawfileCharge"] <- df_index_r_i$PrecursorCharge
      df_res_r_i[,"precursorMz"] <- df_index_r_i$PrecursorMz
      df_res_r_i[,"compensationVoltage"] <- df_index_r_i$CV
      
      
      ## A) Extract reporter intensities of PSM i, identical in code to MSnBase quantify(method="max"). Also extract PIC (Total Peptide Ion Current in MS2 scan), and min Intensity of MS2 scans
      
      # Extract whole ms2 spectrum
      df_ms2spec_i <- extract_single_spectrum(df = df_noise_r, scannr = scannr_i)
      
      # Calulate selection width for quantification
      m <- reporter_ions@mz
      lwr <- m - reporter_ions@width
      upr <- m + reporter_ions@width
      
      # Initiate reporter ion intensity vector for PSM i
      reporterInt_i <- numeric(length(m))
      
      # Extract reporter intensities and save them
      for (l in 1:length(m)){
        region <-  df_ms2spec_i [df_ms2spec_i$mz > lwr[l] & df_ms2spec_i$m < upr[l], ]
        if(nrow(region) == 0){
          reporterInt_i[l] <- NA
        } else reporterInt_i[l] <- max(region$i)
      }
      reporter_intensities[i,] <- reporterInt_i
      
      # Extract PIC (Total Peptide Ion Current in the MS2-scan)
      try(df_res_r_i[,"PIC"] <- sum(df_ms2spec_i$i) - sum(reporter_intensities[i,], na.rm=TRUE))
      
      # Extract minIntensity_MS2 (minimum intensity found in the MS2-scan)
      df_res_r_i[,"minIntensity_MS2"] <- min(df_ms2spec_i$i)
      
      
      ## B) Check for parent MS1 scan; extract noise value in isolation window; extract total isolation window make-up of previous and following ms1 scans
      
      # Find scan index of previous and following MS1 spectra within raw file with identical CV
      scannr_ms1_i <- setNames(numeric(2), c("MS1_previous","MS1_following"))
      scannr_ms1_i["MS1_previous"] <- max(scannr_ms1_r[scannr_ms1_r < scannr_i & cv_ms1_r == df_index_r_i$CV])
      scannr_ms1_i["MS1_following"] <- min(scannr_ms1_r[scannr_ms1_r > scannr_i & cv_ms1_r == df_index_r_i$CV])
      scannr_ms1_i <- scannr_ms1_i[!is.infinite(scannr_ms1_i)]
      
      # Calculate distances in retention time between ms2 and adjacent ms1 spectra
      if (length(scannr_ms1_i) == 2){
        rt_ms1_i <- setNames(numeric(2), c("MS1_previous","MS1_following")) 
        rt_ms1_i["MS1_previous"] <- df_index_r$RetentionTime[df_index_r$ScanNumber == scannr_ms1_i["MS1_previous"]]
        rt_ms1_i["MS1_following"] <- df_index_r$RetentionTime[df_index_r$ScanNumber == scannr_ms1_i["MS1_following"]]
        rt_ms2_i <- df_index_r$RetentionTime[df_index_r$ScanNumber==scannr_i]
        dist_ms1_ms2_i <- abs(rt_ms1_i - rt_ms2_i)
        relative_weights_i <- 1 - dist_ms1_ms2_i/sum(dist_ms1_ms2_i)
      } else {relative_weights_i <- 1}
      
      # Calculate isolation window range 
      target_mz <- df_index_r_i$PrecursorMz
      window_lwr <- target_mz - as.numeric(df_index_r_i$IsolationWidth)/2
      window_upr <- target_mz + as.numeric(df_index_r_i$IsolationWidth)/2
      
      # Go over both preceding and following MS1 scan 
      for (p in 1:length(scannr_ms1_i)){
        df_window_p <- extract_single_spectrum(df = df_noise_r, scannr = scannr_ms1_i[p], window=c(window_lwr, window_upr))
        
        # Handle the preceding MS1 scan
        if (p == 1){
          
          # Catch the exception where the precursor peak is not found in the last MS1 spectrum with identical CV
          if (!any(abs(df_window_p$mz - target_mz) < mass_error_tolerance)){
            warning(paste0("Scan number",scannr_i, ": precursor peak was not found in the isolation window of the last MS1 scan."))
            scannr_ms1_secondlast <- max(scannr_ms1_r[scannr_ms1_r < scannr_ms1_i["MS1_previous"] & cv_ms1_r == df_index_r_i$CV])
            df_window_secondlast <- extract_single_spectrum(df = df_noise_r, scannr = scannr_ms1_secondlast, window=c(window_lwr, window_upr))
            
            # Check if the precursor ion peak can be found in the second last MS1 spectrum
            if (any(abs(df_window_secondlast$mz - target_mz) < mass_error_tolerance)){
              warning("However, it was found in the second last MS1 scan. \n ")
              
              # If TRUE, instead use information of second last spectrum instead. 
              df_window_p <- df_window_secondlast
              df_res_r_i[,"parent_MS1"] <- "second_last"
              
            } else {
              warning("Nor was it found in the second last MS1 scan.\n")
              
              # If FALSE, impute precursor intensity  
              df_ms1spec <- extract_single_spectrum(df = df_noise_r, scannr = scannr_ms1_i[p], window=c(window_lwr - 1000, window_upr + 1000))
              noise_i <- try(mean(df_ms1spec$n))
              min_i <- try(min(df_ms1spec$i))
              df_window_p <- rbind(df_window_p, data.frame(mz=target_mz, i = min_i, n = noise_i))
              df_res_r_i[,"parent_MS1"] <- "not found at all"
            }
          }
          
          ## Weight intensities by relative distance. Also save noise values
          df_window_p$i <- df_window_p$i*relative_weights_i[p]
          df_res_r_i[,"noiseValue"] <- mean(df_window_p$n)
          df_window <- df_window_p
        }
        
        # Handle the following MS1 scan
        if (p == 2){
          
          # Check if isolation window is totally empty (i.e. does not show any ion peaks)
          if (nrow(df_window_p) == 0){
            
            # If TRUE: skip it. This improved results in my experience.
            df_window_p <- NULL
            
          } else {
            # If FALSE, check if precursor can be found within a certain mass tolerance
            proposed_precursor_mass <- df_window_p$mz[which.min(abs(df_window_p$mz - target_mz))]
            proposed_correction <-  target_mz - proposed_precursor_mass
            if (abs(proposed_correction) < mass_error_tolerance){
              # If TRUE: slightly adjust mz
              df_window_p$mz <- df_window_p$mz + proposed_correction
            } else {
              # If FALSE, further decision depends on whether precursor ion was already found in the last MS1 scan
              if (df_res_r_i[, "parent_MS1"] == "preceding"){
                # If TRUE: impute with precursor at minimum intensity in vicinity. Precursor intensity probably too low to be detected
                df_ms1spec <- try(extract_single_spectrum(df = df_noise_r, scannr = scannr_ms1_i[p], window=c(window_lwr - 5, window_upr + 5)))
                min_i <- try(min(df_ms1spec$i))
                df_imputed_prec <- data.frame(mz=target_mz, i = min_i*relative_weights_i[p], n = NA)
                df_window_p <- rbind(df_window_p, df_imputed_prec)
              } else {
                # If FALSE: better skip. Unreliable data
                df_window_p <- NULL
              }
            }
          }
          
          ## Merge info of both isolation windows (previous and following MS1 scan)
          df_window <- rbind(df_window, df_window_p)
        }
      }
      
      
      ## C) Calculate the total intensity in the isolation window (TIW), infer +1/-1 isotopes and finally calculate the precursor purity fraction (PPF)
      
      # Calculate TIW (total intensity in the isolation window)
      df_res_r_i[, "TIW"] <- sum(df_window$i)
      
      # Calculate which of the peaks corresponds to the precursor peptide ion (i.e. the mass that triggered the MS/MS)
      ind_prec <- which(abs(df_window$mz - target_mz) < mass_error_tolerance)
      
      # Find out which of the peaks correspond to +1 & -1 isotopes of the precursor peptide ion, then combine this index with the precursor index
      heavy_mz <- target_mz + 1/charge_i
      light_mz <- target_mz - 1/charge_i
      ind_isotope <- which( (abs(df_window$mz - heavy_mz ) <  mass_error_tolerance/2) | (abs(df_window$mz - light_mz) <  mass_error_tolerance/2) )
      if (length(ind_isotope) > 0) {ind_prec <- c(ind_prec, ind_isotope)}
      
      # Calculate PFF
      PPF_i <- sum(df_window$i[ind_prec]) / sum(df_window$i, na.rm=TRUE)
      df_res_r_i[, "PPF"] <- PPF_i 
      
      ## Return an intermediate result for each PSM i in rawfile r 
      df_res_r_i_return <- cbind(df_res_r_i, reporter_intensities)
    }
    
    ## Delete memory-intensive objects
    rm(list=c("df_noise_r", "df_index_r", "scannr_r", "scannr_ms1_r", "cv_ms1_r", "key_msms_r"))
    gc()
    
    ## Return an intermediate result for each rawfile r
    df_res_r_return <- df_res_r
  }
  return(df_res)
}  






### Write function that calculates raw-file specifc peptide density (x = mz, y = retention time) for each PSM in the PSM-table.
calculate_peptide_density <- function(msms,
                                      rawfile = "Raw.file__CV",
                                      retentionTime = "Retention.time",
                                      precursorMZ = "precursorMz", 
                                      modifiedSequence = "Modified.sequence",
                                      charge = "Charge",
                                      precursorIntensity = "precursorIntensity"){
  
  ## This function requires:
  ## 1) a PSM-table (for example MQ search result "msms.txt"); specified by the parameter msms as data frame. 
  ## note: This function calculates raw-file specific 2D kernel density estimates using unique peptide IDs in the PSM-table (unique in combination of modified sequence and charge).
  ##       This density is first evaluated on a 300 x 300 grid, with dimensions x = precursorMz, y = retention time, and z=density; 
  ##       and subsequently interpolated at distinct x and y coordinates for each observation (PSM) in the PSM-table (msms) to generate a PSM-wise density estimate.
  ## note: Per default, the "proximity" is defined as the top 20 closest unique peptides in the 2D euclidian plane of x = precursorMz, y = retention time.
  
  ## Ensure that NAs in the column precursorIntensity are set to 0
  msms[is.na(msms[,precursorIntensity]), precursorIntensity] <- 0
  
  ## Create unique rownames for msms
  rownames(msms) <- as.character(1:nrow(msms))
  
  ## Create unique peptide characteristics column by concatenating modified sequence and charge for each PSM
  msms$unique_peptide <- paste0(msms[,modifiedSequence],"_",msms[,charge])
  
  ## Initiate density estimates as an empty vectors 
  peptideDensity  <- numeric(nrow(msms))
  
  ## Calculate min_rt/max_rt and min_mz/max_mz of the experiment. This will determine the span of the x and y-dimension of the density grid. 
  min_mz <- min(msms[,precursorMZ])
  max_mz <- max(msms[,precursorMZ])
  min_rt <- min(msms[,retentionTime])
  max_rt <- max(msms[,retentionTime])
  
  
  ## Loop through all raw files to calculate the unique peptide density, and interpolate for measured PSMs.
  for (r in unique(msms[,rawfile])){
    print(r)
    
    ## Extract subdataframe for each rafile.
    msms_rawfile_r <- msms[msms[,rawfile] == r , , drop=FALSE]
    rownames(msms_rawfile_r) <- 1:(nrow(msms_rawfile_r))
    
    ## Reduce subdataframe msms_rawfile_r to only unique sequence + charge PSMs (when more than one entry, take PSM with highest precursor intensity). This non-redundant subdataframe will be used to calculate densities. 
    ind_keep_msms_rawfile_r <- numeric(0)
    for (j in unique(msms_rawfile_r$unique_peptide)){
      msms_rawfile_r_j <- msms_rawfile_r[msms_rawfile_r$unique_peptide == j,,drop=FALSE]
      if (nrow(msms_rawfile_r_j) > 1 ){
        ind_keep_msms_rawfile_r <- append(ind_keep_msms_rawfile_r , values = rownames(msms_rawfile_r_j)[which.max(msms_rawfile_r_j[,precursorIntensity])]  )
      } else {
        ind_keep_msms_rawfile_r <- append(ind_keep_msms_rawfile_r , values = rownames(msms_rawfile_r_j))
      }
    }
    msms_rawfile_r <- msms_rawfile_r[ind_keep_msms_rawfile_r,]
    
    ## Estimate 2d kernel joint density (using gaussian kernel function with rule-of-thumb bandwidth estimation), 300*300 grid points, for rawfile r.
    kernel_grid_r <- kde2d(x = msms_rawfile_r[,precursorMZ],  y = msms_rawfile_r[,retentionTime], n = 200 , lims = c(x1 = min_mz -1 , xu = max_mz +1, y1 = min_rt -1, yu = max_rt +1 ))
    kernel_grid_r$z <- kernel_grid_r$z/(sum(kernel_grid_r$z))*1   # rescale density to sum up to 1
    
    ## Extract precursorMZ and rententionTime of all PSMs from sample/rawfile i in msms. For these x and y coordinates, we want peptide density estimates.
    bool_rawfile_r_in_msms <- msms[,rawfile] == r
    precursorMZ_r <- msms[bool_rawfile_r_in_msms, precursorMZ]
    retentionTime_r <- msms[bool_rawfile_r_in_msms, retentionTime]
    
    ## Interpolate joint density grid at the points observed for each PSM in rawfile i using the "fields" package.
    densities_r <- fields::interp.surface(obj=kernel_grid_r, loc=cbind(precursorMZ_r,retentionTime_r))
    
    ## Take square root of density to make differences in z-axis values less extreme across the density. 
    ## (Why? It worked better in my experience, i.e. results in better R^2 in the model. After all, it is just an approximation of the true hidden complexity density along the LCMS run).
    peptideDensity[bool_rawfile_r_in_msms] <- sqrt(densities_r)
    
    ## Plot 2D-densities of PSMs in rawfile i.
    par(mfrow=c(1,2))
    scatter3D(x=precursorMZ_r, y=retentionTime_r, z=peptideDensity[bool_rawfile_r_in_msms], xlab="precursorMz", ylab="retention time", main=paste0("density ", r))
    scatter3D(x=precursorMZ_r, y=retentionTime_r, z=peptideDensity[bool_rawfile_r_in_msms], xlab="precursorMz", ylab="retention time", main=paste0("density ", r), phi=90)
  }  
  return(list(peptideDensity = peptideDensity))
}  






### Write function that builds a decision tree on peptide characteristics in order to find empirical peptide classes that explain differences in Y best. Splits resulting in n < 50 will not be conducted. 
determine_pepChar_classes <- function(msms = df_msms, 
                                      min_number = 100,             ## this determines the minimum acceptable number of observations to generate a new split in the tree 
                                      Y_var = "Y",
                                      X_var = c("factorCharge_2",
                                                "factorCharge_3",
                                                "factorCharge_4",
                                                "factorLabels_1",
                                                "factorLabels_2",
                                                "factorLabels_3",
                                                "seqChar_R",
                                                "seqChar_K",
                                                "seqChar_H",
                                                "factorExtra")){
  ## This function requires:
  ## 1) A PSM-table (for example  MQ search result "msms.txt") saved as dataframe; specified by the parameter msms. 
  ##    This table should have:
  ##    a variable denoting Y (the total reporter ion intensity of the PSM); its column name specified by the parameter Y_var.
  ##    several binary variables reflecting characteristics of peptides that influence the fragmentation efficiency of the PSM; their column names specified by a string X_var. 
  
  ## Rename Y_var to Y
  names(msms)[names(msms) == Y_var] <- "Y"
  
  ## Create a unique PSM index in the PSM-table (important for reassembling results later)
  msms$PSM_id <- 1:nrow(msms)
  
  ## Create a column "pepChar" in the PSM-table
  msms$pepChar <- character(length(nrow(msms)))
  
  
  ## Write a helper function that tests an individual split and returns sum of squared residuals, res = (Y - fitted)^2
  calculate_model_res <- function(df, splitVariable){
    
    # Rename corresponding variable in function
    names(df)[names(df) == splitVariable] <- "splitVariable"
    
    # Check if splitVariable has 2 or more levels
    if (length(unique(df$splitVariable))  <= 1 ){
      return(NA)
    }
    
    # Check if splitVariable would create a split larger than the specified minimum number
    if (min(table(df$splitVariable)) < min_number){
      return(NA)
    }
    
    # Estimate model after the split
    X <- model.matrix(data = df,
                      object = ~ 0 + precursor:splitVariable + nonprecursor:rawfileCharge + noiseEstimate)
    
    # Check for all constant columns and remove them from the model matrix to prevent singular X  matrix
    bool_keep <- apply(X, MARGIN = 2, FUN=function(x){length(unique(x)) > 1})
    X <- X[,bool_keep]
    
    # Estimate the model
    model <- rlm(y=df$Y,
                 x=X,
                 psi="psi.bisquare",
                 maxit=10000)
    
    # Calculate sum of squared residuals in log2 space 
    df$log2Y <- log2(df$Y)
    df$log2fitted <- log2(model$fitted.values)
    
    # Return sum of squared residuals of log variables
    return(sum((df$log2Y - df$log2fitted)^2)  )
  }   
  
  
  ## Create classification tree as a list (although itself not list-structed, but linear. Easier to program this way)
  tree_data <- vector(mode="list", length=length(X_var)+1)
  tree_data[[1]] <- list(df=msms)
  
  ## Create the result object. All terminal nodes (leafes) will be saved here.
  res_list <- data.frame(PSM_id=NULL, pepChar=NULL)
  
  ## Run the algorithm
  iter <- length(X_var) + 1
  writeLines("Stepwise calculation of splits that best explain the variance in Y: \n")
  for (i in 1:iter){   
    list_of_splits <- tree_data[[i]]
    
    for (j in 1:length(list_of_splits)){
      split_j <- list_of_splits[[j]]
      
      ## For split_j, test all possible splitting options in variable that have not been done so far
      split_j_variables <- X_var[X_var %in% names(split_j)]
      
      ## If all variables have already been used, skip the loop and save results 
      if (length(split_j_variables)==0){
        writeLines(paste0("ran out of variables for new split", " (", unique(split_j$pepChar), ")"))
        res_list <- list.append(res_list, data.frame(PSM_id = split_j$PSM_id, pepChar = split_j$pepChar))
        next()
      }
      
      ## Calculate which potential split variable reduces the remaining model variance the most
      ssres <- setNames(numeric(length(split_j_variables)), nm = split_j_variables)
      for (s in split_j_variables){
        ssres[s] <- calculate_model_res(df=split_j, splitVariable = s)
      }
      var <- split_j_variables[which.min(ssres)]
      
      ## If there is no further split possible because of too little dichotomy, skip the loop and save results
      if(length(var) == 0){
        writeLines(paste0("not enough observations to split further", " (", unique(split_j$pepChar), ")"))
        res_list <- list.append(res_list, data.frame(PSM_id = split_j$PSM_id, pepChar = split_j$pepChar))
        next()
      }
      
      ## Create the two new splits, and update split history in the variable "pepChar"
      split_j_1 <- split_j[ split_j[,var] == unique(split_j[,var])[1], ]
      split_j_2 <- split_j[ split_j[,var] == unique(split_j[,var])[2], ]
      split_j_1$pepChar <- paste0(split_j_1$pepChar, "_", split_j_1[,var])
      split_j_2$pepChar <- paste0(split_j_2$pepChar, "_", split_j_2[,var])
      split_j_1[[var]] <- NULL
      split_j_2[[var]] <- NULL
        
      writeLines(unique(split_j_1$pepChar))
      writeLines(unique(split_j_2$pepChar))
        
      ## Save the new splits as new nodes in the list tree_data
      if (is.null(tree_data[[i + 1]])){
        tree_data[[i + 1]]  <- list(split_j_1)
      } else {  tree_data[[i + 1]] <- list.append( tree_data[[i + 1]], split_j_1)  }
      tree_data[[i + 1]] <- list.append( tree_data[[i + 1]], split_j_2)

      writeLines("\n")  
    }
    
    ## Stop condition
    if (i < (iter ) && is.null(tree_data[[i+1]])  ){ break("no more new splits possible")}
    writeLines("----------------------------------")
    cat("\n")
  }
  
  cat("\n")
  cat("\n")
  cat("\n")
  df_res <- do.call(what=rbind,res_list)
  pepChar <- df_res[order(df_res$PSM_id),]$pepChar
  return(pepChar)
  
}




### Write interference correction function
interference_correction <- function(reporterI_matrix = reporterI_matrix, 
                                    EIL = EIL, 
                                    max_interference = 0.8,
                                    min_intensity = min_intensity_MS2){
  
  ## Ensure that reporterI_matrix is in fact a matrix
  reporterI_matrix <- as.matrix(reporterI_matrix)
  
  ## Plot colSums before correction
  par(mfrow=c(1,2))
  ymax <- max(colSums(reporterI_matrix, na.rm=TRUE))
  barplot(colSums(reporterI_matrix, na.rm=TRUE), las=2, main="total intensities before correction", border="grey", ylim=c(0, ymax), cex.names = 0.7)
  
  ## Calculate mean reporter intensity
  mean_reporterIntensity <- apply(reporterI_matrix, FUN=mean, MARGIN=1, na.rm=TRUE)
  
  ## Calculate interference signal that is to be subtracted per channel
  subtraction_intensity <- mean_reporterIntensity*pmin(EIL, max_interference)
  
  ## Do correction by subtracting estimated interference intensity from normalized data 
  reporterI_matrix[reporterI_matrix==0] <- NA
  reporterI_matrix_corrected <- sweep(reporterI_matrix, MARGIN = 1, FUN="-", STATS=subtraction_intensity)
  
  ## Replace corrected intensities smaller than the minimum observed intensity (of the MS2-spectrum) with the minimum observed intensity. Also replace NAs with minimum intensity
  min_intensity_matrix <- as.matrix(min_intensity) %*% matrix(data=1, nrow=1, ncol=ncol(reporterI_matrix))
  bool_replace_matrix <- reporterI_matrix_corrected < min_intensity_matrix | is.na(reporterI_matrix_corrected)
  reporterI_matrix_corrected[bool_replace_matrix] <- min_intensity_matrix[bool_replace_matrix]
  
  ## Plot after correction
  barplot(colSums(reporterI_matrix_corrected, na.rm=TRUE), main="total intensities after correction", col="#20217E", border="#20217E", ylim=c(0, ymax), las=2, cex.names = 0.7,)
  
  ## Add suffix "_norm_corrected" reporter ion to column names
  colnames(reporterI_matrix_corrected) <- paste0(colnames(reporterI_matrix_corrected), "__interference_corrected")
  
  ## Return corrected reporterIntensity matrix
  return(reporterI_matrix_corrected)
}





### Write LOESS normalization function
loess_norm <- function(m){
  m[m==0] <- NA
  m_log <- log2(m)  
  m_norm <- 2^normalizeBetweenArrays(m_log, method="cyclicloess", cyclic.method = "fast")
  return(m_norm)
}






### Write DESeq normalization function (using DESeq2's size factor estimation)
DESeq_norm <- function(m, sizefactors=NULL){
  
  # create counts from intensity data in the required range
  m_copy <- m
  m_copy[is.na(m_copy)] <- 0
  m_counts <- round(log2(m_copy+1)*1000,digit=0)
  library(DESeq2)
  
  
  # if no sizefactors are supplied, calculate them based on m
  if (is.null(sizefactors)){
    # create an object summarized experiment class
    dds <- DESeqDataSetFromMatrix(countData = m_counts,
                                  colData = data.frame(condition=rep("group",times=ncol(m_counts))),
                                  design =  ~ 1)
    
    # calculate normalization factors via DESeq's estimateSizeFactors(). Save them in working directory (so they can be used later on a different table of the same experiment) 
    sizefactors <- estimateSizeFactors(dds)$sizeFactor
    if (!file.exists("Results")){
      dir.create("Results")
    }
    save(sizefactors, file=paste0(getwd(),"/Results/sizefactors.Rdata"))
  } else {
    sizefactors=sizefactors
  }
  
  # perform normalization by column-wise multiplication with size-factors
  m_counts_norm  <- sweep(m_counts, STATS=1/sizefactors, FUN="*", MARGIN = 2)
  
  # retransform to intensity range 
  m_norm <- 2^(m_counts_norm/1000) - 1
  
  # return normalized intensity matrix
  return(m_norm)
}

