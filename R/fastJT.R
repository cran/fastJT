fastJT <- function(Y, X, outTopN = 15L, numThreads = 1L, standardized = TRUE) {
	markerNames <- colnames(Y)
	SNPNames <- colnames(X)
	
	if(!is.numeric(X))
	{
		stop("Non-numerical matrix for X encounted!")
	}
	# provide default col and row names for the data frame if not provided.
    if (is.null(markerNames))
        markerNames <- paste("Mrk:", 1:ncol(Y), sep="")
    if (is.null(SNPNames))
        SNPNames <- paste("SNP:", 1:ncol(X), sep="")
	
	if(numThreads == 1) # calling single treaded function
	{	
	# compute the statistics
		if(is.na(outTopN))
    	{
			res <- .Call('fastJT_fastJT', PACKAGE = 'fastJT', Y, X, !is.na(outTopN), 15, standardized)
			rownames(res$J) <- SNPNames
			res$XIDs <- SNPNames
    	}
		else
		{	
			res<-.Call('fastJT_fastJT', PACKAGE = 'fastJT', Y, X, !is.na(outTopN), outTopN, standardized)
		
			# using the SNP names instead of there matrix index for the result
			temp <- (res$J)	
			rows <- min(outTopN,ncol(X))
			for(i in 1:rows)
			{		for(j in 1:ncol(Y))
					{temp[i,j] <- SNPNames[res$XIDs[i,j]+1]}
			}
			res$XIDs <-temp
			colnames(res$XIDs) <- markerNames
		}
	}
	
	if(numThreads != 1) # calling parallel version function
	{	
		# compute the statistics
		if(is.na(outTopN))
    	{
			res<-.Call('fastJT_fastJTmp', PACKAGE = 'fastJT', Y, X, !is.na(outTopN), numThreads, 15, standardized)
        	rownames(res$J) <- SNPNames
			res$XIDs <- SNPNames
		}
		else
		{	
			res<-.Call('fastJT_fastJTmp', PACKAGE = 'fastJT', Y, X, !is.na(outTopN), numThreads, outTopN, standardized)
			# using the SNP names instead of there matrix index for the result
			temp_mp <- (res$J)
			rowsmp <- min(outTopN,ncol(X))	
			for(i in 1: rowsmp)
			{	for(j in 1:ncol(Y))
					{temp_mp[i,j] <- SNPNames[res$XIDs[i,j]+1]}
			}
			res$XIDs <-temp_mp
			colnames(res$XIDs) <- markerNames
		}
	}	
	# set the col and row names for the results matrix.
	colnames(res$J) <- markerNames
	
	# set class for the result
	class(res) <- "fastJT"	
	
	# add attribute for the result
	attr(res, 'outTopN') <- outTopN
	attr(res, 'standardized') <- standardized	
	return(res)
}




