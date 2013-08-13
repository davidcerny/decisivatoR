isdecisive <- function(filename, unrooted=T, fflag=F,format="csv") {
  
  if(format == "nexus") {
    #print(paste("Reading input file ",filename,sep=''))
    print("Reading input file...")
    data <- read_nexus(filename)
    #======Taxa==========
    X <- rownames(data)
    #======Partitions of data set======================
    S <- lapply(seq_len(ncol(data)), function(i) data[,i])
    S <- lapply( seq_len(length(S)), function(i) lapply( seq_len(length(S[[i]])), function(j) if ( S[[i]][j] == 1 ) { S[[i]][j] <- X[j] } else { S[[i]][j] <- "NA"  } ) )
    S <- lapply( seq_len(length(S)), function(i) S[[i]] <- S[[i]][which(S[[i]] != "NA")]  )
    S <- lapply( seq_len(length(S)), function(i) unlist(S[[i]]))
    #print(paste("Done reading input file ",filename,sep=''))
    print("Done reading input file.")
  }else if(format == "csv") {
    #print(paste("Reading input file ",filename,sep=''))
    print("Reading input file...")
    #================Loading data==========================#
    data <- read.csv(filename,header=T)
    colnames(data)[1] <-"T"
    #======Taxa==========
    X <- as.character(lapply(seq_len(dim(data)[1]), function(i) toString(data[i,1])))
    #======Partitions of data set======================
    ans <- data.frame(data[,2:length(data)])
    S <- lapply(seq_len(ncol(ans)), function(i) ans[,i])
    S <- lapply( seq_len(length(S)), function(i) lapply( seq_len(length(S[[i]])), function(j) if ( S[[i]][j] == 1 ) { S[[i]][j] <- X[j] } else { S[[i]][j] <- "NA"  } ) )
    S <- lapply( seq_len(length(S)), function(i) S[[i]] <- S[[i]][which(S[[i]] != "NA")]  )
    S <- lapply( seq_len(length(S)), function(i) unlist(S[[i]]))
    #print(paste("Done reading input file ",filename,sep=''))
    print("Done reading input file.")
  } else {
    print(paste("Unknown format option: ",format,sep=''))
    return
  }
  
  
  #============Functions to compute decisiveness and fix the data set==================
  if(unrooted==F) {
    print("Rooted tree choice...")
    ans <- .Call("IsDecisiveRooted",X, S, length(X),length(S),fflag)
  } else {
    print("Unrooted tree choice...")
    ans <- .Call("IsDecisiveUnrooted",X, S, length(X),length(S),fflag)
  }
  #============Return value: fixed or not fixed matrix==============
  final <- matrix(0,ncol=length(ans)-1,nrow=length(X))
  if(format == "nexus") {
    rownames(final) <- rownames(data)
    colnames(final) <- colnames(data)
  } else if(format == "csv") {
    rownames(final) <- c(levels(data[,1]))
    colnames(final) <- colnames(data)[2:length(colnames(data))]
  }
  final<-as.table(final)
  for(i in seq_len(length(ans)-1)) {
    for(j in seq_len(length(ans[[i]]))) {final[ ans[[i]][j], colnames(final)[i] ] <- 1}
  }
  
  if(fflag == TRUE) {
    #Add a '*' and compute the list of suggested characters in form of {taxon, character}:
    suggested <- c()
    for(i in 1:dim(final)[1]) {
  	  for(j in 1:dim(final)[2]) {
  		  if(final[i,j] != data[i,j+1]) {
          final[i,j] = '*'
          pair <- c(rownames(final)[i],colnames(final)[j])
          suggested <-rbind(suggested, pair)
  		  }
  	  }
    }
    colnames(suggested) <- c("T","C")
    rownames(suggested) <- seq_len(dim(suggested)[1])
    out <- list(final,ans[[length(ans)]],suggested)
  } else {
    out <- list(final,ans[[length(ans)]])
  }
}

read_nexus <- function(filename) {
  library(ape)
  
  ans<-read.nexus.data(filename)
  
  t<-read.delim(filename,sep='\n')
  genes <- matrix(nrow=1, ncol=3)
  for(i in seq(1,dim(t)[1])) {
    if(length(grep('CHARSET', t$X.NEXUS[i])) > 0)  {
      if( (length(grep('coding', t$X.NEXUS[i])) > 0) || (length(grep('noncoding', t$X.NEXUS[i])) > 0) )  next
      
      trw <- strsplit(toString(t$X.NEXUS[i]),' ')
      pos <- strsplit(trw[[1]][4],'-')
      genes<-rbind( genes, c( trw[[1]][2],pos[[1]][1],substring(pos[[1]][2],1,nchar(pos[[1]][2])-1) ) )
      
    }
  }
  
  genes <- data.frame(genes[(2:dim(genes)[1]),])
  
  data_matrix <- matrix(0,nrow=length(ans),ncol=dim(genes)[1])
  colnames(data_matrix) <- genes[,1]
  rownames(data_matrix) <- names(ans)
  
  for(i in seq(1,length(ans))) {
    for(j in seq(1,dim(genes)[1])) {
      str <- paste(ans[i][as.numeric(as.character(genes$X2[j])):as.numeric(as.character(genes$X3[j]))],collapse='')
      aaa<-gregexpr('-',str)
      if( (length(aaa[[1]])/(as.numeric(as.character(genes$X3[j])) - as.numeric(as.character(genes$X2[j])))) <= 0.1) {
        data_matrix[i,j] <- 1
      }
    }
  }
  
  out <- data_matrix
}

#==============Other functions==================
