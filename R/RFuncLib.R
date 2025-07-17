List_NA = function(DM){
  ColVec = colnames(DM)

  cs = colSums(is.na(DM))

  ncs = length(cs)

  IDVec = c()
  NameVec = c()
  for(i in 1:ncs){
    val = cs[[i]][1]
    if(val!=0){
      IDVec = c(IDVec, i)
      NameVec = c(NameVec, ColVec[i])
    }
  }

  nVec = length(IDVec)

  lst = list()

  if(nVec!=0){
    for(i in 1:nVec){
      idx = IDVec[i]
      Vec = DM[, idx]
      wh = which(is.na(Vec)==TRUE)
      lst[[NameVec[i]]]=wh
    }


  }

  return(lst)

}


#### Get the upper directory
GetUpperDir = function(CurDir, option=1){

  tmpstr = strsplit(CurDir, "/")

  CurFolder = tmpstr[[1]][length(tmpstr[[1]])]
  if(option!=1){
    CurFolder=paste("/", CurFolder, sep="")
  }

  tmpstr2 = strsplit(CurDir, CurFolder)
  UpperDir = tmpstr2[[1]][1]

  return(UpperDir)

}



Check_Variable = function(DM, Required_Vars){

  DM_Vars = colnames(DM)
  n_required = length(Required_Vars)

  Common_Vars = intersect(Required_Vars, DM_Vars)

  n_common = length(Common_Vars)

  if(n_common<n_required){

    for(i in 1:n_required){

      var_name = Required_Vars[i]
      wh = which(DM_Vars==var_name)

      if(length(wh)==0){
        tmpstr = paste(var_name, " is missing.")
        print(tmpstr)
        print("")
        return(0)
      }

    }
  }else{

    DM2 = DM[, Required_Vars]

    #DM2[1, ] = rep(NA, times=ncol(DM2))

    lst = List_NA(DM2)
    len_lst = length(lst)

    if(len_lst!=0){
      print("There is NA in the dataset.")
      print("Variables and indices of NA's in the variables")
      print(lst)
      return(0)
    }

  }

  print("Good to go!")
  return(1)
}







Get_Quantile=function(Vec, alpha){

  quan = quantile(Vec, (1-alpha))
  QuanVal = quan[[1]][1]

  return(QuanVal)

}


RFInvUVec = function(xVec, N){
  n = length(xVec)
  #X_sort = sort(xVec)
  IndVec = floor(runif(N)*n)+1
  xGen = xVec[IndVec]
  return(xGen)
}



# Kernerl density probability
kden_prob = function(fkden, fxval){
  width <- fkden$x[2] - fkden$x[1]
  start <- min(fkden$x)
  x_index <- as.integer(( fxval - start) / width) + 1
  prob <- fkden$y[x_index]
  #print("width is yahooo!!!")
  #print(width)
  return(prob)
}


Compute_AUC = function(Mat){

  nLen = nrow(Mat)
  AUC = 0

  FPVec = rep(0, times=nLen)
  TPVec = rep(0, times=nLen)

  for(i in 1:nLen){

    nTP=Mat[i,1]
    nFP=Mat[i,2]
    nFN=Mat[i,3]
    nTN=Mat[i,4]

    TPVec[i] = nTP/(nTP+nFN)
    FPVec[i] = nFP/(nFP+nTN)


  }

  TPVec2 = c(0, TPVec)
  FPVec2 = c(0, FPVec)


  for(i in 1:nLen){

    SFP = FPVec2[i]
    EFP = FPVec2[i+1]

    STP = TPVec2[i]
    ETP = TPVec2[i+1]

    Inc = abs(SFP-EFP)*(STP+ETP)/2
    AUC = AUC+Inc


  }

  return(AUC)


}




