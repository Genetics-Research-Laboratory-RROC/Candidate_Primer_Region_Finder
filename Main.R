#############################################################################
gc()
rm(list = ls())

#############################################################################
if (! "tcltk" %in% rownames(installed.packages()))  install.packages("tcltk")

#############################################################################
OBorOT <- "OT"
DataDir <- paste(tcltk::tk_choose.dir(default = getwd(), 
                                      caption = "Select Data Directory"),"/",sep = "")
FuncDir <- MainDir <- paste(tcltk::tk_choose.dir(default = getwd(), 
                                                 caption = "Select Code Directory"),"/",sep = "")
source(paste(FuncDir, "my_intersect.R", sep = ""))
source(paste(FuncDir, "update_list.R", sep = ""))

CancerDataName <- list.files(paste(DataDir,"/",OBorOT,"/Cancer/",sep = ""))
NOfCancer <- length(CancerDataName)
NormalDataName <- list.files(paste(DataDir,"/",OBorOT,"/Normal/",sep = ""))
NOfNormal <- length(NormalDataName)
#############################################################################
NOfChr <- length(list.files(paste(DataDir,"/",OBorOT,"/Cancer/",CancerDataName[1],"/",sep = "")))

limit_len_small = 80 # minimum region length with one CpG
limit_len_big = 150  # maximum region length with one CpG
CovThr = 19
# mr_lower = 40
# mr_upper = 60
# per = 0.75

#############################################################################
Coverage = list()
MR = list()
Desired_parts_idx = list()
All_Chr = matrix(, ncol = 32, nrow = 0)

for (ChrIdx in 1:NOfChr) {
    cat("\n\nChromosome", ChrIdx, "processing start ...")
    for (IdxFile in 1:(NOfCancer+NOfNormal)) {
        cat("\n    Sample", IdxFile,"th processing ...")
        File = ifelse (IdxFile <= NOfCancer,
                       paste(DataDir,OBorOT,"/Cancer/",CancerDataName[IdxFile],"/",CancerDataName[IdxFile],"_CPGInfo",OBorOT,"_Chr",toString(ChrIdx),".txt",sep=""),
                       paste(DataDir,OBorOT,"/Normal/",NormalDataName[IdxFile-NOfCancer],"/",NormalDataName[IdxFile-NOfCancer],"_CPGInfo",OBorOT,"_Chr",toString(ChrIdx),".txt",sep=""))
        tmpData = read.table(File, header = TRUE)
        
        tmpData <- tmpData[order(tmpData$SP),]   # added (Assume data ordered by start position)
        
        MyData = t(tmpData$SP)
        LenData = length(MyData)
        Cord_len <- abs(c(MyData[1] - MyData[2], 
                          MyData[1:(LenData-2)]-MyData[3:LenData], 
                          MyData[length(MyData)-1]-MyData[length(MyData)]))  # Example: x,y,z,w  ->   abs(x-y),abs(x-z),abs(y-w),abs(z-w)
        Desired_parts_idx[[IdxFile]] = which(Cord_len>limit_len_small & Cord_len<limit_len_big)
        if (IdxFile == 1){  # for first sample
            Desired_parts = MyData[Desired_parts_idx[[1]]]
            Coverage[[1]] = tmpData$Coverage[Desired_parts_idx[[1]]]
            MR[[1]] = tmpData$MR[Desired_parts_idx[[1]]]
        } else{
            temp = my_intersect(Desired_parts, MyData[Desired_parts_idx[[IdxFile]]])
            Desired_parts = unlist(temp[1])
            kx = unlist(temp[2])
            ky = unlist(temp[3])
            
            tmp = Desired_parts_idx[1:IdxFile-1];   tmp = lapply(tmp, UpdateList, kx);   Desired_parts_idx[1:IdxFile-1] = tmp;    # update previous processed sample
            Desired_parts_idx[[IdxFile]] = Desired_parts_idx[[IdxFile]][ky]   # add current sample
            
            tmp = Coverage[1:IdxFile-1];   tmp = lapply(tmp, UpdateList, kx);   Coverage[1:IdxFile-1] = tmp;    # pdate previous processed sample
            Coverage[[IdxFile]] = tmpData$Coverage[Desired_parts_idx[[IdxFile]]]   # add current sample
            
            tmp = MR[1:IdxFile-1];   tmp = lapply(tmp, UpdateList, kx);   MR[1:IdxFile-1] = tmp;  # pdate previous processed sample
            MR[[IdxFile]] = tmpData$MR[Desired_parts_idx[[IdxFile]]]   # add current sample
        }
    }
    if(length(Desired_parts)==0){
        next
    }
    coverage_mat = matrix(unlist(Coverage),byrow = FALSE,ncol = NOfCancer+NOfNormal)
    MR_mat = matrix(unlist(MR),byrow = FALSE,ncol = NOfCancer+NOfNormal)
    MC_mat = matrix(unlist(Desired_parts), byrow = FALSE,ncol = 1)
    
    Chr = cbind(MC_mat,MR_mat,coverage_mat)
    remove_idx = which(coverage_mat < CovThr, arr.ind = TRUE)
    remove_idx = unique(remove_idx[,1])
    Chr = Chr[-remove_idx,,drop=FALSE]
    
    if (dim(Chr)[1] == 0){
        next
    }
    # nn = round(NOfNormal * per)
    # nc = round(NOfCancer * per)
    # select_idx = which((rowSums(Chr[,2:(1+NOfCancer)]>mr_upper) >= nc & rowSums(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal)]<mr_lower) >= nn) 
    #                    | (rowSums(Chr[,2:(1+NOfCancer)]<mr_lower) >= nc & rowSums(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal)]>mr_upper) >= nn))
    # Chr = Chr[select_idx,]
    
    mean_c1 = rowMeans(Chr[,2:(1+NOfCancer),drop=FALSE])
    mean_c2 = rowMeans(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal),drop=FALSE])
    tmpSd1 = lapply(data.frame(t(Chr[,2:(1+NOfCancer),drop=FALSE])),sd)
    tmpSd1 = matrix(unlist(tmpSd1), byrow = FALSE, ncol = 1)
    tmpSd2 = lapply(data.frame(t(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal),drop=FALSE])),sd)
    tmpSd2 = matrix(unlist(tmpSd2), byrow = FALSE, ncol = 1)
    fr = (mean_c1 - mean_c2)^2 / (tmpSd1^2 + tmpSd2^2)
    
    max_cancer = apply(Chr[,2:(1+NOfCancer),drop=FALSE], 1, max)
    min_cancer = apply(Chr[,2:(1+NOfCancer),drop=FALSE], 1, min)
    max_normal = apply(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal),drop=FALSE], 1, max)
    min_normal = apply(Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal),drop=FALSE], 1, min)
    
    min_overlap = apply(cbind((max_cancer - min_normal), (max_normal - min_cancer)),1, min)
    
    
    Chr = cbind(Chr[,1:(1+NOfCancer),drop=FALSE], cbind(max_cancer, min_cancer), 
                  Chr[,(2+NOfCancer):(1+NOfCancer+NOfNormal),drop=FALSE], cbind(max_normal, min_normal, min_overlap), 
                  fr,Chr[,(2+NOfCancer+NOfNormal):(dim(Chr)[2]),drop=FALSE])
    Chr = Chr[order(Chr[,(NOfCancer+NOfNormal+7)], decreasing = TRUE),,drop=FALSE]
    
    colnames(Chr) = c("SP",paste("c",sapply(1:NOfCancer, toString),"_mc",sep = ""),"max_ca","min_ca",
                         paste("n",sapply(1:NOfNormal, toString),"_mc",sep = ""),
                         "max_nc","min_nc","overlap","fisher",
                         paste("c",sapply(1:NOfCancer, toString),"_cov",sep = ""),
                         paste("n",sapply(1:NOfCancer, toString),"_cov",sep = ""))
    
    Chr = cbind("Chr" = (numeric(length(Chr[,1])) + ChrIdx),  Chr)
    All_Chr = rbind(All_Chr, Chr)
    cat("\n Number of cpg finded: ", dim(Chr)[1])
    cat("\nChromosome", ChrIdx, "processing finish.")
}

All_Chr = All_Chr[order(All_Chr[,(NOfCancer+NOfNormal+8)], decreasing = TRUE),]
write.table(All_Chr, file = paste(FuncDir,"result.txt",sep = ""), sep = "\t", row.names = FALSE)




