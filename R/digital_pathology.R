#' @title Match Feature with label in MIL situation in digital_pathology
#' @description It can directly match features with label in MIL situation in digital_pathology save the feature_lab files in "save file" for each patient.For each patient has one feature file of ten thousands tile samples and one label. The id of patient is present in file name.
#' @param feature_path the path of file which contains feature csv files. 
#' @param Label a label data frame with column names 'id','label','CancerType'.
#' @param cancer_types a vector contains all cancer types you want to focus.
#' @param save_path the file path you want to save the matched 'fealabel' data frame.
#' @param id_len one-dimension numeric.the length of the patient id in feature file name. The default value is 12 for TCGA patients.
#' @param ifclass a bool value for whether attach a 0-1 label to feature,TRUE for yes,FALSE for no. The default value is TRUE.
#' @param threshold one-dimension numeric.the threshold for 0 or 1, when is.class==TRUE. The default value is 0.4 for MSI
#' @return a list of patient id after match
#' @examples 
#' \dontrun{
#'     feature_path ="/data/cenmin/TCGA_CRC_Tumorfeaposinocsv"
#'     label_path = "/data/cenmin/statR/MSI.csv"
#'     Label = read.csv(label_path,sep="\t")
#'     save_path = "/data/cenmin/statR/TCGA"
#'     cancer_types = c("COAD","READ")
#'     colnames(Label) <- c('id','label','CancerType')
#'     save_id = matchFeaLab(feature_path,Label,cancer_types,save_path)
#' }
#' @export
matchFeaLab  <- function(feature_path,Label,cancer_types,save_path,id_len=12,ifclass = TRUE,threshold =0.4){
  id_all <- c()
  
  # feature file_name
  file_name = list.files(feature_path)
  # number of patient
  n = length(file_name)
  # feature file path
  dir = paste(feature_path,"/",file_name,sep="")
  
  # get people id name
  id_name = c() 
  for(i in 1:length(file_name)){
    id_name = c(id_name,substring(file_name[[i]],1,id_len))
  }
  
  # get label:id,MSI,Type
  Label_need = Label[Label$CancerType %in% cancer_types,]
  Label_need$id <- substring(Label_need$id,1,id_len)
  
  # match
  if(ifclass == TRUE){ # 0-1 label
    for(i in 1:n){
      if(id_name[i] %in% Label_need$id){
        id_all <- append(id_all,id_name[i])
        #print(file_name[i])
        fea <- read.csv(dir[i],header = FALSE)
        la <- Label_need[Label_need$id==id_name[i],'label']
        lab <- rep(ifelse(la>=threshold,1,0),nrow(fea))
        fealab <- cbind(fea,lab)
        write.csv(fealab,file = paste0(save_path,"/",file_name[i]),row.names = FALSE)
      }
    }
  }else{ # continuous label
    for(i in 1:n){
      if(id_name[i] %in% Label_need$id){
        id_all <- append(id_all,id_name[i])
        #print(file_name[i])
        fea <- read.csv(dir[i],header = FALSE)
        la <- Label_need[Label_need$id==id_name[i],'label']
        lab <- rep(la,nrow(fea))
        fealab <- cbind(fea,lab)
        write.csv(fealab,file = paste0(save_path,"/",file_name[i]),row.names = FALSE)
      }
    }
  }
  return(id_all)
}


#' @title train models for continuous label
#' @description this function can automatically split the train test datasets, and train lm lasso and svr models for continuous label and gain the MSE on test set and models. Also plot the prediction vs true y of test data. 
#' @param data a data frame with feature and label, the column name of  label should be set at the last column.
#' @param frac a one-dimension numeric data between (0,1). The proportion of train_data.The default value is 0.8.
#' @param seed the seed of split dataset into train and test data, the default value is 12345.
#' @param models a vector of model names want to be trained.the default value is c("lm","lasso","svr"),you can input any subset of default value expect null.
#' @param kernel the kernel in svm model.choose one of "linear","polynomial","radial basis","sigmoid".
#' @return a list of MSE and models.
#' @import e1071
#' @import caret
#' @import glmnet
#' @examples 
#' \dontrun{
#' #create dataset
#' set.seed(12)
#' x1 <- 1:100
#' x2 <- seq(2,18,length.out = 100)+rnorm(100,3,4)
#' y <- 4x1 - 3x2+ rnorm(100,0,2)
#' data <- data.frame(x1=x1,x2=x2,y=y)
#' results <- fitcontimodels(data)
#' }
#' @export
fitcontimodels <- function(data,frac = 0.8,seed=12345,models = c("lm","lasso","svr"),kernel="linear"){
  # split train-test data
  y = 0 #initial for R package
  n <- nrow(data)
  m <- ncol(data)
  colnames(data) <- c(paste0('v',1:(m-1)),'y')
  set.seed(seed)
  train_index <- sample(c(1:n),floor(frac*n))
  train_data <- data[train_index,]
  test_data <- data[-train_index,]
  
  
  # define MSE
  MSE <- function(pred,true){
    mean((pred-true)^2)
  }
  MSEs <- c()
  row_names <- c()
  model = list()
  
  # lm
  if("lm" %in% models){
    model_lm <- lm(y~., data = train_data)
    pred_lm <- predict(model_lm,test_data)
    MSE_lm <- MSE(pred_lm,test_data$y)
    plot(test_data$y,pred_lm,main="LM",xlab="True_y in test_data",ylab = "Pred_y in test_data")
    abline(0,1,col="red")
    MSEs <- c(MSEs,MSE_lm)
    row_names <- c(row_names,"lm")
    model1 <- list(model_lm)
    model <- append(model,model1)
  }
  # lasso
  if("lasso" %in% models){
    cv_lasso <- cv.glmnet(as.matrix(subset(train_data,select=-y)),train_data$y,family="gaussian",alpha=1,type.measure = "mse",grouped=FALSE)
    model_lasso <- glmnet(as.matrix(subset(train_data,select=-y)),train_data$y,ffamily = "gaussian", alpha = 1)
    pred_lasso <- predict(model_lasso, newx = as.matrix(subset(test_data,select=-y)), s=cv_lasso$lambda.min, type="response")
    MSE_lasso <- MSE(pred_lasso,test_data$y)
    plot(pred_lasso,test_data$y,main="LASSO",xlab="True_y in test_data",ylab = "Pred_y in test_data")
    abline(0,1,col="red")
    MSEs <- c(MSEs,MSE_lasso)
    row_names <- c(row_names,"lasso")
    model2 <- list(model_lasso)
    model <- append(model,model2)
  }
  # svr
  if("svr" %in% models){
    model_svr <- svm(y~., train_data,type = "eps-regression",kernel = kernel)
    pred_svr <- predict(model_svr,test_data)
    MSE_svr <- MSE(pred_svr,test_data$y)
    plot(pred_svr,test_data$y,main="SVR",xlab="True_y in test_data",ylab = "Pred_y in test_data")
    abline(0,1,col="red")
    MSEs <- c(MSEs,MSE_svr)
    row_names <- c(row_names,"svr")
    model3 <- list(model_svr)
    model <- append(model,model3)
  }
  
  # gbm and so on
  MSES <- data.frame(model=row_names,MSE=MSEs)
  return(list(MSES,model))
}