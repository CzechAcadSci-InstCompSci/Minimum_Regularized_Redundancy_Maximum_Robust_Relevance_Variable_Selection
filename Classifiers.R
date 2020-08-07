### auxiliary file comparing various classifiers suitable for high-dimensional data
### the same data as in DimReduction.R must be used
### contaminations of the raw data are considered here as well

sesp=function(yvalid, klasif)
  {sp=0; ppv=0; npv=0;
   se= sum((yvalid==1)&(klasif==1))/sum(yvalid==1);
   sm=sum(yvalid==0);
     if (sm>0) sp= sum((yvalid==0)&(klasif==0))/sm;
  sm=sum(klasif==1);
     if (sm>0) ppv= sum((yvalid==1)&(klasif==1))/sm;
  sm=sum(klasif==0);
     if (sm>0) npv= sum((yvalid==0)&(klasif==0))/sm;
#print("y");print(yvalid);print("klasif");print(klasif);print("end in sesp");
  list(se=se,sp=sp)
 }#sesp

### PAM
locdata=noisedata;
  library(pamr);
  mydata <- list(x=locdata,y=y)
  mytrain <- pamr.train(mydata) #overview for various values of the threshold
  moje= 2.234; #to be manually set to the optimal value
  print(c("we choose a particular threshold", moje));
  sm= pamr.predict(mytrain, mydata$x, threshold=moje); #prediction for the given threshold
  sm=as.integer(sm)-1;
  print(sm-y);
  sum(abs(sm-y));
  #mycv = pamr.cv(mytrain,mydata)

### SCRDA
library(rda);
#rows must be variables  
yy=y*1+1;
myrda=rda(noisedata, yy);

########### noise contamination 
nn=dim(globdata)[1]*dim(globdata)[2];
### noise1: normal
 # noise=rnorm(nn, mean=0, sd=sqrt(0.1));
### noise2: contaminated normal
  pom1=rnorm(nn, mean=0, sd=sqrt(0.01));
  pom2=(runif(nn)<0.15)*rnorm(nn, mean=0, sd=1);
  noise=pom1+pom2;
### noise3: cauchy
  pom1=rnorm(nn); pom2=rnorm(nn); #ratio of two N(0,1)
  noise= 0.002 * (pom1/pom2); #c=0.002;  
  noise=matrix(noise, nrow=dim(globdata)[1], ncol=dim(globdata)[2]);
  noisedata=globdata+noise;
### use always
  noise=matrix(noise, nrow=dim(globdata)[1], ncol=dim(globdata)[2]);
  noisedata=globdata+noise;

### PCA + LDA
library(MASS);
pcadata=t(noisedata);
count=20;
  pca=prcomp(pcadata, retx=TRUE); #observations in rows of the dataset
   #ano, vlastni cisla vysvetli 6%, 5% atd.
  xsubset=pca$x[,1:count];
  my=lda(xsubset, y);
  yhat=predict(my, xsubset)$class;
  yhat=as.integer(yhat)-1;
  ff=sesp(y,yhat); se=ff$se; sp=ff$sp; youden=se+sp-1; print(c("se,sp, youden"));
  print(c(se,sp,youden));

### lasso
library(glmnet); #nelze na serveru
tdata=t(noisedata); #48*38590
fit=glmnet(tdata,y, family="binomial");
a=predict(fit,newx=tdata[1:48,],s=0.01); #both positive and negative
sum((a>0)-y)  #number of classification errors

### SVM
library(e1071); 
tdata=t(noisedata); #48*38590
mysvm=svm(tdata, y);
a=predict(mysvm,tdata); #od 0 do 1
sum((a>0.5)*1-y) #number of classification errors

### RRLDA
library(rrlda); 
a=rrlda(noisedata,y);

