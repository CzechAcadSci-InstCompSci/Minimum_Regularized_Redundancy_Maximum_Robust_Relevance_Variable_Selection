### MRMR dimensionality reduction
### relevance and redundancy measures are considered here
### the whole approach performs the variable selection based on the particular choices of relevance and redundancy measures
### the classification performance, evaluated in a cross validation, depends on a single parameter denoted as beta
### such beta should be chosen (and this is implemented here to be performed automatically) which optimizes the classification performance
rm(list=ls())
### iris data; row=sample, column=variable
  # data=iris[1:100,];
  # y=(data[,5]=="setosa")*1; #binary
  # data=data[,1:4]; data=t(data); #u genetickych dat totiz geny v radcich
### the user needs a high dimensional dataset
### we use here response 'y' and two classes of gene expressions denoted as 'cmp' and 'co'
  y=(group == "CMP"); 
  cmp=which(y==1);
  co=which(y==0);

robcor=function(x, y, weightsel) ### robust correlation coefficient
# 1=linear, 2=logistic, 3=adaptive
{if (weightsel==1) 
  {w=length(x):1; w=w/sum(w); optw=slts(x,y,w=w)$weights; 
  }#if
if (weightsel==2) optw=slts(x,y)$optw; #default weights in slts
if (weightsel==3) optw=twostepslws(x,y)$optw; 
  #resulting weights of the LWS
#print("hi in robcor");print(x); print(y);
library(boot);
mat=matrix(c(x,y), ncol=2);
out=corr(mat,w=optw);
#print("in robcor, optw"); print(optw);
list(out=out)
}#robcor

relev=function(x,y)########################################
#relevance between x and y measured
{#print("hi in relev");
#print("x"); print(x); print(length(x));
#print("y"); print(length(y));
  xx=discret(x)$out;
    out=info(xx,y)$out;
 # out=cor(x,y); #correlation coefficient
  #out=t.test(x~y)$statistic; #error corrected
  #out=wilcox.test(x~y)$p.value; #only this is normalized
  #out=cor(x,y,method="spearman");
  #z=cbind(x,y); out=cor.shrink(z)[1,2]; #shrinkage makes no sense 
 ### robust correlation with a specific weight selection
   # weightsel=1; # 1 =linear, 2=logistic, 3=adaptive
   # out=robcor(x,y,weightsel)$out;
###
out=abs(out); #always
#print(c("relev, out", out));
list(out=out)
} #relev

first=function()
{#select one gene with the largest relevance
rm=rep(0,p);
#print(c("p", p));
print("computing the first gene, please wait");
for (i in 1:p)
  {#print("hi in first"); print(i);
   rm[i]=relev(mydata[i,], y*1)$out;
  }#for
#plot(rm);
vybr=which.max(rm); #dim 1; index of first gene; 
#print("first, maximal relev"); print(max(rm)); print("press sth"); readline();
nevybr = rep(1,p); #length p; values 0 or 1; 1 at the beginning 
nevybr[vybr]=0; 
print("selected:"); print(vybr); 
#print("not selected:"); print(nevybr);
#print("end of first, press sth"); readline()
list(vybr=vybr, nevybr=nevybr) 
} #first

runbeta=function(beta) #gene selection for a particular beta
{f=first(); vybr=f$vybr; nevybr=f$nevybr;
for (i in 1:9) #next gene must be selected
  {print(c("computing for i, please wait", i));
   pocetnevybr = sum(nevybr);
   loss=rep(-1,p);
   for (j in 1:p)
     if (nevybr[j]>0.5)
       {#print(c("gene", j));
        pom1=mydata[j,]; pom1=as.matrix(pom1); #j-th gene
        pom3=as.matrix(mydata[vybr,]);
          if (length(vybr)>1) pom3=t(pom3); #we need columns
        pom4 = relev(pom1, y)$out;
        pom5 = redund(pom3, pom1)$out; #order is important
        ### loss function
          loss[j]= pom4-beta*pom5; #combination of relevance, redundancy
          #loss[j]= -pom5;
        #print("index of gene; loss:"); print(c(j,loss[j]));
        #print(j); 
       }#if  
   #print(c("loss functions", loss)); print("press sth"); readline();
   mygene = which.max(loss); #this gene is selected 
   vybr=c(vybr, mygene); #vector longer by 1 value; for example 3;1
   nevybr[mygene]=0;
   ### print  
     print(c("this was i", i));
     print("genes selected:"); print(podm[vybr]); 
     #print("genes not selected:"); print(nevybr); 
     print("****************************************************"); print(" ");
  }#for ############
list(vybr=vybr)
}#runbeta

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

discret=function(z)######################################
#discretization of a continuous variable
{out=(z < mean(z))*1;
list(out=out)
}#discret

info=function(x,y)########################################
#mutual information of two binary variables
{nc=length(x);
a=sum((x==0)&(y==0))/nc;
b=sum((x==0)&(y==1))/nc;
c=sum((x==1)&(y==0))/nc;
d=sum((x==1)&(y==1))/nc;
e=sum(x==0)/nc; f=sum(x==1)/nc;
g=sum(y==0)/nc; h=sum(y==1)/nc;
out1 = a*log(a/(e*g));
out2 = b*log(b/(e*h));
out3 = c*log(c/(f*g));
out4 = d*log(d/(f*h));
out = out1+out2+out3+out4;
list(out=out)
}#info

redund=function(x,y) ##################################
#redundancy between x and y measured; x is a matrix
{#print("ahoj in reduncancy");print(dim(x)); print(length(y));print("ypsilon"); print(y);
nc=dim(x)[2];
### multiple correlation coefficient
 # ryx=rep(0, nc); #row
 # for (i in 1:nc)  ryx[i]=cor(y,x[,i]);
 # rxx=cor(x);
 # out= ryx %*% solve(rxx) %*% ryx;
    #cor() computes correlations between columns of matrices
    #square of multiple coefficient coefficient
### mutual information
  yy=discret(y)$out;
  out=0;
  for (d in 1:nc)
    {xx=discret(x[,d])$out;
     out=out+info(xx,yy)$out;
    }#for
### compromise shrinkage approach 
 # if (nc<4) 
    ### average r
      # {out=0;
      # for (d in 1:nc) out=out+abs(cor(x[,d], y));
      # out=out/nc; #average
   ### multiple r
    #  {ryx=rep(0, nc); #row
    #   for (i in 1:nc)  ryx[i]=cor(y,x[,i]);
    #   rxx=cor(x);
    #   out= ryx %*% solve(rxx) %*% ryx;
   ### common: multiple shrinkage corr.coef.
   # }  else {
   #      z=cbind(y,x);
   #      r=cor.shrink(z);
   #      ryx=r[1,2:(nc+1)];
   #      rxx=r[2:(nc+1), 2:(nc+1)];
   #      out= ryx %*% solve(rxx) %*% ryx; # or sqrt
   #   #    }#else
### average/sum of correlation coefficients
 #  out=0;
 #  for (d in 1:nc)
 #     out=out+abs(cor(x[,d], y));
    #   out=out+abs(cor(x[,d], y, method="spearman"));
    #   #w=cbind(x[,d], y); print(w); print("press sth"); readline()
    # out=out+abs( cor.shrink(cbind(x[,d],y))[1,2] ); #NONSENSE, 
    #   out=out+abs( ks.test(x[,d],y)$p.value );
    #  {signtest=sum( (x[,d]-y) > 0); out=out+binom.test(signtest,n=48)$p.value; }
 # out=out/nc; #average
### common
#print(c("redund, out", out)); 
#print("readline after redund"); readline()
list(out=out)
} #redund

### dimension reduction
mean1= apply(origdata[,cmp],1,mean); mean2=apply(origdata[,co],1,mean);
podm1=((mean1/mean2)>1.03); podm2=((mean1/mean2)<0.97); podm=which(podm1+podm2==1);
mydata=origdata[podm,]; #1261*48
nriter=10;
n=dim(mydata)[2]; #n patients
p=dim(mydata)[1]; #p genes
out=matrix(-1, nrow=10, ncol=nriter); 
for (z in 1:nriter)
  {beta= (z-1)*0.1; 
   #beta=0.3;
   print(c("beta", beta));
   ff=runbeta(beta);
   out[,z]=podm[ff$vybr]; #podm: numbering in the original matrix
   dump("out", "out.R");
  }#for

############ noise contamination & dimension reduction 
nn=dim(origdata)[1]*dim(origdata)[2];
### noise1: normal
 # noise=rnorm(nn, mean=0, sd=sqrt(0.1));
### noise2: contaminated normal
 # pom1=rnorm(nn, mean=0, sd=sqrt(0.01));
 # pom2=(runif(nn)<0.15)*rnorm(nn, mean=0, sd=1);
 # noise=pom1+pom2;
### noise3: cauchy
  pom1=rnorm(nn); pom2=rnorm(nn); #ratio of two N(0,1)
  noise= 0.002 * (pom1/pom2); #c=0.002;  
  noise=matrix(noise, nrow=dim(origdata)[1], ncol=dim(origdata)[2]);
  noisedata=origdata+noise;
###
noise=matrix(noise, nrow=dim(origdata)[1], ncol=dim(origdata)[2]);
noisedata=origdata+noise;
mean1= apply(noisedata[,cmp],1,mean); mean2=apply(noisedata[,co],1,mean);
podm1=((mean1/mean2)>1.03); podm2=((mean1/mean2)<0.97); podm=which(podm1+podm2==1);
mydata=noisedata[podm,]; #1261*48
nriter=10;
n=dim(mydata)[2]; #n patients
p=dim(mydata)[1]; #p genes
out=matrix(-1, nrow=10, ncol=nriter); 
for (z in 1:nriter)
  {beta= (z-1)*0.1; 
   #beta=0.3;
   print(c("beta", beta));
   ff=runbeta(beta);
   out[,z]=podm[ff$vybr]; #podm: numbering in the original matrix
   dump("out", "out.R");
  }#for

mylda=function(xsubset, j) ###################################
#for classification on the selected set of 10 genes
 {#print(" "); print(c("j", j)); 
  radky=xsubset[,j];
  locdata=t(origdata[radky,]); #10*48
  data1=locdata[cmp,]; data2=locdata[co,];
     xmean = apply(data1, 2, mean); 
     ymean = apply(data2, 2, mean);    
  sinv=solve(var(locdata));
  count=48;
  left=rep(-1, count); klasif=rep(-1, count);
  right=0.5*(xmean-ymean) %*% sinv %*% (xmean+ymean);
   for (i in 1:count) 
       {z=locdata[i,]; 
        left[i] =  (xmean-ymean) %*% sinv %*% z;
        klasif[i] = (left[i] > right);
      }
  #print("mylda, klasif"); print(klasif);
  ff=sesp(y,klasif); se=ff$se; sp=ff$sp; print(c("j-1, se, sp")); print(c(j-1, se, sp));
}#mylda

ldastar=function(xsubset, j) ##########################
  {#print(" "); print(c("j", j)); 
  radky=xsubset[,j];
  locdata=t(origdata[radky,]); #10*48
  data1=locdata[cmp,]; data2=locdata[co,];
     xmean = apply(data1, 2, mean); 
     ymean = apply(data2, 2, mean);    
  sinv=invcov.shrink(locdata);
  count=48;
  left=rep(-1, count); klasif=rep(-1, count);
  right=0.5*(xmean-ymean) %*% sinv %*% (xmean+ymean);
   for (i in 1:count) 
       {z=locdata[i,]; 
        left[i] =  (xmean-ymean) %*% sinv %*% z;
        klasif[i] = (left[i] > right);
      }
  #print("mylda, klasif"); print(klasif);
  ff=sesp(y,klasif); se=ff$se; sp=ff$sp; print(c("j-1, se, sp")); print(c(j-1, se, sp));
}#ldastar

### autovalidation
xsubset=source("out.R")$value; for (j in 1:10)
  f=mylda(xsubset, j);
  #f=ldastar(xsubset, j);

### leave-one-out
xsubset=source("out55.R")$value; #only these genes
out=matrix(-1, nrow=48, ncol=10);
for (j in 1:10)
  {radky=xsubset[,j]; #only these genes
   locdata=t(origdata[radky,]); #48*10
   for (i in 1:48)
     {locloc=locdata[1:47,]; yloc=y[1:47];
      if (i<48) {locloc[i,]=locdata[48,]; yloc[i]=y[48];}
      cmploc=which(yloc==1); coloc=which(yloc==0); #preparing the training set
      data1=locloc[cmploc,]; data2=locloc[coloc,];
      xmean = apply(data1, 2, mean);
      ymean = apply(data2, 2, mean);
      ### variance matrix
         #sinv=invcov.shrink(locloc);
         sinv=solve(var(locloc));
      right=0.5*(xmean-ymean) %*% sinv %*% (xmean+ymean); #classif. rule
      z=locdata[i,]; 
      left=(xmean-ymean) %*% sinv %*% z;
      out[i,j] = (left > right);
      #print("output in for cycle"); print(out);
    }#for
  }#for
youden=rep(-1, 10);
for (j in 1:10) 
  {ff=sesp(y,out[,j]); print(c((j-1)/10, ff$se, ff$sp, ff$se+ff$sp-1));
   youden[j]=ff$se+ff$sp-1;
  }#for
print(c("maximal youden", max(youden))); #############

LeaveTwelve=function(xsubset)
{dvanact=c(1,1,2,2,1, 1,2,3,3,3, 2,3,4,4,5, 4,5,4,5,6, 6,5,6,6,10, 10,10,10,11,11, 7,7,8,7,7, 8,11,12,11,12, 12,12,8,9,9, 9,8,9)  
out=matrix(-1, nrow=48, ncol=10);
for (j in 1:10) for (k in 1:12)
  {radky=xsubset[,j]; # 10 genes
   mydata=t(origdata[radky,]); #48*10
   trainvec=which(dvanact != k); # 44 patients 
   locdata=mydata[trainvec,]; #44*10
   validvec=which(dvanact == k); # 4 patients
   yloc=y[trainvec]; # 44
   cmploc=which(yloc==1); coloc=which(yloc==0); # 22; preparing the training set   
   data1=locdata[cmploc,]; data2=locdata[coloc,];
   xmean = apply(data1, 2, mean);
   ymean = apply(data2, 2, mean);
   ### variance matrix
      #sinv=invcov.shrink(locdata);
      sinv=solve(var(locdata));
   right=0.5*(xmean-ymean) %*% sinv %*% (xmean+ymean); #classif. rule
   for (i in 1:4) # only 4 patients
     {z=mydata[validvec[i],];
      left= (xmean-ymean) %*% sinv %*% z;
      out[validvec[i],j] = (left> right);
      #print("output in for cycle"); print(out);
    }#for
  }#for
for (j in 1:10) 
  {ff=sesp(y,out[,j]); print(c((j-1)/10, ff$se, ff$sp));}
}#LeaveTwelve

### call LeaveTwelve
xsubset=source("out10.R")$value; #only these genes
ff=LeaveTwelve(xsubset);

