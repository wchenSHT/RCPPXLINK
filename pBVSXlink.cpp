#include <Rcpp.h> 
#include   <stdlib.h>
#include <math.h>  
using namespace Rcpp; 
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include<vector>


// [[Rcpp::export]]

int getcharge(NumericVector x, NumericVector y,int t){
  
 int charge;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      charge=y[i];
      break;
    }
  }
  
  return charge;
  
}


// [[Rcpp::export]]

double getmass(NumericVector x, NumericVector y,int t){
  
  double mass;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
     mass=y[i];
    break;
    }
  }
  
  return mass;
  
}



// [[Rcpp::export]]

NumericVector getscan(NumericVector x, NumericVector y,int t){
  
  NumericVector z;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      z.push_back(y[i]);
    }
    else
      break;
  }
  
  return wrap(z);
  
}


// [[Rcpp::export]]
DataFrame scanassociate(NumericVector scan2, NumericVector scan3, NumericVector spectrumscan,NumericVector spectrumpepmz,NumericVector spectrumcharge,NumericVector premass, NumericVector charge, NumericVector scan2unique){
  
  std::vector< int >scanindex2;
  std::vector< int > scan3_1st;
  std::vector< int > scan3_2nd;
  std::vector< double > ms2mass;
  std::vector< int > ms2charge;
  std::vector< double >ms3_1mz;
  std::vector< double >ms3_2mz;
  std::vector< int > ms3_1charge;
  std::vector< int > ms3_2charge;
  std::vector< double >residualmass;
  std::vector< int > deltamass;  
  int m=scan2unique.size();
  
  for(int i=0;i<m;i++){
    int scan2i=scan2unique[i];
    NumericVector scan3i=getscan(scan2,scan3,scan2i);
    double premz3i=getmass(scan2,premass,scan2i);
    int charge3i=getcharge(scan2,charge,scan2i);
    
    double premassi=premz3i*charge3i-charge3i*1.0078;
    int n=scan3i.size();
    if(n>1){
      for(int j=0;j<n;j++){
        for(int t=j+1;t<n;t++){
          int scan3ij=scan3i[j];
          int scan3it=scan3i[t];
          double mass3ij=getmass(spectrumscan,spectrumpepmz,scan3ij);
          double mass3it=getmass(spectrumscan,spectrumpepmz,scan3it);
          int charge3ij=getmass(spectrumscan,spectrumcharge,scan3ij);
          int charge3it=getmass(spectrumscan,spectrumcharge,scan3it);
          double Tframassi=mass3ij*charge3ij-charge3ij*1.0078+mass3it*charge3it-charge3it*1.0078;
          double abstdeltamass=fabs(premassi-Tframassi);
          double tresmass=abstdeltamass-(int(abstdeltamass+0.5));
          int absintmass=fabs(int(abstdeltamass+0.5));
          tresmass=fabs(tresmass);
          if((tresmass<0.04)&&((absintmass<4)|((absintmass<325)&&(absintmass>317)))|((absintmass<225)&&(absintmass>219))|((absintmass<101)&&(absintmass>96))){
          scanindex2.push_back(scan2i);
          ms2mass.push_back(premz3i);
          ms2charge.push_back(charge3i);
          scan3_1st.push_back(scan3i[j]);
          scan3_2nd.push_back(scan3i[t]);
          ms3_1mz.push_back(mass3ij);
          ms3_2mz.push_back(mass3it);
          ms3_1charge.push_back(charge3ij);
          ms3_2charge.push_back(charge3it);
          residualmass.push_back(tresmass);
          deltamass.push_back(absintmass);
          }
        }
      }
        
    }
    
  }
return Rcpp::DataFrame::create(Named("scan2")=scanindex2,Named("scan3_1st")=scan3_1st,Named("scan3_2nd")=scan3_2nd,
             Named("ms2mass")=ms2mass,Named("ms2charge")=ms2charge,Named("ms3_1mz")=ms3_1mz,Named("ms3_1charge")=ms3_1charge,
                   Named("ms3_2mz")=ms3_2mz,Named("ms3_2charge")=ms3_2charge, Named("residualmass")=residualmass,Named("deltamass")=deltamass);
  
  }







  /*** R
  
  setwd("D:\\Rfiles")
    library(data.table)
    library(Rcpp)
    library(inline)
    library(ggplot2)
    
  
    psms<-read.csv("Psms.csv")
    psms<-as.data.table(psms)
    
    psms<-psms[,.(AnnotatedSequence,FirstScan,RTmin,TheoMHDa,Modifications,XCorr)]
    sequenceID<-as.character(psms$AnnotatedSequence)
    sequenceID<-toupper(sequenceID)
    score<-as.numeric(psms$XCorr)
    scanindexID<-as.integer(psms$FirstScan)
    
    scanID<-read.csv("QuantData.csv")
    scanID<-as.data.table(scanID)
    setkey(scanID,MS2ScanNumber,MS3ScanNumber)
    Scan1<-scanID$MS1ScanNumber
    scan2<-scanID$MS2ScanNumber
    scan3<-scanID$MS3ScanNumber
    scanIDprecurmass<-scanID$PrecursorMass
    scanIDcharge<-scanID$PrecursorCharge
   
    scan2unique<-unique(scan2)
    
    spectrumID<-read.csv("Spectrum.csv")
    spectrumscan<-spectrumID$FirstScan
    spectrumpepmz<-spectrumID$PrecursormzDa
    spectrumcharge<-spectrumID$PrecursorCharge
   
    
   combinevar<-scanassociate(scan2,scan3,spectrumscan,spectrumpepmz,spectrumcharge,scanIDprecurmass,scanIDcharge,scan2unique) 
   
   leng<-length(combinevar$scan2)
   seq1<-rep(NA,times=leng)
   seq2<-rep(NA,times=leng)
   score1<-rep(0,times=leng)
   score2<-rep(0,times=leng)
   mod<-rep(NA,times=leng)
   for(i in 1:leng){
    scan3_1<-combinevar$scan3_1st[i]
    index1<-which(scanindexID==scan3_1)
    if(length(index1>0)){
      seq1[i]<-sequenceID[index1]
      score1[i]<-score[index1]
    }
   }
   
   for(i in 1:leng){
     scan3_2<-combinevar$scan3_2nd[i]
     index2<-which(scanindexID==scan3_2)
     if(length(index2>0)){
       seq2[i]<-sequenceID[index2]
       score2[i]<-score[index2]
     }
   }
   
  combinevar<-data.frame(seq1,seq2,combinevar,score1,score2)
  combinevar<-as.data.table(combinevar)
  setkey(combinevar,scan2)
  
 index<-which(!is.na(combinevar$seq1))
 combine2<-combinevar[index,]
 index<-which(!is.na(combine2$seq2))
 combine2<-combine2[index,]
 
 
 leng<-length(combine2$scan2)
 
 mod1<-rep(NA,times=leng)
     mod2<-rep(NA,times=leng)
     for(i in 1:leng){
       scan3_1<-combine2$scan3_1st[i]
       index1<-which(psms$FirstScan==scan3_1)
       mod1[i]<-as.character(psms$Modifications[index1])
       
       scan3_2<-combine2$scan3_2nd[i]
       index2<-which(psms$FirstScan==scan3_2)
       mod2[i]<-as.character(psms$Modifications[index2])
       
     }
 
 write.csv(combine2,file="IDlist.csv")
 
  */