#####################################################################
#####################################################################
# ************************************************************** ####
#                       ALGORITHM FOR THE                        ####
#               AUTOMATIC CALCULATION OF RAINFALL 		     ####	
#	         THRESHOLDS FOR LANDSLIDE OCCURRENCE               ####
#                                                                ####
#  AUTHOR: ELISE MONSIEURS                                       #### 
#  CONTACT: elise.monsieurs@africamuseum.be 			     #### 
#  LICENSE: this code is licensed under GPL(version 2 or later)  ####
#  CITATION: Monsieurs, E., Dewitte, O., Demoulin, A.            ####
#  A susceptibility-based rainfall threshold approach for        ####
#  landslide occurrence                   			     ####
# ************************************************************** ####
#  This script was prepared using                         	     ####
#  R Core Team (2017). R: A language and environment for         ####
#  statistical computing. R Foundation for Statistical           ####
#  Computing, Vienna,Austria.                                    ####
#  URL http://www.R-project.org/                                 ####
#  R version: R-3.4.3                                            ####
# ************************************************************** ####
#  We acknowledge Melillo, M., Brunetti, M. T., Peruccacci, S.,  ####
#  Gariano, S. L., Roccati, A., & Guzzetti, F. for their paper   ####
#  on "A tool for the automatic calculation of rainfall          #### 
#  thresholds for landslide occurrence" in: Environmental        ####
#  Modelling & Software (2018, 105, 230-243), on which this      ####
#  script is based regarding the frequentist and                 ####
#  bootstrap statistical technique.  	                       ####
# ************************************************************** ####
#  The script requires the following input data:                 ####
#  1) LS_Inventory.csv, including following fields:              ####
#	   - ID: ID landslides					           ####
#	   - Date: Date landslide occurrence			     ####
#	   - Mining: "Y" if landslide occurred in a mining area    ####
#	     "N" if landslide did not occur in a mining area       ####
#  2) LS_Susceptibility.csv, including following fields:         ####
#	   - ID: ID landslides					           ####
#	   - Susceptibility: Landslide susceptibility [0,1]	     #### 
#  3) SRE_Day_XXX.csv                                            ####
#	   - Date: dd/mm/yyyy						     ####	
#	   - DailyR: Daily rainfall (mm)                           ####
#	This script was prepared with 3-hourly satellite rainfall  ####
#	estimates (SRE) from TMPA 3B42 RT, extracted for each      ####
#	landslides with ID XXX       		                       ####
# ************************************************************** ####
# No Packages are required to be installed.                      ####
# ************************************************************** ####
#####################################################################
#####################################################################


#Set directories
current_dir<-"C:/ARS Threshold"
dir_file_name_in<-paste(current_dir,"/INPUT",sep="")
dir_file_name_out<-paste(current_dir,"/OUTPUT",sep="")
dir.create(dir_file_name_out)


#******************************************** PART ONE *********************************************
#*********Calculate Antecedent Rainfall & prepare dataset for threshold identification**************

#PART ONE A: Calculate Antecedent Rainfall**********************************************************

#Read data landslide inventory
filename<-paste(dir_file_name_in,"/LS_Inventory.csv",sep="")
LSinv<-read.csv(filename, header=TRUE)
TotalLS<-dim(LSinv)[1]

#Initiate dataframe
LS_AntR<-data.frame(matrix(ncol = 4, nrow = TotalLS*3)) 
names(LS_AntR)<-c("ID","Date","Mining","AntR")

#set a and b in weighting funchtion
aa<-(-1.2)
b<-(1.2)
D<-"6W"

#initiate row values
p<-1
d<-2
a<-3

#Start loop landslides
for (i in 1:TotalLS){

ID<-as.character(LSinv$ID[i],stringsAsFactors = FALSE)
Date_LS<-as.character(LSinv$Date[i],stringsAsFactors = FALSE)
Mining<-as.character(LSinv$Mining[i],stringsAsFactors = FALSE)

NameID<-paste("SRE_Day_",ID,".csv", sep="")

#read LS event complete Daily time series 
filename<-paste(dir_file_name_in,"/SRE/",NameID,sep="")
RFile<-read.csv(filename, header=T)
names(RFile)<-c("Date","DailyRain")
RI<-RFile$DailyRain
lRI<-length(RI) 

#Initialize new time series for calculated Antecedent rainfall values
lRW<-list()

#Start loop time series, with D = 6 weeks (42 days)
for (j in 42:lRI){
R0<-RI[j]  #first day is RI[j] with weight 1, according to the antecedent rainfall function
RW<-R0

#Calculate weighted rain 
for (k in 1:41){RW<-RW+(RI[j-k]*(exp((aa*k)/(RI[j-k])^b)))}  
lRW[j]<-RW

#End loop time series
}

#Add dates to new rainfall time series
Dates0<-format(as.Date(RFile$Date),"%d/%m/%Y")
Dates<-Dates0[42:lRI]
Date_R<-data.frame(Date=Dates, AntecR=unlist(lRW[42:lRI]),stringsAsFactors=FALSE)

#Select AR of Day landslide (D), one day Prior to D (P), and one day After D(A)
Ind_LS_D<-which(Date_R$Date==Date_LS)
Ind_LS_P<-Ind_LS_D-1
Ind_LS_A<-Ind_LS_D+1

AntecR_D<-Date_R[Ind_LS_D,2]
AntecR_P<-Date_R[Ind_LS_P,2]
AntecR_A<-Date_R[Ind_LS_A,2]

ID_D<-paste(ID,"_D",sep="")
ID_P<-paste(ID,"_P",sep="")
ID_A<-paste(ID,"_A",sep="")

Date_LS_P<-as.Date(Date_R[Ind_LS_P,1],format = "%d/%m/%Y")
Date_LS_A<-as.Date(Date_R[Ind_LS_A,1],format = "%d/%m/%Y")

#Save Date landslide and antecedent rainfall in overview table LS_AntR   
LS_AntR$ID[p]<-ID_P
LS_AntR$Date[p]<-as.character(Date_LS_P)
LS_AntR$Mining[p]<-Mining
LS_AntR$AntR[p]<-AntecR_P

LS_AntR$ID[d]<-ID_D
LS_AntR$Date[d]<-as.character(as.Date(Date_LS,format = "%d/%m/%Y"))
LS_AntR$Mining[d]<-Mining
LS_AntR$AntR[d]<-AntecR_D

LS_AntR$ID[a]<-ID_A
LS_AntR$Date[a]<-as.character(Date_LS_A)
LS_AntR$Mining[a]<-Mining
LS_AntR$AntR[a]<-AntecR_A

p<-p+3
d<-d+3
a<-a+3

#End loop landslides
}


#PART ONE B: Create metric for sampling probability***************************************************
Probab<-data.frame(matrix(0,dim(LS_AntR)[1],1)) 
names(Probab)<-"Prob"
for (j in 1:dim(LS_AntR)[1]) {
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="D") {Probab[j,1]<-0.66}
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="P") {Probab[j,1]<-0.17}
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="A") {Probab[j,1]<-0.17}
	if (LS_AntR$AntR[j]<5){Probab[j,1]<-0}
	if (LS_AntR$Mining[j]=="Y"){Probab[j,1]<-0}
}
LS_AntR$Prob<-Probab[,1]


#PART ONE C: Add landslide susceptibility data*******************************************************
#read LS Susceptibility data
filename<-paste(dir_file_name_in,"/LS_Susceptibility.csv",sep="")
SFile<-read.csv(filename, header=T)
names(SFile)<-c("ID2","Susceptibility")

#Extract original ID data from LS_AntR
ID<-list()
for (i in 1:dim(LS_AntR)[1]){
ID[i]<-substr(LS_AntR$ID[i], 1,nchar(LS_AntR$ID[i])-2)
}
ID<-unlist(ID)
LS_AntR$ID2<-ID

#Merge landslide susceptibility data with LS_AntR
LS_AntR_S_0<-merge(LS_AntR,SFile,by="ID2")
LS_AntR_S<-LS_AntR_S_0[-match("ID2",names(LS_AntR_S_0))] #Delete ID2 column


#PART ONE D: Create AR-S dataframe for probability sampling********************************************
LS_AntR_S$SampleFreq<-0
for (i in 1:dim(LS_AntR_S)[1]){
	if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="D"){LS_AntR_S$SampleFreq[i]<-4}
if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="P"){LS_AntR_S$SampleFreq[i]<-1}
if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="A"){LS_AntR_S$SampleFreq[i]<-1}
if(LS_AntR_S$Prob[i]==0){LS_AntR_S$SampleFreq[i]<-1}
}

LS_AntR_S_Prob <- LS_AntR_S[rep(row.names(LS_AntR_S), LS_AntR_S$SampleFreq), 1:length(colnames(LS_AntR_S))-1]



#PART ONE A*: Add information on parameters of the AR function*******************************************
param<-data.frame(matrix("",ncol = 1, nrow = dim(LS_AntR_S_Prob)[1]),stringsAsFactors = FALSE)
param[1,1]<-paste("a= ",aa," ; b= ",b," ; D= ", D, sep="")
names(param)<-"ParametersAntR"
LS_AntR_S_Prob<-cbind.data.frame(LS_AntR_S_Prob,param,stringsAsFactors = FALSE)

#PART ONE D: Write resulting AR-S dataframe to OUTPUT****************************************************
filename<-paste(dir_file_name_out,"/LS_AntR_S_Prob.csv",sep="") 
write.csv(LS_AntR_S_Prob, file=filename,quote=F,row.names=F)



#********************************************* PART TWO ********************************************
#*******Procedure for the reconstruction of the multiple conditions responsible of landslide********

#Clear history function
clearhistory <- function() {
  write("", file=".blank")
    loadhistory(".blank")
      unlink(".blank")
                            }

#Identify axis labels
xls<-"LS Susceptibility"
xlr<-"Antecedent Rain (mm)"

#Read Antecedent Rainfall and susceptibility data
Susc_AR<-LS_AntR_S_Prob

#Save Original database
Susc_AR_Orig<-Susc_AR

#Delete Prob=0 out of database
Susc_AR<-Susc_AR[Susc_AR$Prob>0,]

#Prepare bootstrapping
length_LS<-dim(Susc_AR)[1]
list_n_land<-as.character(Susc_AR$ID)
sample.number<-length_LS
Perc10<-ceiling(length_LS*0.1)
Perc20<-ceiling(length_LS*0.2)
number_of_boot<-5000
quant<-c(0.5,0.1,0.05)
lquant<-length(quant)
DB_value_out<-data.frame()

alpha.5.vec<-vector(mode = "numeric", length = number_of_boot)
alpha.10.vec<-vector(mode = "numeric", length = number_of_boot)
alpha.50.vec<-vector(mode = "numeric", length = number_of_boot)
beta.5.vec<-vector(mode = "numeric", length = number_of_boot)
beta.10.vec<-vector(mode = "numeric", length = number_of_boot)
beta.50.vec<-vector(mode = "numeric", length = number_of_boot)
R_Sqr_Perc5.vec<-vector(mode = "numeric", length = number_of_boot)
R_Sqr_Perc10.vec<-vector(mode = "numeric", length = number_of_boot)
R_Sqr_Perc50.vec<-vector(mode = "numeric", length = number_of_boot)

Sig_Alfa_P10<-list()
Sig_Beta_P10<-list()
Sig_Alfa_P20<-list()
Sig_Beta_P20<-list()
Sig_Alfa_P50<-list()
Sig_Beta_P50<-list()
sa1<-1
sb1<-1
sa2<-1
sb2<-1
sa50<-1
sb50<-1


    for (n_boot in 1:number_of_boot)    #*********** open for cycle on n.boot
    {
    
      clearhistory()
      cat("\f")
      print(paste("Progress: sample size = ",sample.number," points)",sep=""),quote=F)
      print(paste("Progress :",round((n_boot*100/number_of_boot),0),"%"),quote=F)

      #************************** Bootstrap *************************
      index_rif_land<-c()
      no.sample_index_landslide<-sort(sample(c(1:sample.number),sample.number,replace = TRUE, prob = NULL))    
	index_sample_id_land<-list_n_land[no.sample_index_landslide]
 
	for (sd_i in 1:sample.number) {
	tank_index<-which(Susc_AR$ID==index_sample_id_land[sd_i]) 
      ifelse(length(tank_index)==1 ,index_rif_land<-c(index_rif_land,tank_index),
      index_rif_land<-c(index_rif_land,sample(tank_index,1)))        
      }


      no.sample<-sort(index_rif_land)
           
      AR<-Susc_AR$AntR[no.sample] 
        Susc<-Susc_AR$Susceptibility[no.sample]
           no.events<-length(AR) 
      
      #The empirical data are first logtransformed
	log.AR<-log10(AR)
      log.susc<-log10(Susc)
      fit.straight.line <- lm(log.AR~ log.susc) 
	
	a.straight<-coef(fit.straight.line)[1]
      b.straight<-coef(fit.straight.line)[2]

	R_Sqr_Perc50<-summary(fit.straight.line)$r.squared
      value.straight.line<-predict(fit.straight.line)
	Sig_Alfa_P50[sa50]<- summary(fit.straight.line)$coeff[1,4]<0.05; sa50<-sa50+1
	Sig_Beta_P50[sb50]<- summary(fit.straight.line)$coeff[2,4]<0.05 ; sb50<-sb50+1

  
      #********** Difference calculation **************************
      difference<-log.AR - value.straight.line 
	#Save as a dataframe with the index
	df_dif<-data.frame(Dif=difference, Index=seq(1, length_LS, by=1))
	#Sort from small to high residual values (difference)
	df_dif_sort<-df_dif[order(df_dif$Dif, decreasing = F), ]
	#Select the 10% and 20% lowest residual values
	Dif_Per10<-df_dif_sort[1:Perc10,]
	Dif_Per20<-df_dif_sort[1:Perc20,]
	ind_perc10<-Dif_Per10$Index
	ind_perc20<-Dif_Per20$Index

	#select the AR-S values based on the indices of the 10% & 20% lowest residual values
	log.AR_Perc10<-log.AR[ind_perc10]
	log.susc_Perc10<-log.susc[ind_perc10]
	log.AR_Perc20<-log.AR[ind_perc20]
	log.susc_Perc20<-log.susc[ind_perc20]

	#Least square fit on subsetted data
	#10% lowest residuals
	fit.straight.line_Perc10 <- lm(log.AR_Perc10~ log.susc_Perc10)	
	a.straight_Perc10 <-coef(fit.straight.line_Perc10)[1]
      b.straight_Perc10 <-coef(fit.straight.line_Perc10)[2]
	R_Sqr_Perc10<-summary(fit.straight.line_Perc10)$r.squared
	Sig_Alfa_P10[sa1]<-summary(fit.straight.line_Perc10)$coeff[1,4]<0.05 ;sa1<-sa1+1
	Sig_Beta_P10[sb1]<-summary(fit.straight.line_Perc10)$coeff[2,4]<0.05 ;sb1<-sb1+1

	#20% lowest residuals
	fit.straight.line_Perc20 <- lm(log.AR_Perc20~ log.susc_Perc20)	
	a.straight_Perc20 <-coef(fit.straight.line_Perc20)[1]
      b.straight_Perc20 <-coef(fit.straight.line_Perc20)[2]
	R_Sqr_Perc20<-summary(fit.straight.line_Perc20)$r.squared
	Sig_Alfa_P20[sa2]<- summary(fit.straight.line_Perc20)$coeff[1,4]<0.05; sa2<-sa2+1
	Sig_Beta_P20[sb2]<- summary(fit.straight.line_Perc20)$coeff[2,4]<0.05 ; sb2<-sb2+1

	#Calculation of the rainfall threshold for different exceedance probabilities)
	alpha.5.vec[n_boot]<-round(10^a.straight_Perc10,1)
	alpha.10.vec[n_boot]<-round(10^a.straight_Perc20,1)
	alpha.50.vec[n_boot]<-round(10^a.straight,1)
	beta.5.vec[n_boot]<-round(b.straight_Perc10,5) 
	beta.10.vec[n_boot]<-round(b.straight_Perc20,5) 
	beta.50.vec[n_boot]<-round(b.straight,5) 
	R_Sqr_Perc5.vec[n_boot]<-R_Sqr_Perc10
	R_Sqr_Perc10.vec[n_boot]<-R_Sqr_Perc20
	R_Sqr_Perc50.vec[n_boot]<-R_Sqr_Perc50      
      
      
    }  #*********** close for cycle on n.boot
    
    #calculate mean alpha and mean beta	
    mean.beta.5<-mean(beta.5.vec, na.rm = T)
    mean.beta.10<-mean(beta.10.vec, na.rm = T)
    mean.beta.50<-mean(beta.50.vec, na.rm = T)
    mean.alpha.5<-mean(alpha.5.vec, na.rm = T)
    mean.alpha.10<-mean(alpha.10.vec, na.rm = T)
    mean.alpha.50<-mean(alpha.50.vec, na.rm = T)
    vector_alfa_mean<-c(mean.alpha.50,mean.alpha.10,mean.alpha.5)
    vector_beta_mean<-c(mean.beta.50,mean.beta.10,mean.beta.5)

    #calculate standard deviation alpha and beta	
    sigma.beta.5<-round(sd(beta.5.vec, na.rm = T),5)
    sigma.beta.10<-round(sd(beta.10.vec, na.rm = T),5)
    sigma.beta.50<-round(sd(beta.50.vec, na.rm = T),5)	
    sigma.alpha.5<-sd(alpha.5.vec, na.rm = T)
    sigma.alpha.10<-sd(alpha.10.vec, na.rm = T)
    sigma.alpha.50<-sd(alpha.50.vec, na.rm = T)
    vector_alfa_sigma<-c(sigma.alpha.50,sigma.alpha.10,sigma.alpha.5)
    vector_beta_sigma<-c(sigma.beta.50,sigma.beta.10,sigma.beta.5)

    #Calculate mean Determination coefficient
    mean.RSqr_Perc5<-mean(R_Sqr_Perc5.vec, na.rm = T)
    mean.RSqr_Perc10<-mean(R_Sqr_Perc10.vec, na.rm = T)
    mean.RSqr_Perc50<-mean(R_Sqr_Perc50.vec, na.rm = T)
    vector_RSqr_mean<-c(mean.RSqr_Perc50,mean.RSqr_Perc10,mean.RSqr_Perc5)

    #calculate standard deviation Determination coefficient
    sigma.RSqr_Perc5<-sd(R_Sqr_Perc5.vec, na.rm = T)
    sigma.RSqr_Perc10<-sd(R_Sqr_Perc10.vec, na.rm = T)
    sigma.RSqr_Perc50<-sd(R_Sqr_Perc50.vec, na.rm = T)
    vector_RSqr_sigma<-c(sigma.RSqr_Perc50,sigma.RSqr_Perc10,sigma.RSqr_Perc5)

    #beta
    DB_value_out<-rbind(DB_value_out,data.frame(sample=sample.number,variable="beta",probability=quant[c(1:lquant)],
                                                mean=vector_beta_mean[1:lquant],sigma=round(vector_beta_sigma[1:lquant],3),
                                                min=round(vector_beta_mean[1:lquant]-vector_beta_sigma[1:lquant],1),
                                                max=round(vector_beta_mean[1:lquant]+vector_beta_sigma[1:lquant],1)       ))
    #alpha
    DB_value_out<-rbind(DB_value_out,data.frame(sample=sample.number,variable="alfa",probability=quant[c(1:lquant)],
                                                mean=vector_alfa_mean[1:lquant],sigma=round(vector_alfa_sigma[1:lquant],3),
                                                min=round(vector_alfa_mean[1:lquant]-vector_alfa_sigma[1:lquant],1),
                                                max=round(vector_alfa_mean[1:lquant]+vector_alfa_sigma[1:lquant],1)       ))

    #R-square
    DB_value_out<-rbind(DB_value_out,data.frame(sample=sample.number,variable="R²",probability=quant[c(1:lquant)],
                                                mean=vector_RSqr_mean[1:lquant],sigma=round(vector_RSqr_sigma[1:lquant],3),
                                                min=c("","",""), max=c("","","")   ))

polygon<-data.frame()
Susc<-Susc_AR$Susceptibility

l_do<-max(length(Susc),max(Susc))
x_vector<-seq.int(from=min(Susc), to=max(Susc),length.out=l_do)
polygon<-data.frame(Susc=x_vector) 
md<-max(x_vector)
mmd<-min(x_vector)
mid<-length(x_vector)

#5% threshold uncertainties ; beta and alfa are 3th and 6th th row in DB_value_out respectively.
beta_new<-DB_value_out$mean[3]
beta_min<-as.numeric(DB_value_out$min[3])
beta_max<-as.numeric(DB_value_out$max[3])
 
vector_data<-matrix(nrow = mid, ncol = 6)
polygon_5<-polygon
vector_data[1:mid,1]<-as.numeric(DB_value_out$min[6])*x_vector^(beta_new)
vector_data[1:mid,2]<-as.numeric(DB_value_out$max[6])*x_vector^(beta_new)  
vector_data[1:mid,3]<-as.numeric(DB_value_out$mean[6])*x_vector^(beta_min)
vector_data[1:mid,4]<-as.numeric(DB_value_out$mean[6])*x_vector^(beta_max)

#for all susc values, ARmin and ARmax are calculated      
for( k in 1:l_do) { 
vector_data[k,5]<-min(vector_data[k,c(1:4)])
vector_data[k,6]<-max(vector_data[k,c(1:4)]) }

#For each susceptibility value, the bounding min and max uncertainty range values are given
polygon_5<-cbind(polygon_5,data.frame(min=vector_data[,5],max= vector_data[,6] )) 
        
#10% threshold uncertainties; beta and alfa are 2th and 5th th row in DB_value_out respectively.
beta_new.10<-DB_value_out$mean[2]
beta_min.10<-as.numeric(DB_value_out$min[2])
beta_max.10<-as.numeric(DB_value_out$max[2])

vector_data_10<-matrix(nrow = mid, ncol = 6)   
polygon_10<-polygon
vector_data_10[1:mid,1]<-as.numeric(DB_value_out$min[5])*x_vector^(beta_new.10)
vector_data_10[1:mid,2]<-as.numeric(DB_value_out$max[5])*x_vector^(beta_new.10)
vector_data_10[1:mid,3]<-as.numeric(DB_value_out$mean[5])*x_vector^(beta_min.10)
vector_data_10[1:mid,4]<-as.numeric(DB_value_out$mean[5])*x_vector^(beta_max.10)
      
#for all susc values, ARmin and ARmax are calculated
for( k in 1:l_do) { 
vector_data_10[k,5]<-min(vector_data_10[k,c(1:4)])
vector_data_10[k,6]<-max(vector_data_10[k,c(1:4)]) }

#For each susceptibility value, the bounding min and max uncertainty range values are given
polygon_10<-cbind(polygon_10,data.frame(min=vector_data_10[,5],max= vector_data_10[,6] )) 

#50% threshold uncertainties ; beta and alfa are 1st and 4th th row in DB_value_out respectively.
beta_new.50<-DB_value_out$mean[1]
beta_min.50<-as.numeric(DB_value_out$min[1])
beta_max.50<-as.numeric(DB_value_out$max[1])
 
vector_data_50<-matrix(nrow = mid, ncol = 6)
polygon_50<-polygon
vector_data_50[1:mid,1]<-as.numeric(DB_value_out$min[4])*x_vector^(beta_new.50)
vector_data_50[1:mid,2]<-as.numeric(DB_value_out$max[4])*x_vector^(beta_new.50)  
vector_data_50[1:mid,3]<-as.numeric(DB_value_out$mean[4])*x_vector^(beta_min.50)
vector_data_50[1:mid,4]<-as.numeric(DB_value_out$mean[4])*x_vector^(beta_max.50)

#for all susc values, ARmin and ARmax are calculated      
for( k in 1:l_do) { 
vector_data_50[k,5]<-min(vector_data_50[k,c(1:4)])
vector_data_50[k,6]<-max(vector_data_50[k,c(1:4)]) }

#For each susceptibility value, the bounding min and max uncertainty range values are given
polygon_50<-cbind(polygon_50,data.frame(min=vector_data_50[,5],max= vector_data_50[,6] )) 
   
#Write threshold parameters in OUTPUT
filename<-paste(dir_file_name_out,"/DB_value_out.csv",sep="") 
write.csv(DB_value_out, file=filename,quote=F,row.names=F)

#Signinficance parameters 
SAP10<-unlist(Sig_Alfa_P10)
SBP10<-unlist(Sig_Beta_P10)
SAP20<-unlist(Sig_Alfa_P20)
SBP20<-unlist(Sig_Beta_P20)
SAP50<-unlist(Sig_Alfa_P50)
SBP50<-unlist(Sig_Beta_P50)

CSAP10<-length(SAP10[SAP10[TRUE]])
CSBP10<-length(SBP10[SBP10[TRUE]])
CSAP20<-length(SAP20[SAP20[TRUE]])
CSBP20<-length(SBP20[SBP20[TRUE]])
CSAP50<-length(SAP50[SAP50[TRUE]])
CSBP50<-length(SBP50[SBP50[TRUE]])

P_CSAP10<-(CSAP10/n_boot)*100
P_CSBP10<-(CSBP10/n_boot)*100
P_CSAP20<-(CSAP20/n_boot)*100
P_CSBP20<-(CSBP20/n_boot)*100
P_CSAP50<-(CSAP50/n_boot)*100
P_CSBP50<-(CSBP50/n_boot)*100

#Plot thresholds - Linear *****************************************************************
#Split landslides out and in Mining areas for plotting
NoMining<-Susc_AR_Orig[Susc_AR_Orig$Mining=="N",]
NoMining_17<-NoMining[NoMining$Prob==0.17,]
NoMining_66<-NoMining[NoMining$Prob==0.66,]
NoMining_0<-NoMining[NoMining$Prob==0,]
Mining<-Susc_AR_Orig[Susc_AR_Orig$Mining=="Y",]

filename<-paste(dir_file_name_out,"/ARSuscept.pdf",sep="")  
pdf(file=filename,paper = "special", width =30,height =21,pointsize=30,colormodel="srgb",
useDingbats=F,fillOddEven=T,version="1.7")

par(mgp=c(2.2,0.8,0),mar=c(5,4,5,2))
plot(NoMining_17$Susceptibility,NoMining_17$AntR,type="p",yaxs="i",xaxs="i",xlim=c(0,1),ylim=c(1,max(Susc_AR_Orig$AntR)),xlab=xls,ylab=xlr,pch=21,col="gray7",bg="blue",cex.lab=1,cex=0.7) 
points(NoMining_66$Susceptibility,NoMining_66$AntR,type="p",yaxs="i",xaxs="i",xlim=c(0,1),ylim=c(0,max(Susc_AR_Orig$AntR)),pch=21,cex=1.2, col="gray7",bg="blue") 
points(NoMining_0$Susceptibility,NoMining_0$AntR,type="p",yaxs="i",xaxs="i",xlim=c(0,1),ylim=c(0,max(Susc_AR_Orig$AntR)),pch=21,cex=0.5, col="red",bg="grey") 
#points(Mining$Susceptibility,Mining$AntR,type="p",yaxs="i",xaxs="i",xlim=c(0,1),ylim=c(0,max(Susc_AR_Orig$AntR)),pch=21,cex=0.5, col="red",bg="grey") 

text(0.01,max(Susc_AR_Orig$AntR)-6, "Probability", lwd=2, lty=1, bty="n", col="gray7",adj=0,cex=1.15)
legend(rep(0,2),c(max(Susc_AR_Orig$AntR)-7,1),c("0.67","0.17"),pt.bg="blue",col="gray7",lwd=NA,pch=21,bty = "n",cex=1,box.lty=0,pt.cex=c(1.2,0.7),pt.lwd=1,y.intersp = 1.25)
#legend(0,max(Susc_AR_Orig$AntR)-21,"0",pt.bg="gray",col="red",lwd=NA,pch=21,bty = "n",cex=1,box.lty=0,pt.cex=c(0.5),pt.lwd=1,y.intersp = 1.25)

	#5%
	polygon(c(x_vector,rev(x_vector)),c(polygon_5[,2],rev(c(polygon_5[,3]))),col=rgb(0.1,0.9,0.1,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[6]*x_vector^(DB_value_out$mean[3]),col="darkgreen")

	#10%
      polygon(c(x_vector,rev(x_vector)),c(polygon_10[,2],rev(c(polygon_10[,3]))),col=rgb(0.9,0.1,0.1,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[5]*x_vector^(DB_value_out$mean[2]),col="darkred")

	#50%
      polygon(c(x_vector,rev(x_vector)),c(polygon_50[,2],rev(c(polygon_50[,3]))),col=rgb(0.01,0.01,0.01,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[4]*x_vector^(DB_value_out$mean[1]),col="black")
      
	#5%
	pq<-5
      sa<-round(DB_value_out$mean[6],1)
      si<-round(DB_value_out$sigma[6],1)
      sg<-round((DB_value_out$mean[3]),2)
      sib<-round(DB_value_out$sigma[3],2)	
	ri<-round(DB_value_out$mean[9],2)	
      
      mylabel = bquote("T" ~ .(format(pq, digits = 1)) ~"%:"  ~ italic(AR)== ~ "(" ~ .(format(sa, digits = 2,nsmall = 1)) ~ "±"~ .(format(si, digits = 1,nsmall = 1)) ~ ")" ~ 
                        italic(S)^( .(format(sg, digits = 3,nsmall = 1)) ~ "±" ~ .(format(sib, digits = 2,nsmall = 1)))            )

  
	 #10%
	pq2<-10
      sa2<-round(DB_value_out$mean[5],1)
      si2<-round(DB_value_out$sigma[5],1)
      sg2<-round((DB_value_out$mean[2]),2)
      sib2<-round(DB_value_out$sigma[2],2)
	ri2<-round(DB_value_out$mean[8],2)
      mylabel2 = bquote("T" ~ .(format(pq2, digits = 1)) ~"%:"  ~ italic(AR)== ~ "(" ~ .(format(sa2, digits = 2,nsmall = 1)) ~ "±"~ .(format(si2, digits = 1,nsmall = 1)) ~ ")" ~
                          italic(S)^( .(format(sg2, digits = 3,nsmall = 1)) ~ "±" ~ .(format(sib2, digits = 2,nsmall = 1)))            )

      
	 #50%
	pq3<-50
      sa3<-round(DB_value_out$mean[4],1)
      si3<-round(DB_value_out$sigma[4],1)
      sg3<-round((DB_value_out$mean[1]),2)
      sib3<-round(DB_value_out$sigma[1],2)
	ri3<-round(DB_value_out$mean[7],2)
      mylabel3 = bquote(italic(AR)== ~ "(" ~ .(format(sa3, digits = 2,nsmall = 1)) ~ "±"~ .(format(si3, digits = 1,nsmall = 1)) ~ ")" ~
                          italic(S)^( .(format(sg3, digits = 3,nsmall = 1)) ~ "±" ~ .(format(sib3, digits = 2,nsmall = 1)))            )

#Add labels
text(0.15,6,mylabel,col="darkgreen")
text(0.15,12,mylabel2,col="darkred")
text(0.15,18,mylabel3,col="black")
#Lables in margins
#mtext(mylabel,side=3, line=0 , lty=1, bty="n", col="darkgreen",adj=1,cex=1)  
#mtext(mylabel2,side=3, line=1 , lty=1, bty="n", col="darkred",adj=1,cex=1)      
#mtext(mylabel3,side=3, line=2 , lty=1, bty="n", col="black",adj=1,cex=1)

	
	#Plot parameter settings
      par(mfrow=c(1,1),mar=c(1, 1,1, 1))
      plot.new()
	text(0,0.99, "Antecedent rainfall", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.94, "Function: AR=sum((e^(a*t/r^b))*r)", lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.89, paste("a = ",aa,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
      text(0,0.84, paste("b = ",b,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.79, paste("Accumulation period=",D,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)

	text(0,0.6, "Significance at 0.05 sign. level for 5000 iteration", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.55, paste("T5 alfa significance: ",P_CSAP10,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.50, paste("T5 beta significance: ",P_CSBP10,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.45, paste("T10 alfa significance: ",P_CSAP20,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.40, paste("T10 beta significance: ",P_CSBP20,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.35, paste("T50 alfa significance: ",P_CSAP50,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.30, paste("T50 beta significance: ",P_CSBP50,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)

	text(0,0.2, "Mean determination coefficient R²", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.15, paste("R² best fit of subset 10% lowest residuals = ",ri,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.10, paste("R² best fit of subset 20% lowest residuals = ",ri2,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.05, paste("R² best fit of subset 50% lowest residuals = ",ri3,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)




 dev.off()




#Plot thresholds - Log-scale *****************************************************************
filename<-paste(dir_file_name_out,"/ARSuscept_log.pdf",sep="")  
pdf(file=filename,paper = "special", width =30,height =21,pointsize=30,colormodel="srgb",
useDingbats=F,fillOddEven=T,version="1.7")

plot(NoMining_17$Susceptibility,NoMining_17$AntR,type="p",yaxs="i",xaxs="i",log="xy",xlab=xls,ylab=xlr,pch=21,col="gray7",bg="blue",cex.lab=1,cex=0.7)
points(NoMining_66$Susceptibility,NoMining_66$AntR,type="p",yaxs="i",xaxs="i",pch=21,cex=1.2, col="gray7",bg="blue")
points(NoMining_0$Susceptibility,NoMining_0$AntR,type="p",yaxs="i",xaxs="i",pch=21,cex=0.5, col="red",bg="grey") 
#points(Mining$Susceptibility,Mining$AntR,type="p",yaxs="i",xaxs="i",pch=21,cex=0.5, col="red",bg="grey") 

text(0.38,8.1, "Probability", lwd=2, lty=1, bty="n", col="gray7",adj=0,cex=1.15)
legend(rep(0.38,2),c(8.1,1),c("0.67","0.17"),pt.bg="blue",col="gray7",lwd=NA,pch=21,bty = "n",cex=1,box.lty=0,pt.cex=c(1.2,0.7),pt.lwd=1,y.intersp = 1.25)
#legend(0.38,8,"0",pt.bg="gray",col="red",lwd=NA,pch=21,bty = "n",cex=1,box.lty=0,pt.cex=c(0.5),pt.lwd=1,y.intersp = 1.25)


	#5%
	polygon(c(x_vector,rev(x_vector)),c(polygon_5[,2],rev(c(polygon_5[,3]))),col=rgb(0.1,0.9,0.1,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[6]*x_vector^(DB_value_out$mean[3]),col="darkgreen")
           
	#10%
      polygon(c(x_vector,rev(x_vector)),c(polygon_10[,2],rev(c(polygon_10[,3]))),col=rgb(0.9,0.1,0.1,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[5]*x_vector^(DB_value_out$mean[2]),col="darkred")
     
	#50%
      polygon(c(x_vector,rev(x_vector)),c(polygon_50[,2],rev(c(polygon_50[,3]))),col=rgb(0.01,0.01,0.01,0.1),border=NA) 
      lines(x_vector,DB_value_out$mean[4]*x_vector^(DB_value_out$mean[1]),col="black")

#Add labels in margins
mtext(mylabel,side=3, line=0 , lty=1, bty="n", col="darkgreen",adj=1,cex=1)  
mtext(mylabel2,side=3, line=1 , lty=1, bty="n", col="darkred",adj=1,cex=1)      
mtext(mylabel3,side=3, line=2 , lty=1, bty="n", col="black",adj=1,cex=1)


	#Plot parameter settings
      par(mfrow=c(1,1),mar=c(1, 1,1, 1))
      plot.new()
	text(0,0.99, "Antecedent rainfall", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.94, "Function: AR=sum((e^(a*t/r^b))*r)", lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.89, paste("a = ",aa,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
      text(0,0.84, paste("b = ",b,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.79, paste("Accumulation period= ",D,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)

	text(0,0.6, "Significance at 0.05 sign. level for 5000 iteration", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.55, paste("T5 alfa significance: ",P_CSAP10,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.50, paste("T5 beta significance: ",P_CSBP10,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.45, paste("T10 alfa significance: ",P_CSAP20,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.40, paste("T10 beta significance: ",P_CSBP20,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.35, paste("T50 alfa significance: ",P_CSAP50,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.30, paste("T50 beta significance: ",P_CSBP50,"%",sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)


	text(0,0.2, "Mean determination coefficient R²", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1.5)
	text(0,0.15, paste("R² best fit of subset 10% lowest residuals = ",ri,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.10, paste("R² best fit of subset 20% lowest residuals = ",ri2,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)
	text(0,0.05, paste("R² best fit of subset 50% lowest residuals = ",ri3,sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=1)

dev.off()



#Calculate N° landslides below the %curves*************************************************************************************
AR<-Susc_AR$AntR
susc<-Susc_AR$Susceptibility

T5_AR_Predict<-sa*susc^sg
T10_AR_Predict<-sa2*susc^sg2

T5_difference<-AR-T5_AR_Predict
T10_difference<-AR-T10_AR_Predict

Below_T5<-length(T5_difference[T5_difference<0])
Below_T10<-length(T10_difference[T10_difference<0])

Below<-data.frame(Threshold=c("T5","T10"), DataPoints_Below=c(Below_T5,Below_T10))

filename<-paste(dir_file_name_out,"/BelowThresholdPoints.csv",sep="") 
write.csv(Below, file=filename,quote=F,row.names=F)



