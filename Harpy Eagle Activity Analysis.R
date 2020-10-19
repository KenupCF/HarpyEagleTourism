## Loading necessary R packages
library("plyr")
library("dplyr")
library("magrittr")
library("lubridate")
library("circular")
library("StreamMetabolism")
library("RPushbullet")
library("ggplot2")
library("maptools")
library("devtools")
library("lme4")


setwd("C:/Users/caiok/Dropbox/03-Work/01-Science/01-Papers/00-Sendo Escritos/HarpyActivity")
setwd("C:/Users/ckenup/Dropbox/03-Work/01-Science/01-Papers/00-Sendo Escritos/HarpyActivity")

# Load custom functions from GIThub
source_url("https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/Functions/Independent_records%20v3.0.R")
source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")
	
# Define custom functions

na.rm<-function(x){x<-x[!is.na(x)];return(x)}
empty2na<-function(x){if(length(x)==0){x<-NA};return(x)}
inf2nan<-function(x){x[x%in%c(Inf,-Inf)]<-NaN;return(x)}
nan2zero<-function(x){x[is.nan(x)]<-0;return(x)}
na2zero<-function(x){x[is.na(x)]<-0;return(x)}
na2false<-function(x){x[is.na(x)]<-FALSE;return(x)}
cap<-function(x,max=Inf,min=-Inf){x[x>max]<-max;x[x<min]<-min;return(x)}

qsummary<-function(x,vec=c(0,.025,.1,.25,.5,.75,.95,.975,1)){
  	  y<-quantile(x,vec,na.rm=T)
	  y<-c(y,Mean=mean(x,na.rm=T))
	  return(y)
  	  }
	  
bootstrap_time<-function(y,n.simu){
  
  ####Y é o vetor de sttime
  resu<-list()
  for (i in 1:n.simu){
    run=FALSE
    while( run == FALSE ) {  #gambiarra para reamostrar caso tenha valores que darão erro na modal.region
      resampled<-sample(y,replace=TRUE)
      run<-tryCatch(
        expr={class(modal.region(x=resampled,q=.5, bw=5))=="modal.region.circular"},
        error=function(e){F},
        finally=function(e){F})
    } #fechando a função
    resu[[i]]<-list(kernel95=modal.region(x=resampled,q=.95, bw=5))	
    resu[[i]]$kernel50<-modal.region(x=resampled,q=.5, bw=5)	
    resu[[i]]$range95<-sum(apply(resu[[i]]$kernel95$zeros,1,diff))
    resu[[i]]$range50<-sum(apply(resu[[i]]$kernel50$zeros,1,diff))
  }	#fechando o for i
  return(resu)
} #fechando a função

##Generating (smoothed) predictions from GLMMs
smoothPredMer<-function(
	  #Model object
	  mod,
	  #Wether values are scaled
	  scaled=TRUE,
	  #List with mean and SD values for each covariable
	  scaleList=NULL,
	  #Number of simulations to run
	  y.multiplier=1,
	  nsims=1e3){
	  
	  #define which parameters are present in the model
		pars<-colnames(mod@frame)[-c(1,ncol(mod@frame))] 
		par.cols<-pars
		#define which is the random-effects variable
		random.var<-colnames(mod@frame)[ncol(mod@frame)]
		
		#If there are any explanatory variables
		if(length(pars)>0){
		  
		#define the response variable 
		resp.var<-colnames(mod@frame)[1] #
		
		#define the inverse of the link function of the model
		ilink<- mod@resp$family$linkinv 
		
		
		dat<-mod@frame #extract the underlying data
		dat<-dat[!is.na(dat[,resp.var]),] #remove NA entries from data
		
		dat[,resp.var]%<>%as.numeric #convert response variables to numeric data (T/F becomes 1/0)
		
		# Generate random sample of random random effects variable
		randomEffs<-sample(x = dat[,random.var],size = 100,replace=TRUE)
		
		###Simulating data
		dat.sim<-lapply(pars,function(x){ #for each parameters
			np<-pars[pars!=x] #define which are the other parameters, to be kept constant
			
		#Create temporary data
		temp.dat<-cbind(data.frame(
		    #100 values in the range of the current variable
				a = seq(min(dat[,x]), max(dat[,x]),length = 100)), 
			# bird,
			 #100 values with the mean of the fixed variables
		  	matrix(
		 	    apply(t(t(dat[,np])),2,mean),
				  nrow=100,ncol=length(np),byrow=TRUE))
			
			#Rename colnames accordingly
			colnames(temp.dat)<-c(x,np)
			
			#Add birds data
			temp.dat[,random.var]<-randomEffs
			
			####Return predicted values from bootstrap
			mySumm <- function(.) {
        predict(., newdata=temp.dat, re.form=NULL)
			}
			
      ####Collapse bootstrap results into median, 95% PI
      sumBoot <- function(merBoot) {
        return(
          data.frame(Fitted = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))),
                     Lower = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))),
                     Upper = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))))
          )
        )
      }
			
      ##Bootstrap glmer model using simulated data
			boot1<-bootMer(mod, mySumm, nsim=nsims, use.u = FALSE, type="parametric",parallel = "multicore",ncpus = 4)
			##Summarise bootstrap
			PI.boot1<-sumBoot(boot1)
			##Added fited bootstrap data to temp.dat
			temp.dat%<>%cbind(PI.boot1)
			
			#Column names management
		  colnames(temp.dat)[1]<-'value'
			temp.dat$coef<-x
			
			#Return actual values, if the values were previously scaled
			if(scaled){
			  temp.dat$unscaled<-(temp.dat$value * scaleList$scale[x]) + scaleList$center[x]
			}else{
			temp.dat$unscaled<-temp.dat$value}
			
		return(temp.dat)		
		})
		
		#Name each data frame of simulated data
		names(dat.sim)<-pars
		
		#If values are scaled, transform them back to original values
		if(scaled){
		unsc<-sapply(par.cols,function(x){(dat[,x]*scaleList$scale[x])+scaleList$center[x]})
		colnames(unsc)<-par.cols
		dat[,par.cols]<-unsc}
		

		resu<-list()
		
		#For each set of simulated data (i.e., each covariable in the model)
		for (n in 1:length(dat.sim)){
		  
			#Base ggplot, to generate raw and smoothed plots
			gbase<-ggplot(dat,aes_string(
			  x=names(dat.sim)[n],
			  y=as.character(formula(mod))[2])) +
      scale_x_continuous(expand = c(0, 0)) + 
			geom_point(alpha=.4,size=3)
			
			
			#Raw plot, without smoothed means and intervals
			graw<-gbase +
			geom_ribbon(data = dat.sim[[n]],
			            aes(ymin = Lower*y.multiplier, 
			                ymax = Upper*y.multiplier, x = unscaled),
			                fill = 'black', alpha = .4, inherit.aes = FALSE) + 
			geom_line(data = dat.sim[[n]],color= 'black',
			          aes(y = Fitted*y.multiplier, x = unscaled))

			#Smoothed plot
			gsmooth0<-ggplot_build(gbase+
			    geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,
			          y = Fitted*y.multiplier, colour = "mean"), method = "auto",se=F)+ 
			  geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,
			          y = Lower*y.multiplier, colour = "ci"), method = "auto",se=F)+ 
			    geom_smooth(data=dat.sim[[n]],
			      aes(x = unscaled,
			          y = Upper*y.multiplier, colour = "ci"), method = "auto",se=F))
			# +

			#Creating a dataframe from the smoothed results
			smoothDF<-data.frame(x = gsmooth0$data[[2]]$x,
			            y    = gsmooth0$data[[2]]$y,
                  ymin = gsmooth0$data[[3]]$y,
                  ymax = gsmooth0$data[[4]]$y) 
			
			

			#Dataframe keeping only the explanatory and response variables
			dataTrue<-dat[,c(resp.var,names(dat.sim)[n])]
			colnames(dataTrue)<-c("X","Y")
			
			
			#Return, for each explanatory variable a list containg
			resu[[n]]<-list(plotRaw=graw, #the raw plot of the confidence intervals calculated
			                # plotSmooth=gsmooth,
			                respVar=resp.var,         #a character value expressing the response variable
			                expVar=names(dat.sim)[n], #a character value expressing the current explanatory variable
			                model=mod,                #the model used to generate predictions
			                dataSim=dat.sim[[n]],     #the simulated data for the current explanatory variable
			                dataTrue=dataTrue,        #the true, underlying data of the model
			                dataSmoothModel=smoothDF) #a data.frame containing the smoothed values for the predictions' means and confidence intervals
		}
		return(resu)}else{
		return(NULL)}
	}
	
	
# Importing and managing data

dat<-read.csv("all_data2.csv")

# dat%<>%filter(!Nest_ID=="Gordo")

dat$Date<-paste(dat$Year,
                zero_pad(dat$Month,2),
                zero_pad(dat$Day,2),sep="-")
dat$DTime<-paste(dat$Date,dat$Time,sep= " ")
dat$Camera_ID%<>%tolower
dat$DTime<-as.POSIXct(dat$DTime,format="%Y-%m-%d %I:%M:%S %p")


dat2<-dat%>%filter(!is.na(DTime))%>%
	dplyr::mutate(NestSample=Nest_ID)
# Remove dependent records (20 minute criteria)

  # dplyr::filter(Removed.from.Analysis=="Yes")%>%

dat3<-ind_record(
  sp=dat2$Nest_ID, #vetor de espécies = coluna espécie da planilha
  indiv=dat2$Nest_ID, ##vetor de individuos = coluna da planilha
  station=dat2$Nest_ID, ##vetor de cameras
  time=dat2$DTime, #vetor com o horario/data da foto, em POSIXct
  z=20 , #critério de independencia, em minutos
  var= dat2[,!colnames(dat2)%in%c('Nest_ID','DTime')]
  )  
dat3<-dat3[dat3$independent==T,]
dat3$Dec_hour<-hour(dat3$time)+
  (minute(dat3$time)/60)

# Object containing nest names
nests<-unique(dat3$station)

# Formatting data for Activity analysis
dat3$circ<-circular(dat3$Dec_hour, units="hours")

###Harpy Eagle Activity Patterns
{

categs<-list()

categs$All<-dat3
categs$Young<-dat3%>%filter(tolower(Nestling_age)=="fledgling")
categs$Adults<-dat3%>%filter(Number_of_adults>0)
categs$Adults_w_prey<-dat3%>%filter(Number_of_adults>0,!is.na(Prey_ID))

sapply(categs,nrow)

d=categs[[3]]

resuActivia0<-lapply(categs,function(d){


  ###Calculando modal.region.circular------------------------------------
  resu95  <- modal.region(d$circ, bw=5, q=0.95) #Cria o kernel com isoplet(q) de 95%
  resu50  <- modal.region(d$circ, bw=5, q=0.50) #Cria o kernel com isoplet(q) de 50%
  resu10  <- modal.region(d$circ, bw=5, q=0.10) #Cria o kernel com isoplet(q) de 10%
  resu100 <- modal.region(d$circ, bw=5, q=1) #Cria o kernel com isoplet(q) de 05%

  #calculando activity ranges--------------------------------------
  
  bootstrap<-bootstrap_time(d$circ,1)
  range50<-sum(apply(resu50$zeros,1,diff))
  range10<-sum(apply(resu10$zeros,1,diff))
  range95<-sum(apply(resu95$zeros,1,diff))
  
  #Plotting

   # par(bty="o",cex.axis=1,cex.lab=1.1,cex.main=1.4, family="sans", las=1,
			# mgp= c(2.5,0.5,0), xaxs="i",yaxs="i",
			# tcl=-0.3 )  #define os parâmetros gráficos        
		  
   ylim<-c(0,max(resu95$density$y)+.01)
   
   ylab="Kernel density probability"
   
   gg<-ggplot(NULL,aes(as.numeric(resu95$density$x), as.numeric(resu95$density$y)))+
     coord_cartesian(xlim = c(0, 24),ylim=ylim)+
     scale_x_continuous(
			  breaks=seq(from=0,to=24,by=2),
			  expand=c(0,0))+
     scale_y_continuous(expand = c(0, 0))+
     	xlab("Hours of the day") +
			ylab("Kernel Density Probability") +
     geom_polygon(aes(
       x=c(0,as.numeric(resu100$density$x),24),
       y=c(0,as.numeric(resu100$density$y),00)),fill="gray75")+
 			theme(
			  # text = element_text(size=10*cex.),
			  panel.grid.major = element_blank(), 
			  panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"))
   
 ####Plotar CORE 50%
		  ###Mexer aqui!! Melhorar a linha 
		  {
		  x <- list()
		  y <- list()
		   for (j in 1:nrow(resu50$zeros)) {
			x[[j]]=  as.numeric(resu50$density$x[
			  as.numeric(resu50$density$x<resu50$zeros[j,2]) & 
				as.numeric(resu50$density$x>resu50$zeros[j,1])])
			
			y[[j]]= as.numeric(resu50$density$y [
			  as.numeric(resu50$density$x < resu50$zeros[j,2]) & 
				as.numeric(resu50$density$x > resu50$zeros[j,1])])
			y[[j]][c(1,length(y[[j]]))] <- 0
		  color=adjustcolor("gray50",alpha=1)
		  

			# polygon(x, y, col=color,lty='dotted')
		  } #fechando o for J
 			}
      x50<-unlist(x)
      y50<-unlist(y)
		  gg<-gg+
      geom_polygon(aes(x=x50,y=y50),fill=color)
 

######Plotar Core 10%
		  {
		  x <- list()
		  y <- list()
		   for (j in 1:nrow(resu10$zeros)) {
			x[[j]]=  as.numeric(resu10$density$x[
			  as.numeric(resu10$density$x<resu10$zeros[j,2]) & 
				as.numeric(resu10$density$x>resu10$zeros[j,1])])
			
			y[[j]]= as.numeric(resu10$density$y [
			  as.numeric(resu10$density$x < resu10$zeros[j,2]) & 
				as.numeric(resu10$density$x > resu10$zeros[j,1])])
			y[[j]][c(1,length(y[[j]]))] <- 0
		  color=adjustcolor("gray25",alpha=1)
		  

			# polygon(x, y, col=color,lty='dotted')
		  } #fechando o for J
 			}
      x10<-unlist(x)
      y10<-unlist(y)
		  # gg<-gg+
      # geom_polygon(aes(x=x10,y=y10),fill=color)
		  
		     	  
    
    resu<-list(plot=gg,range=list("95"=range95,"50"=range50,"10"=range10),
	bootstrap=bootstrap,
	resu95=resu95,
	resu50=resu50,
	resu10=resu10,
	dataset=d)
    return(resu)
			
})


}

RPushbullet::pbPost("note",title="Harpy Activity Finished",
                    apikey="o.8v2VQ88kZXsMdUrk86PtaxXfbsaPOwMm")

save(resuActivia0,file="Activity Results.RData")

####Probabilidade de uma camera detectar uma harpia
{


nests2discard<-c("NB2","Gordo","NB","Valdeci","ApiacÃ¡s_1_II")  
  
###Data.frame defining different scenarios for analysis  
simulation.pars<-expand.grid(
    # time.bins=c(1,2,3,4,5),
    time.bins=c(1),
	exp.var=c("DateBin"),
	#How long after the first visit is it feasible for an harpy eagle to be detected (in days)
  # maxTimeAllocation=c(30,120),
	maxTimeAllocation=c(30,120),
	#Which age class to analyse
	age.class=c("Adult","Fledgling"))

binned.data.list<-list()
models<-list()
Plot<-list()
bootstrapList<-list()
nSims<-1e4

# False omission rate for each age class
FOR = c(Adult=.269,Fledgling=0)


i=1
resu<-list()
bootstrapped<-list()


      ####Return predicted values from bootstrap (for a )
			mySumm <- function(.) {
        predict(., newdata=new.df, re.form=~0)
			}
			
      ####Collapse bootstrap results into median, 95% PI
      sumBoot <- function(merBoot) {
        return(
          data.frame(Fitted = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))),
                     Lower = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))),
                     Upper = ilink(apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))))
          )
        )
      }

      sumBoot <- function(merBoot) {
         return(
            data.frame(mean = apply(merBoot$t, 2, function(x) mean(x)),
                       sd = apply(merBoot$t, 2, function(x) sd(x))))
    
      }
	###

simulation.pars$exp.var<-as.character(simulation.pars$exp.var)
i=3
for(i in 1:nrow(simulation.pars)){
  
  
  binLoop<-simulation.pars$time.bins[i]
  timeAllocLoop<-simulation.pars$maxTimeAllocation[i]
    
  binned.data.list[[i]]<-dat2%>%
  dplyr::filter(!NestEvent%in%nests2discard)%>%
  dplyr::group_by(Nest_ID)%>%
  mutate(DateBin=cut(as.Date(Date),breaks=paste(binLoop,"days", sep=" "),labels=FALSE))%>%
  mutate(DateSampledBin=as.numeric(as.factor(DateBin)))%>%  
  group_by(DateBin,NestEvent)%>%
  summarise(Nest_ID=unique(Nest_ID),
            Adult=sum(na.rm(Number_of_adults>0))>0,
            Fledgling=sum(na.rm(tolower(Nestling_age)=="fledgling"))>0,
            AdultOrFledgling=Adult | Fledgling,Date=min(Date),
            DateSampledBin=unique(DateSampledBin))%>%
  filter(!is.na(Adult))%>% #@bug
  ungroup()
  
 ### Filling days without records with FALSE
  
splitted<-split(binned.data.list[[i]],binned.data.list[[i]]$NestEvent)

binned.data.list[[i]]<-rbind.fill(lapply(splitted,function(df){
  fill.list<-list(
    Adult=FALSE,Fledgling=FALSE,AdultOrFledgling=FALSE,
    NestEvent=unique(df$NestEvent),
    Nest_ID=unique(df$Nest_ID))
  df2<-df%>%
    ungroup()%>%
  mutate(MaxDate=max(DateBin),MinDate=min(DateBin),
         DateBin=factor(DateBin,levels=unique(MinDate):unique(MaxDate)))%>%
  tidyr::complete(DateBin,fill=fill.list)%>%
  mutate(DateBin=as.numeric(as.character(DateBin)))
  return(df2)
  }))%>%
  #filter(DateBin<=120)%>%
  ungroup()
  
  class(binned.data.list[[i]])<-"data.frame"

  #GLMM solution
  {
  
  Plot[[i]]<-ggplot(data=binned.data.list[[i]],
	aes_string(x=simulation.pars$exp.var[i],y=simulation.pars$age.class[i]))+
    geom_jitter(alpha=.3,height=0,width = .5)+
    facet_wrap(vars(Nest_ID),scales = "free_x")
  
  
  form<-as.formula(paste(simulation.pars$age.class[i],"~",simulation.pars$exp.var[i],"+(1|Nest_ID)"))
  
  models[[i]]<-glmer(data=binned.data.list[[i]],formula = form,family="binomial")

  new.df<-expand.grid(X=1:(max(binned.data.list[[i]][,as.character(simulation.pars$exp.var[i])])*binLoop))
  
  colnames(new.df)[1]<-as.character(simulation.pars$exp.var[i])
  
  ilink<-models[[i]]@resp$family$linkinv    
      
  boot1<-bootMer(models[[i]], mySumm, nsim=nSims, use.u = F, type="parametric",parallel = "multicore",ncpus = 4)
  
  trueProb<-ilink(boot1$t)/(1-FOR[as.character(simulation.pars$age.class[i])])
  probUndect<-1-trueProb
  

  firstVisit<-round(runif(min = 1,max=ceiling(timeAllocLoop/binLoop),n =nSims))
  
  resu[[i]]<-sapply(1:length(firstVisit),function(v){
  
    timeToFirstDetect<-min(which(sapply(trueProb[v,firstVisit[v]:ncol(trueProb)],function(x){rbinom(n = 1,size=1,prob=x)})==1))*binLoop
    
    return(timeToFirstDetect)
  
  })
  quantile(resu[[i]],.95)
  }
  
  ####BOOTSTRAPPING SOLUTION
  
  bootstrapped[[i]]<-rbind.fill(lapply(1:nSims,function(z){
    
    binned.data.list[[i]][sample(nrow(binned.data.list[[i]]),replace=F),]%>%
      ungroup()%>%
      filter(!duplicated(DateBin))%>%
      arrange(DateBin)->y
    
    class(y)<-"data.frame"
    y$Exp<-y[,as.character(simulation.pars$age.class[i])]

    y<-y%>%
      filter(DateBin>=firstVisit[z])%>%
      summarise(DayFirstDetection=min(DateBin[Exp])-firstVisit[z]+1)%>%
      mutate(Iteration=z)
    return(y)
    }))%>%
    mutate(AgeClass=simulation.pars$age.class[i],maxTimeAllocation=timeAllocLoop)
  
  
  print(i)
  }

}




resu<-lapply(resu,function(x){
  x[x==Inf]<-NA
  return(x)
  })

names(resu)<-apply(simulation.pars,1,paste,collapse = "-")
names(models)<-names(resu)

save(resu,models,bootstrapped,bootstrapList,file="Detection Probability v2.RData")
RPushbullet::pbPost("note",title="Harpy Prob Detection Finished",
                    apikey="o.8v2VQ88kZXsMdUrk86PtaxXfbsaPOwMm")


rbind.fill(bootstrapped)%>%
group_by(AgeClass,maxTimeAllocation)%>%
summarise(mean=mean(DayFirstDetection),
          q50=quantile(DayFirstDetection,.5),
          q75=quantile(DayFirstDetection,c(.75)),
          q95=quantile(DayFirstDetection,c(.95)))



##################################
###Plotting these probabilities###
##################################
{

smoothPred<-lapply(models,smoothPredMer,scaled=FALSE,nsims=nSims)
	cex.<-1.2

gg<-lapply(1:length(smoothPred),function(i){
 x<-smoothPred[[i]]
 multiplier<-1/(1-FOR[as.character(simulation.pars$age.class[i])])
 
 
 ggplot() +
	    #Plot jittered actual data
	  	geom_jitter(data = x$dataTrue,
	                mapping = aes(x=Y,y=X*multiplier),
	                height = 0,shape=3,width=.4,alpha=.4,size=2*cex.)+
	    #Plot smoothed trend line
	    geom_line(data=x$dataSmoothModel,
	              mapping = aes(x=x,y=y*multiplier))+
	    #Plot smoothed confidence intervals
	    geom_ribbon(data=x$dataSmoothModel,
	                
	                mapping = aes(x=x,ymin=ymin*multiplier,ymax=ymax*multiplier),alpha=.4)
})
save(gg,"Figures.RData")

DetectProdAdultPlot<-gg[[4]]+ylab("Probability of Adult Detection")+xlab("Days ellapsed")+
     scale_x_continuous(expand = c(0, 0))+
	 theme_bw()

DetectProdJuvenPlot<-gg[[8]]+ylab("Probability of Fledgling Detection")+xlab("Days ellapsed")+
     scale_x_continuous(expand = c(0, 0))+
	 theme_bw()
	 
	 

	svg(
	   file='.\\Figure_2.svg',width=12,height=7)

grid.arrange(ggplotGrob(DetectProdAdultPlot),ggplotGrob(DetectProdJuvenPlot),ncol=2,nrow=1)

dev.off()

	 

}

save(resu,models,gg,smoothPred,file="Detection Probability v2.RData")

####################################################
###Tempo médio de permanencia de adulto na camera###
####################################################
##

oversampledCameras<-c("João_Bocão_cam_i","João_Bocão_cam_ii")
oversampledNests<-c("JoÃ£o_BocÃ£o","NB2","Gordo","Cerrado")

z<-3 #interval below which where we can assume the harpy is around the nest

dat2$format<-gsub(x=dat2$ï..Photo_ID,pattern = ".*\\.",replacement = "")

VisitTimeDF_adult<-dat2%>%
  dplyr::mutate(NC=paste(Nest_ID,Camera_ID,sep="_"))%>%
  dplyr::filter(Nest_ID%in%oversampledNests)%>%
  dplyr::filter(Number_of_adults>0)%>%
  dplyr::group_by(NC)%>%
  dplyr::arrange(DTime)%>%
  # dplyr::filter(tolower(Adult_on_photo)=='yes')%>%
  dplyr::mutate(TimeSinceLastDetection=c(0,
    abs(as.numeric(
      DTime[-1]-DTime[-length(DTime)]))))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(Nest_ID=as.character(Nest_ID))%>%
  dplyr::arrange(Nest_ID,DTime)

VisitTimeDF_fledged<-dat2%>%
  dplyr::mutate(NC=paste(Nest_ID,Camera_ID,sep="_"))%>%
  dplyr::filter(Nest_ID%in%oversampledNests)%>%
  dplyr::filter(Nestling_age=="Fledgling")%>%
  dplyr::group_by(NC)%>%
  dplyr::arrange(DTime)%>%
  # dplyr::filter(tolower(Adult_on_photo)=='yes')%>%
  dplyr::mutate(TimeSinceLastDetection=c(0,
    abs(as.numeric(
      DTime[-1]-DTime[-length(DTime)]))))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(Nest_ID=as.character(Nest_ID))%>%
  dplyr::arrange(Nest_ID,DTime)



#### For adults

detectionDuration<-0
Nest<-VisitTimeDF_adult$Nest_ID[1]
NC<-VisitTimeDF_adult$Nest_ID[1]
lastValid<-1
for(i in 2:nrow(VisitTimeDF_adult)){
# for(i in 2:8){
  
  TimeSinceLastDetection<-difftime(
      VisitTimeDF_adult$DTime[i],VisitTimeDF_adult$DTime[i-1],units = "mins")%>%
      as.numeric%>%abs
  
  if(TimeSinceLastDetection<=z & 
     (VisitTimeDF_adult$Nest_ID[i]==VisitTimeDF_adult$Nest_ID[lastValid])){

   detectionDuration[length(detectionDuration)]<-abs(as.numeric(
     difftime(VisitTimeDF_adult$DTime[lastValid],
              VisitTimeDF_adult$DTime[i],units = "mins")))
      
    
  }else{
    lastValid<-i
    detectionDuration<-c(detectionDuration,0)
    Nest<-c(Nest,VisitTimeDF_adult$Nest_ID[i])
    }
}
detection.df_adult<-data.frame(Nest,detectionDuration)%>%
  mutate(Age="Adult")

### For Fledglings
detectionDuration<-0
Nest<-VisitTimeDF_fledged$Nest_ID[1]
NC<-VisitTimeDF_fledged$Nest_ID[1]
lastValid<-1
for(i in 2:nrow(VisitTimeDF_fledged)){
# for(i in 2:8){
  
  TimeSinceLastDetection<-difftime(
      VisitTimeDF_fledged$DTime[i],VisitTimeDF_fledged$DTime[i-1],units = "mins")%>%
      as.numeric%>%abs
  
  if(TimeSinceLastDetection<=z & 
     (VisitTimeDF_fledged$Nest_ID[i]==VisitTimeDF_fledged$Nest_ID[lastValid])){

   detectionDuration[length(detectionDuration)]<-abs(as.numeric(
     difftime(VisitTimeDF_fledged$DTime[lastValid],
              VisitTimeDF_fledged$DTime[i],units = "mins")))
      
    
  }else{
    lastValid<-i
    detectionDuration<-c(detectionDuration,0)
    Nest<-c(Nest,VisitTimeDF_fledged$Nest_ID[i])
    }
}
detection.df_fledged<-data.frame(Nest,detectionDuration)%>%
  mutate(Age="Fledgling")



detection.df<-rbind(detection.df_adult,detection.df_fledged)


detection.df%>%
  group_by(Age)%>%
  summarise(Mean=mean(detectionDuration),
            NoVisits=n(),
            MeanNoInstant=mean(detectionDuration[detectionDuration>0]))


propFast<-sum(detection.df$detectionDuration==0)/nrow(detection.df)


save(detection.df,file="Harpy Time Spent on Nest.RData")
RPushbullet::pbPost("note",title="Harpy Prob Detection Finished",
                    apikey="o.8v2VQ88kZXsMdUrk86PtaxXfbsaPOwMm")



save.image(file="Full.RData")


















detectionDuration<-0
Nest<-VisitTimeDF$Nest_ID[1]
NC<-VisitTimeDF$NC[1]
lastValid<-1
TimeSinceLastDetection<-0
for(i in 2:nrow(VisitTimeDF)){
# for(i in 2:8){
  
  
  if(
     (VisitTimeDF$NC[i]==VisitTimeDF$NC[i-1])){
  TimeSinceLastDetection[i]<-difftime(
      VisitTimeDF$DTime[i],VisitTimeDF$DTime[i-1],units = "mins")%>%
      as.numeric%>%abs}else{TimeSinceLastDetection[i]<-0}
  
    NC<-c(NC,VisitTimeDF$NC[i])
    }
  
y<-data.frame(NC,TimeSinceLastDetection,VisitTimeDF$DTime)

y%>%group_by(NC)%>%
  summarise(propFast=sum(TimeSinceLastDetection>0 & TimeSinceLastDetection<=8)/n())%>%
  arrange(desc(propFast))

RPushbullet::pbPost("note",title="Harpy Time At Nest Finished")


dat3%>%select(Date,DateBin2,DateBin3,DateBin5)%>%View


categs$Adults%>%



AdultsByNest<-split(categs$Adults,categs$Adults$station)
AdultsByNest<-AdultsByNest[sapply(AdultsByNest,function(x){nrow(x)>1})]










timeToAdultDays<-lapply(AdultsByNest,function(d){
  
  timeSecs<-as.numeric(d$time)
  range<-range(timeSecs)
  simulatedStart<-runif(10e3,min=range[1],max=range[2])
  timeToAdultDays<-sapply(simulatedStart,function(x){
      y<-timeSecs - x
      y<-y[y>0]
      return(min(y))
      })/(60*60*24)
  timeToAdultDaysq95=quantile(timeToAdultDays,.95)
  densPlot<-plot(density(timeToAdultDays))
  
  return(list(v=timeToAdultDays,q95=timeToAdultDaysq95,plot=densPlot))
  
  
  })









for(i in 1:length(resu)){
  tiff(
	   file=paste('.\\Figure1_',i,'.tiff',sep=""),
	   # units='cm',
	  res=200,
	  pointsize=12,
	  width=750,height=750)
  plot(resu[[i]]$plot)
  dev.off()
  }


resu2<-lapply(names(resu),function(x){
  r<-resu[[x]]
  boots<-r$bootstrap
  r95=sapply(boots,function(x){x$range95})
  r50=sapply(boots,function(x){x$range50})
  r95summ=data.frame(mean=mean(r95),
                     median=quantile(r95,.5),
                     lcl=quantile(r95,.025),
                     ucl=quantile(r95,.975),
                     KernellCore=95)
  r50summ=data.frame(mean=mean(r50),
                     median=quantile(r50,.5),
                     lcl=quantile(r50,.025),
                     ucl=quantile(r50,.975),
                     KernellCore=50)  
  summ<-rbind.fill(r50summ,r95summ)%>%mutate(Category=x
  r$bootSumm<-summ
  return(r)
  })

temp<-lapply(resu2,function(x){x$bootSumm})%>%rbind.fill


