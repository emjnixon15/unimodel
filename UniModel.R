library(Rcpp)
library(tidyverse)
library(gplots)
source("plotfunctions.R")

starttime=Sys.time()


params <- list()
params <- within(params, {

  SIMTIME=300
  iperiod=5.0;
  load("data/betamat.Rdata")#the transmission matrix
  load("data/Npop.Rdata") #the  number of students in each group j (school/year)
  years=read.csv("data/yearnumbers1.csv") #a row for each group j saying which year they are - e.g row 5-32 are groups for all schools by year 1 
  nages=length(Npop) #number of groups
  numtests = read.csv("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/unimodel_work/testsperday_extractLateralFlowtesting2020.12.149am.csv", header = TRUE)  #test data - real uni booking data from the end of last term - need to decide when in the term you are going to do these tests- can we assume uptake will be similar?
  numtests = numtests[-1,]#remove na row
  
  
  gam_p=1/2;
  beta=betamat/(iperiod+1/gam_p)
  gamma=1/iperiod;  sigma=1/(3.2);  gam_a=1/(iperiod+1/gam_p); gam_h=1/3; gam_q=1/14 ##recovery rates
  testrate=1/2; testrate_a=rep(0,nages); testrate_a2=0; testrate_a3=0; ##testing parameters
  eps_h_2=c(rep(0.2,7), rep(0.4,7), rep(0.6,7), rep(0.8,7), rep(1,SIMTIME-(4*7))); #a vector of values for eps_h, one for each timestep. For now I've just said 20% will be there in the first week of term, 40% in the seconf etc., until at 5 weeks all students are back
  eps_testrate_2= c(numtests$Total.Swipes, rep(0, SIMTIME-nrow(numtests))); #a vector of values for eps_testrate, one for each timestep. Based on the booking data
  eps=0.3; eps_q=0.5; eps_h = 0; eps_testrate = 0; ##added in eps_h to scale the FOI according to the proportion of students back from holiday, eps_h is 0 when you are specifying it to be different at different times using eps_h_2, eps_testrate is for scaling the testing based on the proportion of students being tested
  h=0.002 #hospitalisation rate
  f=0.75 #fraction asymptomatic
  mrate=0.038 #mortality rate
  backgroundrate = 1e-4 #background infection rate

  ## initial conditions, list within list
  init <- within(list(), {
    S=Npop
    E=rep(0,nages)
    E[sample(1:nages,10)]=1
    A=rep(0,nages)
    A[sample(1:nages,10)]=1
    I=rep(0,nages)
    I[sample(1:nages,10)]=1
    H=rep(0,nages)
    R=rep(0,nages)
    R[sample(1:nages,10)]=1
    P=rep(0,nages)
    D=rep(0,nages)
    Hinc=rep(0,nages)
    Q=rep(0,nages)
    QA=rep(0,nages)
    N=S+E+A+I+H+R+P+Q+QA
  })
})
attach(params)
set.seed(params$seed)
sourceCpp('studentmodel_reactive2.cpp')

nsim=10
par(mar=c(5,5,1,1),mfrow=c(1,2))

result.rep=c()
for(nn in 1:nsim)
{
  set.seed(nn)
  result <- c19uni(params)
  result$nsim <- nn
  result$time[params$SIMTIME]=(params$SIMTIME-1)
  result.rep=rbind(result.rep,result)
}


Infecteds=getinfecteds(result.rep,years)
allinf=getasymp(result.rep,years)

Infecteds %>%
  group_by(time) %>%
  summarise_all(mean) -> meanout

maxbaseline=max(allinf$total)*1.2
plotcases(allinf,Infecteds,maxbaseline,'baseline')
plotyears(meanout)


end.time = Sys.time()

print(end.time-starttime)
