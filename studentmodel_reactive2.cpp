#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List c19uni(List params) {

  int j=0, J=0;
  double k=0.0, m=0.0, i=0.0;
  //printf("1\n");
  // pull out list w/in list into its own object
  List init = params["init"];
  List S0 = init["S"];
  List E0 = init["E"];
  List A0 = init["A"];
  List P0 = init["P"];
  List I0 = init["I"];
  List R0 = init["R"];
  List H0 = init["H"];
  List Q0 = init["Q"];
  List QA0 = init["QA"];
  List D0 = init["D"];
  List N0 = init["N"];

  //printf("2\n");
  // use Rcpp as() function to "cast" R vector to cpp scalar
  int nsteps = as<int>(params["SIMTIME"]);
  int nages = as<int>(params["nages"]);

  //printf("3\n");
  // initialize each state vector in its own vector
  NumericMatrix SS(nsteps,nages);
  NumericMatrix EE(nsteps,nages);
  NumericMatrix AA(nsteps,nages);
  NumericMatrix PP(nsteps,nages);
  NumericMatrix II(nsteps,nages);
  NumericMatrix HH(nsteps,nages);
  NumericMatrix RR(nsteps,nages);
  NumericMatrix QQ(nsteps,nages);
  NumericMatrix QA(nsteps,nages);
  NumericMatrix DD(nsteps,nages);
  NumericMatrix NN(nsteps,nages);

  //printf("4\n");
  for(j=0; j<nages; j++)
 {
     SS(0,j) = S0[j];
    EE(0,j) = E0[j];
  AA(0,j) = A0[j];
 PP(0,j) = P0[j];
  II(0,j) = I0[j];
  HH(0,j) = H0[j];
 RR(0,j) = R0[j];
 QQ(0,j) = Q0[j];
 QA(0,j) = QA0[j];
 DD(0,j) = D0[j];
 NN(0,j) = N0[j];
  }
  //printf("5\n");
  // fill time w/zeros
  NumericVector time(nsteps);
  //printf("5B\n");
  // pull out params for easy reading
  //[1] "init"       "mrate"      "f"          "h"          "eps"        "testrate_a" "testrate"   "gam_q"      "gam_p"      "gam_h"      "gam_a"
  //[12] "sigma"      "gamma"      "beta1"      "rzero"      "nages"      "years"      "groupnames" "Npop"       "iperiod"    "SIMTIME"    "seed"
  double sigma = params["sigma"];
  double mrate = params["mrate"];
  double f1 = params["f"];
  double eps = params["eps"];
  double eps_q = params["eps_q"];
  double eps_h = params["eps_h"]; //new holiday parameter (proportion of students back)
  double eps_testrate = params["eps_testrate"]; //new testrate parameter (proportion of students tested)
  double h = params["h"];
  double gamma = params["gamma"];
  double gam_p = params["gam_p"];
  double gam_h = params["gam_h"];
  double gam_a = params["gam_a"];
  double gam_q = params["gam_q"];
  double backgroundrate = params["backgroundrate"];

  double testrate = params["testrate"];
  double testrate_a2 = params["testrate_a2"];
  double testrate_a3 = params["testrate_a3"];
  //printf("5C\n");
  NumericVector testrate_a = params["testrate_a"];
  NumericVector eps_h_2 = params["eps_h_2"]; // Added this in, meant to be a vector with a value for eps_h at each timestep
  NumericVector eps_testrate_2 = params["eps_testrate_2"]; // Added this in, meant to be a vector with a value for eps_testrate at each timestep
  NumericMatrix beta = params["beta"];
  NumericVector testrate_a_use = testrate_a;

  //printf("6\n");
  // declare transition vectors
  NumericVector SE(nages); //set later
  NumericVector EA(nages,f1 * sigma);
  NumericVector EP(nages,(1.0 - f1) * sigma);
  NumericVector PInf(nages,gam_p);
  NumericVector AR(nages,gam_a);
  NumericVector IH(nages,h * gamma);
  NumericVector IR(nages,(1.0 - h) * gamma);
  NumericVector HR(nages,(1 - mrate) * gam_h);
  NumericVector HD(nages,mrate * gam_h);
  NumericVector IQ(nages,eps_testrate *testrate);
  NumericVector AQ(nages);//set later
  //NumericVector AQ testrate_a;//set later
  NumericVector PQ(nages);//set later
  NumericVector QR(nages,gam_q);

  double prob = 0.0;

  int flag=0,daycount=1,newcases=0,newcases_then=0,delta=0;
  double Etot=0.0,Ntot=0.0,Ntot2=0.0;
  for (int t= 0; t < (nsteps-1); t++) {
    time[t] = t;
    eps_h=eps_h_2[t];//Updating the eps_h value so that it is correct for this current time step- getting the value for this timestep - not sure if this works
    eps_testrate=eps_testrate_2[t];//Updating the eps_testrate value so that it is correct for this current time step- getting the value for this timestep - not sure if this works

    for(j = 0; j < nages; j++)
    {
      SS(t+1,j) = SS(t,j);
      EE(t+1,j) = EE(t,j);
      AA(t+1,j) = AA(t,j);
      PP(t+1,j) = PP(t,j);
      II(t+1,j) = II(t,j);
      RR(t+1,j) = RR(t,j);
      QQ(t+1,j) = QQ(t,j);
      QA(t+1,j) = QA(t,j);
      HH(t+1,j) = HH(t,j);
      DD(t+1,j) = DD(t,j);

      if(EE(t+1,j)<0){printf("set up E[%i,%i]=%f\n",t+1,j,EE(t+1,j)); break;}
    }
    Etot=0;
    for(j = 0; j < nages; j++) Etot = Etot + EE(t+1,j);
    //printf("t=%i Etot=%f\n",t+1,Etot);
    if(isnan(Etot)){break;}

    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("1 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}

    flag=0;
    //printf("t = %i\n",t);
    for(J = 0; J < nages; J++)
      {
        delta=0;
        SE[J] = 0.0;
        prob=0.0;
        for(j = 0; j < nages; j++)
        {
          if(j==J){delta=1;}
          if(NN(t,j)>0)
          {
            //prob += beta(J, j) * (eps_q*eps*QA(t,j) + eps_q*QQ(t,j) + eps*AA(t,j)+PP(t,j)+II(t,j)) / NN(t,j);
            if(eps_q==0)
            {
              if(t>=70 & t<=84)
              {
              prob += beta(J, j) * eps_h *(eps_q*delta*eps*QA(t,j) + eps_q*delta*QQ(t,j) + eps*AA(t,j)+PP(t,j)+II(t,j)) / NN(t,j); //added in eps_h(t,j) to scale foi due to students going on holiday, not sure if this was correct
              }
              else{
                prob += beta(J, j) * eps_h *(eps_q*delta*eps*QA(t,j) + eps_q*delta*QQ(t,j) + delta*eps*AA(t,j)+PP(t,j)+delta*II(t,j)) / NN(t,j); //added in eps_h(t,j) to scale foi due to students going on holiday, not sure if this was correct
              }
            }
            else{
              prob += beta(J, j) * eps_h *(eps_q*delta*eps*QA(t,j) + eps_q*delta*QQ(t,j) + eps*AA(t,j)+PP(t,j)+II(t,j)) / NN(t,j);//added in eps_h(t,j) to scale foi due to students going on holiday, not sure if this was correct
            }

          }

          //printf("J %i j %i beta %f prob %f NN(t,j) %f\n",J,j,beta(J,j),prob,NN(t,j));
        }

        SE[J] = 1 - exp(-prob - backgroundrate);
    }

    // set AQ
    for(j = 0; j < nages; j++) AQ[j] = eps_testrate * testrate_a_use[j];
    // set PQ
    for(j = 0; j < nages; j++) PQ[j] = eps_testrate * testrate_a_use[j];


    /////////////////////////
    // State Equations
    /////////////////////////
    // discrete-time model
    // SE
    for(j = 0; j < nages; j++)
      {
        i=0.0;
        prob = SE[j];
        //printf("num susceptibles S[%i,%i]=%f prob %f SE %f\n",t,j,SS(t,j),prob,SE[j]);
        i = R::rbinom(SS(t,j), prob);
        //printf("new infections i=%f\n",i);
        SS(t+1,j) = SS(t+1,j) - i;
        //EE(t+1,j) = EE(t+1,j) + i;
        EE(t+1,j) += i;

      }

    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("2 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}

    // EA or EP
    for(j = 0; j < nages; j++) {
      i=0.0; k=0;
      if(EE(t,j)>0)
      {
        prob = 1 - exp(-(EA[j] - EP[j]));
        i = R::rbinom(EE(t,j), prob);
        if(i > 0) {
          prob = EA[j] / (EA[j] + EP[j]);
          k = R::rbinom(i, prob);
        }
      }
      if(isnan(EE(t+1,j))){break;}
      //EE(t+1,j) = EE(t+1,j) - i;
      EE(t+1,j) -= i;
      AA(t+1,j) = AA(t+1,j) + k;
      PP(t+1,j) = PP(t+1,j) + i - k;
      if(EE(t+1,j)<0){printf("i= %f E[%i,%i]=%f, E[%i,%i]=%f\n",i,t,j,EE(t,j),t+1,j,EE(t+1,j)); break;}
      //printf("2. E[%i,%i]=%f\n",t+1,j,EE(t+1,j));
    }


    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("3 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}


    // AR or A-QA
    for(j = 0; j < nages; j++) {
      i=0; k=0.0;
      if(AA(t,j)>0)
      {
        prob = 1 - exp(-AR[j] - AQ[j]);
        i = R::rbinom(AA(t,j), prob);
        if(i>0)
        {
          prob = AR[j] / (AR[j]  + AQ[j]);
          k = R::rbinom(i, prob);
          //printf("k=%f\n",k);
        }
        AA(t+1,j) = AA(t+1,j) - i;
        RR(t+1,j) = RR(t+1,j) + k;
        //QQ(t+1,j) = QQ(t+1,j) + i - k;
        QA(t+1,j) = QA(t+1,j) + i - k;
        newcases = newcases +i -k;
      }

    }


    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("4 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}



    // PI or PQ
    for(j = 0; j < nages; j++) {
      i=0; k=0;
      prob = 1 - exp(-PInf[j] - PQ[j]);
      i = R::rbinom(PP(t,j), prob);
      if(i>0)
      {
        prob = PInf[j] / (PInf[j]  + PQ[j]);
        k = R::rbinom(i, prob);
      }
      PP(t+1,j) = PP(t+1,j) - i;
      II(t+1,j) = II(t+1,j) + k;
      QQ(t+1,j) = QQ(t+1,j) + i - k;
      newcases = newcases + i - k;
    }


    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("5 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}


    // IH or IR or IQ
    for(j = 0; j < nages; j++) {
      i=0; k=0; m=0;
      prob = 1 - exp(-(IH[j] + IR[j] + IQ[j]));
      i = R::rbinom(II(t,j), prob);
      if(i > 0) {
        prob = IH[j] / (IH[j] + IR[j] + IQ[j]);
        k = R::rbinom(i, prob);
        if((i-k) > 0)
        {
          prob = IR[j] / (IR[j] + IQ[j]);
          m = R::rbinom(i-k, prob);
          //printf("m=%f\n",m);
        }
      }
      II(t+1,j) = II(t+1,j) - i;
      HH(t+1,j) = HH(t+1,j) + k;
      RR(t+1,j) = RR(t+1,j) + m;
      QQ(t+1,j) = QQ(t+1,j) + i - k - m;
      newcases = newcases + i -k -m;
    }


    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("6 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}

    // HR or HD
    for(j = 0; j < nages; j++) {
      i=0; k=0;
      prob = 1 - exp(-(HR[j] + HD[j]));
      i = R::rbinom(HH(t,j), prob);
      if(i > 0) {
        prob = HR[j] / (HR[j] + HD[j]);
        k = R::rbinom(i, prob);
      }
      HH(t+1,j) = HH(t+1,j) - i;
      RR(t+1,j) = RR(t+1,j) + i;
      DD(t+1,j) = DD(t+1,j) +i -k;
    }

    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("7 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}

    // QR
    for(j = 0; j < nages; j++) {
      if(QQ(t,j)>0)
      {
        prob = 1 - exp(-QR[j]);
        i = R::rbinom(QQ(t,j), prob);
        QQ(t+1,j) = QQ(t+1,j) - i;
        RR(t+1,j) = RR(t+1,j) + i;
      }
    }

    for(j = 0; j < nages; j++) {
      if(QA(t,j)>0)
      {
        prob = 1 - exp(-QR[j]);
        i = R::rbinom(QA(t,j), prob);
        QA(t+1,j) = QA(t+1,j) - i;
        RR(t+1,j) = RR(t+1,j) + i;
      }
    }


    Ntot=0.0; Ntot2=0.0;
    for(j = 0; j < nages; j++)
    {
      Ntot2 = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);
      Ntot = SS(t,j) + EE(t,j) + AA(t,j) + PP(t,j) + II(t,j) + RR(t,j) + HH(t,j) + QQ(t,j) + QA(t,j);
    }
    if(Ntot!=Ntot2){printf("8 t=%i N1 %f N2 %f\n",t, Ntot,Ntot2); break;}

    for(j = 0; j < nages; j++)
    {
      NN(t+1,j) = SS(t+1,j) + EE(t+1,j) + AA(t+1,j) + PP(t+1,j) + II(t+1,j) + RR(t+1,j) + HH(t+1,j) + QQ(t+1,j) + QA(t+1,j);

      if(isnan(NN(t+1,j)) || EE(t+1,j)<0)
      {
      printf("S=%f E=%f A=%f P=%f I=%f R=%f Q=%f H=%f N=%f\n",SS(t+1,j),EE(t+1,j),AA(t+1,j),PP(t+1,j),II(t+1,j),RR(t+1,j),QQ(t+1,j),HH(t+1,j),NN(t+1,j));
      break;
      }
    }


    // for(j = 0; j < nages; j++){
    //   if(II(t+1,j)>(II(t,j)+5)){
    //     testrate_a_use[j] = testrate_a2;
    //     printf("t=%i j=%i II[t+1,j]=%f II[t,j]=%f, increase testing\n",t+1,j,II(t+1,j),II(t,j));
    //     }
    //   else{testrate_a_use[j] = testrate_a3;}
    // }

    daycount++;
    if(daycount==7)
    {
      if(newcases>0)
      {
        if((newcases/(newcases_then+1)) >= 1)
        {
            for(j = 0; j < nages; j++){
              eps_testrate *testrate_a_use[j] = eps_testrate *testrate_a2;
            }
            //printf("t=%i newcases %i newcasesold %i, increase testing\n",t+1,newcases,newcases_then);
        }
        else{
          for(j = 0; j < nages; j++){
            eps_testrate * testrate_a_use[j] = eps_testrate * testrate_a3;
          }
          //printf("t=%i newcases %i newcasesold %i, decrease testing\n",t+1,newcases,newcases_then);
        }
      }
      newcases_then = newcases;
      newcases=0;
      daycount=1;
    }

  };//end of the t loop

  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = SS,
    Named("E") = EE,
    Named("A") = AA,
    Named("I") = II,
    Named("R") = RR,
    Named("H") = HH,
    Named("P") = PP,
    Named("D") = DD,
    Named("Q") = QQ,
    Named("QA") = QA,
    Named("N") = NN
  );
  return sim;
};
