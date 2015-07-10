DATA_SECTION

//  !!CLASS ofstream post("eval.csv")  
!!CLASS ofstream evalout("atfgoa2013.mcmc.out");
//!!CLASS ofstream post("atfgoa2013.csv");
  init_int styr 	//(1) start year of model
//uses separate M for males and females read in from the data file
//changed wt at age, and transition matrix to include all age data 84-96.
  init_int endyr 	//(2) end year
  init_int styr_fut	//(3) start year of projections (endyr+1)  
  init_int endyr_fut //(4) end year of projections
//read only female M when doing profile like on male M
//  init_number Mf
  init_vector M(1,2) //(5)
  init_int phase_F40 //(6)
  init_number median_rec  //(7)
  init_int median_rec_yrs //(8)
  init_int nages_read //(9)
  init_int nages     //(10)
 //selectivity is set to the selectivity at nselages-1 after age nselages 
  init_int nselages  //(11)
  init_int nselages_srv1 //(12)
  init_int monot_sel     //(13)
  init_int monot_sel_srv1 //(14)
  init_int phase_logistic_sel   //(15)
  init_int phase_selcoffs   //(16)
  init_int phase_logistic_sel_srv1  //(17)
  init_int phase_selcoffs_srv1   //(18)
  init_vector wt_like(1,8)    //(19)
//  !! phase_logistic_sel = -2;
//  !! if (phase_logistic_sel<0) phase_selcoffs=2; else phase_selcoffs=-2;
 //sample size for length comps for weighting likelihoods  
  init_int nlen_r    //(20)
//reduce nlen by one to cut off 1st length bin 
  init_int nobs_fish //(21)
  //!!cout<<nobs_fish<<endl;
  init_ivector yrs_fish(1,nobs_fish) //(22)
  init_matrix nsamples_fish(1,2,1,nobs_fish) //(23)
  init_int nobs_srv1  //(24)
  init_ivector yrs_srv1(1,nobs_srv1)  //(25)
  init_int nobs_srv1_length  //(26)
  init_ivector yrs_srv1_length(1,nobs_srv1_length)    //(27)
  init_matrix nsamples_srv1_length(1,2,1,nobs_srv1_length)  //(28)
  init_int nobs_srv1_age  //(29)
  init_ivector yrs_srv1_age(1,nobs_srv1_age)  //(30)
  init_vector like_wght(1,5)    //(31)
  init_matrix nsamples_srv1_age(1,2,1,nobs_srv1_age) //(32)
  init_3darray obs_p_srv1_len_r(1,2,1,nobs_srv1_length,1,nlen_r) //(33) 2x2x21
  //!!cout<<"obs_p_srv1_len_r"<<endl;
  //!!cout<<obs_p_srv1_len_r<<endl;
  init_3darray obs_p_srv1_age_read(1,2,1,nobs_srv1_age,1,nages_read)  //(34)  2x11x21  
  init_3darray obs_p_fish_r(1,2,1,nobs_fish,1,nlen_r)  //(35) 2x31x21
 // !!cout<<"obs_p_fish"<<endl;
 // !!cout<<obs_p_fish<<endl;
  init_vector catch_bio(styr,endyr)  //(36)
  !!cout<<catch_bio<<endl;
  init_vector obs_srv1(1,nobs_srv1) //(37)
  init_vector obs_srv1_sd(1,nobs_srv1) //(38)
 //need wt vector by length for split sex?
 //init_vector wt(1,nlen) 
  init_matrix wt(1,2,1,nages)   //(39)
  init_vector maturity(1,nages)   //(40)
 //this is how to change the name of the data file
 //can put at the beginning of the file to read everything or
 //part way down
 // LOCAL_CALCS
 //  ad_comm::change_datafile_name("atka.ctl")
 // END_CALCS
 //Local_calcs is indented 1 space
 // LOCAL_CALCS
 //  cout<<maturity<<endl;
 // END_CALCS
//length age transition matrix
  init_3darray lenage(1,2,1,nages,1,nlen_r-1)   //(41)
 //LOCAL_CALCS
 //  cout<<nages<<endl;
 //END_CALCS
   int styr_rec; 
   vector cv_srv1(1,nobs_srv1);
 //year
  int i
//age
  int j
//sex
  int k
//
  int ii
  int m
  int nlen  
 

 LOCAL_CALCS
  nlen=nlen_r-1;
  styr_rec=styr-nages+1;
  if(nselages>nages) nselages=nages;
  if(nselages_srv1>nages) nselages_srv1=nages;
//calculate cv for survey
  cv_srv1=elem_div(obs_srv1_sd,obs_srv1);
//change weights to tons
  wt=wt*.001;
//  nages=nages-10
// cout<<nobs_srv1<<endl;
//  cout<<wt<<endl;
 END_CALCS

  3darray obs_p_srv1_age(1,2,1,nobs_srv1_age,1,nages)
  3darray obs_p_srv1_age_r(1,2,1,nobs_srv1_age,1,nages)
  3darray obs_p_srv1_length(1,2,1,nobs_srv1_length,1,nlen) 
  3darray obs_p_fish(1,2,1,nobs_fish,1,nlen)
  vector obs_sexr(1,nobs_fish)
  vector obs_sexr_srv1(1,nobs_srv1_age)
  vector obs_sexr_srv1_l(1,nobs_srv1_length)

INITIALIZATION_SECTION
//can have different mortality for males and females
//  Mm     .35
  mean_log_rec 18
  log_avg_fmort -5.
  q1 1.
  fmort_dev 0.00001
//  fish_slope_f .4
//  fish_sel50_f  5
//  fish_slope_m  .1
//  fish_sel50_m  8
//  srv1_slope_f  .8
//  srv1_sel50_f  4.
//  srv1_slope_m   .4
//  srv1_sel50_m   8
PARAMETER_SECTION
 //parameters to be estimated are all ones that begin with init_ and have a positive
 //phase, negative phase means are fixed.
 //phase of 8 is greater than last phase so does q1 in last phase  
  // init_bounded_number q1(.5,2,8)
 //fix q1 to be 1 otherwise it went to lower bound of .5
  init_bounded_number q1(0.01,20.0,-8)
 //phase of -1 means M is fixed   
//  init_bounded_number Mm(.1,1.0,8)
  init_number mean_log_rec(1)
  init_bounded_dev_vector rec_dev(styr_rec,endyr,-15,15,2)
  
//  init_vector rec_dev_future(styr_fut,endyr_fut,phase_F40);
  init_number log_avg_fmort(1)
  init_bounded_dev_vector fmort_dev(styr,endyr,-5.0,3.5,1)

  init_matrix log_selcoffs_fish(1,2,1,nselages,phase_selcoffs)
  init_matrix log_selcoffs_srv1(1,2,1,nselages_srv1,phase_selcoffs_srv1)

  init_bounded_number fish_slope_f(.1,5.,phase_logistic_sel)
  init_bounded_number fish_sel50_f(1.,8.,phase_logistic_sel)
  init_bounded_number fish_slope_m(.05,2.0,phase_logistic_sel)
  init_bounded_number fish_sel50_m(1.,25.,phase_logistic_sel)
  init_bounded_number srv1_slope_f(.1,5.,phase_logistic_sel_srv1)
  init_bounded_number srv1_sel50_f(1.,10.,phase_logistic_sel_srv1)
  init_bounded_number srv1_slope_m(.01,5.,phase_logistic_sel_srv1)
  init_bounded_number srv1_sel50_m(1.,10.,phase_logistic_sel_srv1)

  init_bounded_number sexr_param_fish(1.0,1.0,-5)  //this was hitting bound of 1.0 so fixed it - should free up to check
  init_bounded_number sexr_param_srv(.25,1.0,phase_selcoffs_srv1)
// Parameters for computing SPR rates 
  init_bounded_number F40(0.05,.5,phase_F40)
  init_bounded_number F35(0.05,.5,phase_F40)
 // init_bounded_number F30(0.01,1.,phase_F40)

  matrix log_sel_fish(1,2,1,nages)
  matrix log_sel_srv1(1,2,1,nages)
  matrix sel(1,2,1,nages)
  matrix sel_srv1(1,2,1,nages)
  vector avgsel_fish(1,2)
  vector avgsel_srv1(1,2)
  matrix popn(1,2,styr,endyr)
  matrix totn_srv1(1,2,styr,endyr)
  vector explbiom(styr,endyr)
  vector pred_bio(styr,endyr)
  vector fspbio(styr,endyr)
  vector pred_srv1(styr,endyr)
  3darray pred_p_fish(1,2,styr,endyr,1,nlen)
//  3darray pred_p_fish_1(1,2,styr,endyr,1,nlen)
  3darray pred_p_srv1_age(1,2,styr,endyr,1,nages)
//  3darray pred_p_srv1_age_1(1,2,styr,endyr,1,nages)
  3darray pred_p_srv1_len(1,2,styr,endyr,1,nlen)
//  3darray pred_p_srv1_len_1(1,2,styr,endyr,1,nlen)
  vector pred_catch(styr,endyr)
  3darray natage(1,2,styr,endyr,1,nages)
  3darray catage(1,2,styr,endyr,1,nages)
  3darray natlength(1,2,styr,endyr,1,nlen)
  //matrix u(styr,endyr,1,nages)
  3darray Z(1,2,styr,endyr,1,nages)
  3darray F(1,2,styr,endyr,1,nages)
  3darray S(1,2,styr,endyr,1,nages)
  vector fmort(styr,endyr)
  number rbar
  vector surv(1,2)
  vector offset(1,3)
  number rec_like
  number rec_like2
  number catch_like
  vector age_like(1,3)
  vector sel_like(1,4)
  number fpen    
  number surv_like
  sdreport_vector recruits(styr,endyr)
  sdreport_vector biomassrep(styr,endyr)
  sdreport_vector fspbiorep(styr,endyr)
  sdreport_number endbiom
  sdreport_number depletion
  objective_function_value f
  number tmp
  vector pred_sexr(styr,endyr)
  vector preds_sexr(styr,endyr)
 // Stuff for SPR and yield projections
  number sigmar
  number ftmp
  number SB0
  number SBF40
  number SBF35
 // number SBF30
  number sprpen
  matrix Nspr(1,3,1,nages)
  3darray nage_future(1,2,styr_fut,endyr_fut,1,nages)
//  matrix fspbiom_fut(1,5,styr_fut,endyr_fut)
  3darray F_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray Z_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray S_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray catage_future(1,2,styr_fut,endyr_fut,1,nages)
  number avg_rec_dev_future
  vector avg_F_future(1,3)
//  sdreport_matrix catch_future(1,5,styr_fut,endyr_fut)// Note, don't put projection for F=0 (it messes up the hessian matrix)
  matrix catch_future(1,5,styr_fut,endyr_fut) // use this not the sdreport when projecting with F=0 otherwise hessian is screwed up
  sdreport_matrix fspbiom_fut(1,5,styr_fut,endyr_fut)
  sdreport_matrix future_biomass(1,5,styr_fut,endyr_fut) 
  vector explbiom_fut(styr_fut,endyr_fut)
  number maxsel_fish
  vector maxsel_srv1(1,2)
  number B0
  number B40
  number B35
  number AMeanRec
  number like_natm
  number like_q

//do likelihood profile on male M - need to estimate M in model to do this
//run the model with the likelihood profile switch 
//  likeprof_number lp_Mm;
//  vector M(1,2)

//    likeprof_number lp_q1;

PRELIMINARY_CALCS_SECTION
 //cout<<"to prelim calcs"<<endl;
//  M(1)=Mf;
//chop lower ages off and accumulate older ages 
 //sex loop
    for(k=1; k<=2; k++)  
    { 
       for (i=1; i <= nobs_srv1_age; i++)
       {
            for(j=1; j < nages; j++)
            {
                //ages go from 3 to 15
                obs_p_srv1_age_r(k,i,j)=obs_p_srv1_age_read(k,i,j+2);
            }
            for(m=nages; m<=nages_read; m++)
            {
                obs_p_srv1_age_r(k,i,nages)+=obs_p_srv1_age_read(k,i,m);
            }
        }
    }
// cout << " at beg of prel calcs  " <<endl; 
//compute sex ratio in  catch
//  cout<< " sum operation "<< sum(obs_p_fish_r(1,1)(2,nlen+1)) <<endl;
      for(i=1; i<=nobs_fish;i++)
        {
         obs_sexr(i)=(sum(obs_p_fish_r(1,i)(2,nlen+1)))/(sum(obs_p_fish_r(1,i)(2,nlen+1))+sum(obs_p_fish_r(2,i)(2,nlen+1)));
        }
//age obs sex ratio in survey
      for(i=1; i<=nobs_srv1_age;i++)
      {
       obs_sexr_srv1(i)=(sum(obs_p_srv1_age_r(1,i)))/(sum(obs_p_srv1_age_r(1,i))+sum(obs_p_srv1_age_r(2,i)));
      }
//length obs sex ratio in survey
      for(i=1; i<=nobs_srv1_length;i++)
      {
       obs_sexr_srv1_l(i)=(sum(obs_p_srv1_len_r(1,i)(2,nlen+1)))/(sum(obs_p_srv1_len_r(1,i)(2,nlen+1))+sum(obs_p_srv1_len_r(2,i)(2,nlen+1)));
      }
 // cout<< " thru sex ratio "<<endl;
 //Compute offset for multinomial
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_
  //k is sex loop
 offset=0.; 
   for(k=1; k<=2; k++)
    {
       for (i=1; i <= nobs_fish; i++)
       {
         //make observations proportions by year      
         //fishery offset
           for (j=1; j<=nlen; j++)
           {
            obs_p_fish(k,i,j)=(obs_p_fish_r(k,i,j+1))/(sum(obs_p_fish_r(1,i)(2,nlen_r))+sum(obs_p_fish_r(2,i)(2,nlen_r)));
             if (obs_p_fish(k,i,j)>0.0)
              {
               offset(1)-=nsamples_fish(k,i)*obs_p_fish(k,i,j)*log(obs_p_fish(k,i,j));
              }
           }
       }
    }
 //survey length offset 
 for(k=1; k<=2;k++)
 {  
     for (i=1; i <= nobs_srv1_length; i++)
     {
  
     //   cout<< " to obs_p_srv1_length "<<endl;
          for (j=1; j<=nlen; j++)
          {
            obs_p_srv1_length(k,i,j)=obs_p_srv1_len_r(k,i,j+1)/sum(obs_p_srv1_len_r(1,i)(2,nlen_r)+obs_p_srv1_len_r(2,i)(2,nlen_r));
             if (obs_p_srv1_length(k,i,j)>0.0)
              {
                 offset(2)-=nsamples_srv1_length(k,i)*obs_p_srv1_length(k,i,j)*log(obs_p_srv1_length(k,i,j));
              }
          }
     }
 }
 //survey age offset 
    for(k=1; k<=2;k++)
    {  
       for (i=1; i <= nobs_srv1_age; i++)
       {
         obs_p_srv1_age(k,i)=obs_p_srv1_age_r(k,i)/(sum(obs_p_srv1_age_r(1,i))+sum(obs_p_srv1_age_r(2,i)));
         //cout<<obs_p_srv1(i)<<endl;
         for (j=1; j<=nages; j++)
         {
            if (obs_p_srv1_age(k,i,j)>0.0)
             {
               offset(3)-=nsamples_srv1_age(k,i)*obs_p_srv1_age(k,i,j)*log(obs_p_srv1_age(k,i,j));
             }
         }
       }
    }
  //cout<<endl<<"to end of offset"<<endl;


PROCEDURE_SECTION
//   q1=qrun(irun);
//use this when doing like profile on male M
//    M(2)=Mm;
//   lp_Mm=Mm;

//this is if doing likelihood profile
//   lp_q1=q1;

   get_selectivity();
 //cout<<"sel"<<endl;
   get_mortality();
 //cout<<"mort"<<endl;
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));
   get_numbers_at_age();
 //cout<<"numbers at age"<<endl;
   get_catch_at_age();
 //cout<<"catch at age"<<endl;
  if (active(F40))
    compute_spr_rates();
  if (last_phase())
  {
    Future_projections();
  }
  if (sd_phase() || mceval_phase())
  { 
   if (mceval_phase())
   {
      evalout << f << " " ;
      // loop over years and print in one long row.
      for (i=styr;i<=endyr;i++)
      evalout<<  fspbio(i) << " ";
      for (i=styr;i<=endyr;i++)    
      evalout<< natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" ";
      for (i=styr;i<=endyr;i++)
      evalout << 2*natage(1,i,1) <<" ";
      // hit carriage return on file
      evalout <<  endl;
   }  
  }

   evaluate_the_objective_function();

//FUNCTION WriteMCMC
// post<<
//// Mm<<","<<
//// q1<<","<<
//// srv1_slope_f <<","<<
//// srv1_sel50_f <<","<<
//// srv1_slope_m <<","<<
//// srv1_sel50_m <<","<<
////  fish_slope_f <<","<<
////  fish_sel50_f <<","<<
////  fish_slope_m <<","<<
////  fish_sel50_m <<","<<
//  pred_bio <<","<<
//  fspbio <<","<<
//  F40 <<","<<
//  fspbiom_fut <<","<<
//  catch_future <<","<<
//  B0 <<","<<
//  B40 <<","<<
//  B35 <<","<<
//
//  endl;

FUNCTION get_growthmatrix

//  for(ilen=1;ilen<=nlenm;ilen++)
//  {
//    for(il2=1;il2<=nlenm;il2++)
//    {
//
//cumd_norm does the cumulative normal distribution for the value length_bins(il2) with
//a mean of mean_length(ilen) and sd sd_mean_length(ilen)
// add 2.5 to the length bin value to get the upper bound of the length bin since 
//length_bins is the mean of the interval
//      if(il2<2)
//      {
//        len_len(1,ilen,il2)=cumd_norm((length_bins(il2)+2.5-mean_length(1,ilen))/sd_mean_length(1,ilen));
 //       len_len(2,ilen,il2)=cumd_norm((length_bins(il2)+2.5-mean_length(2,ilen))/sd_mean_length(2,ilen));
 //     }
 //     if(il2>1)
 //     {
 //       rec_len(il2)=cumd_norm((length_bins(il2)-length_rec)/var_rec)-cumd_norm((length_bins(il2-1)-length_rec)/var_rec);
 //       len_len(1,ilen,il2)=cumd_norm((length_bins(il2)+2.5-mean_length(1,ilen))/sd_mean_length(1,ilen))-cumd_norm((length_bins(il2-1)+2.5-mean_length(1,ilen))/sd_mean_length(1,ilen));          
 //       len_len(2,ilen,il2)=cumd_norm((length_bins(il2)+2.5-mean_length(2,ilen))/sd_mean_length(2,ilen))-cumd_norm((length_bins(il2-1)+2.5-mean_length(2,ilen))/sd_mean_length(2,ilen));         
     //       cout<<(length_bins(il2)-mean_length(1,ilen))/sd_mean_length(1,ilen)<<endl;
     //        cout<<cumd_norm((length_bins(il2)-mean_length(1,ilen))/sd_mean_length(1,ilen))<<endl;
     //        cout<<len_len(1,ilen,il2-1)<<endl;
  //    }
 //     sum_len(1,ilen)+=len_len(1,ilen,il2);
 //     sum_len(2,ilen)+=len_len(2,ilen,il2);
 //   }
 //   len_len(1,ilen,nlenm)=len_len(1,ilen,nlenm)+1-sum_len(1,ilen);
 //   len_len(2,ilen,nlenm)=len_len(2,ilen,nlenm)+1-sum_len(2,ilen);
 // }

FUNCTION get_selectivity
//cout<<"to begin of sel"<<endl;
//fishery selectivities
  if(active(log_selcoffs_fish))
  {
//turn off logistic curve
   phase_logistic_sel=-2;
    for(k=1;k<=2;k++)
    {
      for (j=1;j<=nselages;j++)
      {
        log_sel_fish(k,j)=log_selcoffs_fish(k,j);
      }
    }
   //sets selectivity of ages older than nselages to selectivity at nselages  
   //cout<<"to nselages loop"<<endl;
   for(k=1;k<=2;k++)
   {
     for (j=nselages+1;j<=nages;j++)
     {
       log_sel_fish(k,j)=log_sel_fish(k,j-1);
     }
   }

 for(k=1;k<=2;k++)
   {
     avgsel_fish(k)=log(mean(mfexp(log_selcoffs_fish(k))));
    }
 //vector=vector-scalar same as  vector-=scalar  
 //scaling selectivities by subracting the mean so exp(mean(s))=1.   
 //selectivities can be greater than 1 but mean is 1.
 //cout<<"calc log_sel_fish"<<endl;
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
      if(k==2)
       {
         sel(k)=sel(k)*sexr_param_fish;
       }
      //cout<<"sel survey"<<sel_srv1<<endl;
    }
 }
  else
    {
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
            if(j<=nselages)
             {
               sel(1,j)=1./(1.+mfexp(-1.*fish_slope_f*(double(j)-fish_sel50_f)));
               sel(2,j)=1./(1.+mfexp(-1.*fish_slope_m*(double(j)-fish_sel50_m)));
             }
            else
            {
             sel(1,j)=sel(1,j-1);
             sel(2,j)=sel(2,j-1);
            }
 
          } 
     }
//survey selectivities
  if(active(log_selcoffs_srv1))
  {
//turn off logistic curve
   phase_logistic_sel_srv1=-2;
   for(k=1;k<=2;k++)
    { 
      for (j=1;j<=nselages_srv1;j++)
      {
        log_sel_srv1(k,j)=log_selcoffs_srv1(k,j);
      }
    }
 //sets selectivity of ages older than nselages_srv1 to selectivity at nselages_srv1
   for(k=1;k<=2;k++)
    {
      for (j=nselages_srv1+1;j<=nages;j++)
       {
         log_sel_srv1(k,j)=log_sel_srv1(k,j-1);
       }
    }

 for(k=1;k<=2;k++)
   {
     avgsel_srv1(k)=log(mean(mfexp(log_selcoffs_srv1(k))));
   }
 //vector=vector-scalar same as  vector-=scalar  
 //scaling selectivities by subracting the mean so exp(mean(s))=1.   
 //selectivities can be greater than 1 but mean is 1.
  for(k=1;k<=2;k++)
    {
      log_sel_srv1(k)-=log(mean(mfexp(log_sel_srv1(k))));
      sel_srv1(k)=mfexp(log_sel_srv1(k));
      if(k==2)
       {
         sel_srv1(k)=sel_srv1(k)*sexr_param_srv;
       }
      //cout<<"sel survey"<<sel_srv1<<endl;
    }
 }
  else
    {
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
          if(j<=nselages_srv1)
            {
              sel_srv1(1,j)=1./(1.+mfexp(-1.*srv1_slope_f*(double(j)-srv1_sel50_f)));
              sel_srv1(2,j)=1./(1.+mfexp(-1.*srv1_slope_m*(double(j)-srv1_sel50_m)));
            }
           else
            {
                sel_srv1(1,j)=sel_srv1(1,j-1);
                sel_srv1(2,j)=sel_srv1(2,j-1);
            }
          } 
     }

FUNCTION get_mortality
  maxsel_fish=max(sel(1));
  if(maxsel_fish<max(sel(2)))
      maxsel_fish=max(sel(2));
//  if(active(log_selcoffs_srv1))
//   {
   maxsel_srv1(1)=max(sel_srv1(1));
   maxsel_srv1(2)=maxsel_srv1(1);

    fmort=mfexp(log_avg_fmort+fmort_dev);
 for(k=1;k<=2;k++)
 {
   for (i=styr;i<=endyr;i++)
   {
        F(k,i)=(sel(k)/maxsel_fish)*fmort(i);
        Z(k,i)=F(k,i) + M(k);
   }
 }
    S=mfexp(-1.0*Z);
 // cout<<"to end of get_mortality"<<endl;

FUNCTION get_numbers_at_age
//  maxsel_fish=max(sel(1));
//  if(maxsel_fish<max(sel(2)))
//      maxsel_fish=max(sel(2));
//  if(active(log_selcoffs_srv1))
//   {
//   maxsel_srv1(1)=max(sel_srv1(1));
//   maxsel_srv1(2)=maxsel_srv1(1);
  if(maxsel_srv1(1)<max(sel_srv1(2)))
     {
      maxsel_srv1(1)=max(sel_srv1(2)); 
      maxsel_srv1(2)=maxsel_srv1(1);
     }
//   }
//  else
//   {
 //   maxsel_srv1(1)=max(sel_srv1(1));
//    maxsel_srv1(2)=max(sel_srv1(2)); 
//   }

 // cout<<"begin numbers at age"<<endl;
  int itmp;

 //calc initial population  
  
  for (j=1;j<nages;j++)
    {
      itmp=styr+1-j;
      natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
      natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));

      //cout<<"natage"<<natage(1,styr,j)<<endl;
    }
    itmp=styr+1-nages;
 //cout<<"initial 1"<<endl;
 //cout<<"to last age "<<endl; 
//last age    
    natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
    natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));

 //cout<<"to next years"<<endl;
 // Now do for next several years----------------------------------
  for (i=styr+1;i<=endyr;i++)
    {
      natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
      natage(2,i,1)=natage(1,i,1);
 //cout<<natage(1,i,1)<<endl;
 //cout<<i<<endl;
    }
//     else
//     {
 //     natage(1,i,1)=median_rec;
 //     natage(2,i,1)=natage(1,i,1);
 //     }
//    }
  //cout<<"initial 2"<<endl;
 //numbers at age
  for(k=1;k<=2;k++)
    {
      for (i=styr;i< endyr;i++)
       {
         //cout<<"to subvector op"<<endl;
         //subvector - avoids writing a j loop  =++ increments the right side 
         //(1,nages-1) to 1+1 to nages-1+1 then does the assignment x(i)(1,n) 
         //takes the ith row of x the columns 1 to n
         //      natage(k,i+1)(2,nages)=++elem_prod(natage(k,i)(1,nages-1),S(k,i)(1,nages-1));
        for(j=1;j<nages;j++)
          {
            natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
          }
         //accumulates oldest ages
         // cout<<"done with j loop"<<endl;
         natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
        // cout<<"done with natage nages"<<endl;
        //popn is exploitable numbers
         popn(k,i)= natage(k,i)*sel(k);
        // cout<<"popn "<<endl;
        // cout<<popn(k,i)<<endl;
     }
     // cout<<"to popn"<<endl; 
     popn(k,endyr)=natage(k,endyr)*sel(k);
    }
  //cout<<" to end of 1st loop"<<endl;
   for (i=styr;i<=endyr;i++)
   {
      for(k=1;k<=2;k++)
      {
      //population numbers at length 
         natlength(k,i)=natage(k,i)*lenage(k); 
      //numbers to rescale age comps and length comps
         totn_srv1(k,i)=(natage(k,i)*sel_srv1(k));
       }
    }
 //cout<<"2nd loop"<<endl;
   //predicted survey values
  fspbio=0.; 
   for (i=styr;i<=endyr;i++)
   {
      fspbio(i)+=natage(1,i)*elem_prod(wt(1),maturity);
      explbiom(i)=0.;
      pred_bio(i)=0.;
      pred_srv1(i)=0.;
      for(k=1;k<=2;k++)
      {
          pred_srv1(i)+=q1*(natage(k,i)*elem_prod(sel_srv1(k)/maxsel_srv1(k),wt(k)));
//         pred_srv1(i)+=q1*(natage(k,i)*elem_prod(sel_srv1(k),wt(k))); 
       //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
        // are set different in the tpl file the program will take to value from the bin file and use that 
       //   pred_srv1(i)=1.0*(natage(i)*elem_prod(sel_srv1,wt));
       explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
//      explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k));
       pred_bio(i)+=natage(k,i)*wt(k);
 // cout<<" to lenage calc"<<endl;
       pred_p_srv1_len(k,i)=elem_prod(sel_srv1(k),natage(k,i))*lenage(k)/(totn_srv1(1,i)+totn_srv1(2,i));
  //   pred_p_srv1_len(k,i)=pred_p_srv1_len_1(k,i)/sum(pred_p_srv1_len_1(1,i)+pred_p_srv1_len_1(2,i)); 
       pred_p_srv1_age(k,i)=elem_prod(sel_srv1(k),natage(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
//       if(i==61 && k==1)
//       {
//        for(j=1;j<=nages;j++)
//        {
//         cout<<(sel_srv1(k,j)*natage(k,i,j))/(totn_srv1(1,i)+totn_srv1(2,i))<<endl;
//        }
 //      }
 //   pred_p_srv1_age(k,i)=pred_p_srv1_age_1(k,i)/sum(pred_p_srv1_age_1(1,i)+pred_p_srv1_age_1(2,i));
      }
    }
//cout<<"3rd loop"<<endl;
//       cout<<" predicted age proportions sex 1 year 61 "<<endl; 
//       cout<<pred_p_srv1_age(1,61)<<endl;
//       cout<<" selectivities "<<endl;
//       cout<<sel_srv1<<endl;
//       cout<<" sel params "<<endl;
//       cout<<srv1_slope_f<<"  "<<srv1_sel50_f<<endl;
//       cout<<srv1_slope_m<<"  "<<srv1_sel50_m<<endl;
//       cout<<" numbers at age sex 1 year 61 "<<endl;
//       cout<<natage(1,61)<<endl;
//       cout<<"total numbers for survey year 61"<<endl;
//       cout<<totn_srv1(1,61)<<endl;
//       cout<<totn_srv1(2,61)<<endl;
//variables for standard dev report
    biomassrep=pred_bio;
    fspbiorep=fspbio;
    for(i=styr;i<=endyr;i++)
    {
       recruits(i)=mfexp(mean_log_rec+rec_dev(i));
    } 
    depletion=pred_bio(endyr)/pred_bio(styr);
    endbiom=pred_bio(endyr);
  //  cout<<"end get numbers at age"<<endl;
FUNCTION get_catch_at_age
  for (i=styr; i<=endyr; i++)
  {
     pred_catch(i)=0.;
      for(k=1;k<=2;k++)
      {      
        //--Baranov's equation here-----------------------------------
         for (j = 1 ; j<= nages; j++)
         {
            catage(k,i,j) = natage(k,i,j)*F(k,i,j)*(1.-S(k,i,j))/Z(k,i,j);
            pred_catch(i)+=catage(k,i,j)*wt(k,j);
         }
        pred_p_fish(k,i)=elem_prod(sel(k),natage(k,i))*lenage(k)/(popn(1,i)+popn(2,i));
//        pred_p_fish(k,i)=pred_p_fish_1(k,i)/(sum(pred_p_fish_1(1,i))+sum(pred_p_fish_1(2,i)));
      }
   }
 // cout<<"end catch at age"<<endl;

FUNCTION Future_projections
 //cout<<"to future proj"<<endl;
  for(k=1;k<=2;k++)
  {
    nage_future(k,styr_fut)(2,nages)=++elem_prod(natage(k,endyr)(1,nages-1),S(k,endyr)(1,nages-1));
    nage_future(k,styr_fut,nages)+=natage(k,endyr,nages)*S(k,endyr,nages);
   }
    future_biomass=0.;
    catch_future=0.;
    fspbiom_fut=0.;
    for (int l=1;l<=5;l++)
    {
      switch (l)
      {
        case 1:
          ftmp=F40;
          break;
        case 2:
          ftmp=F35;
          break;
        case 3:
          ftmp=0.0;
          break;
        case 4:
          ftmp=mean(mfexp(log_avg_fmort+fmort_dev(endyr-4,endyr)));
          break;
//        cout<<ftmp<<endl;
        case 5:
          ftmp=mean(mfexp(log_avg_fmort+fmort_dev(endyr-4,endyr)));
//        cout<<ftmp<<endl;
          break;
      }

      // Get future F's
     for(k=1;k<=2;k++)
     {
      for (i=endyr+1;i<=endyr_fut;i++)
      {
       if(l>3){ftmp=mean(mfexp(log_avg_fmort+fmort_dev(endyr-4,endyr)));}
       if(i>(endyr+1) && l==5) {ftmp=F35;}
       if(i>(endyr+1) && l==4) {ftmp=F40;}
//      cout<<ftmp<<endl;
        for (j=1;j<=nages;j++)
        {
          F_future(k,i,j) = (sel(k,j)/maxsel_fish)*ftmp;
          Z_future(k,i,j) = F_future(k,i,j)+M(k);
          S_future(k,i,j) = mfexp(-1.*Z_future(k,i,j));
        }
      }
    // Future Recruitment (and spawners)
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(k,i,1)  = AMeanRec;
       // Now graduate for the next year....
        nage_future(k,i+1)(2,nages) = ++elem_prod(nage_future(k,i)(1,nages-1),S_future(k,i)(1,nages-1));
        nage_future(k,i+1,nages)   += nage_future(k,i,nages)*S_future(k,i,nages);
      }
      nage_future(k,endyr_fut,1)  = AMeanRec;
      // Now get catch at future ages
        fspbiom_fut(l)=0.;
//        catch_future(l)=0.;
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        for (j = 1 ; j<= nages; j++)
        {
          catage_future(k,i,j) = nage_future(k,i,j) * F_future(k,i,j) * ( 1.- S_future(k,i,j) ) / Z_future(k,i,j);
         if(k==1)
          {
          fspbiom_fut(l,i) += nage_future(1,i,j)*wt(1,j)*maturity(j);
          }
        }
        if (l!=3){
            catch_future(l,i)   += catage_future(k,i)*wt(k);
           }
//         if(l==3){
//            catch_future(l,i) = 0.0;
//            }
        future_biomass(l,i) += nage_future(k,i)*wt(k);
 
      }   //end loop over future years
     }   //end loop over sex
     fspbiom_fut(l)=0.;
     for(i=styr_fut;i<=endyr_fut;i++)
     {
      for(j=1;j<=nages;j++)
      {
        fspbiom_fut(l,i)+= nage_future(1,i,j)*wt(1,j)*maturity(j);
      }
     }
    }   //End of loop over F's
  
FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0=0.;
  SBF40=0.;
  SBF35=0.;
  //SBF30=0.;
//average recruitment from 1981 to last year of estimated recruitment to make it consistent
//with projection model (recruits start at 1978 age 1, so age 3 1981)
  AMeanRec=mean(mfexp(mean_log_rec+rec_dev(81,endyr)));
  // Initialize the recruit (1) for each F  (F40 etc)
  for (i=1;i<=3;i++)
    Nspr(i,1)=1.;

  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*mfexp(-1.*M(1));
    Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(M(1)+F40*sel(1,j-1)/maxsel_fish));
    Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(M(1)+F35*sel(1,j-1)/maxsel_fish));
  }
  //cout<<F40<<" "<<F30<<" "<<Nspr<<endl; 
 // cout<<"spr calc"<<endl;
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*mfexp(-1.*M(1))/(1.-mfexp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/(1.-mfexp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/(1.-mfexp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
 // Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-mfexp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));
 //cout<<"plus group"<<endl;
  for (j=1;j<=nages;j++)
  {
   // Kill them off till- use fraction of the year e.g. april=0.25 - not for atf they spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*mfexp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
  //  SBF30  += Nspr(3,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F30*sel(1,j)/maxsel_fish));
   }
 // cout<<"kill thenm off"<<endl;
  sprpen    = 300.*square((SBF40/SB0)-0.4);
//  sprpen   += 300.*square((SBF30/SB0)-0.30);
  sprpen   += 300.*square((SBF35/SB0)-0.35);
//AMeanRec is mean of recruitment(2*recruitment is total recruitment(male+female))
  B0  = AMeanRec*SB0  ;
  B40 = AMeanRec*SBF40  ;
  B35 = AMeanRec*SBF35  ;

 //cout<<"sbr/sno "<<endl;
 //cout<<SBF40/SB0<<" "<<SBF30/SB0<<endl;
FUNCTION evaluate_the_objective_function
  age_like=0.;
  sel_like=0.;
  fpen=.0;
  rec_like=.0;
  rec_like2=0.;
  surv_like=.0;
  catch_like=.0;
  f=.0;
  like_natm=0.;
  like_q=0.;
 if (active(rec_dev))
   {
    age_like=0.;
    int ii;
    //recruitment likelihood - norm2 is sum of square values   
    //cout<<"to rec_like"<<endl;
    //cout<<" rec_devs = "<<rec_dev<<endl; 
    rec_like=norm2(rec_dev(styr_rec,endyr));
//    rec_like2=norm2(rec_dev(endyr-median_rec_yrs+1,endyr))/(2*norm2(rec_dev(81,endyr-median_rec_yrs))/(endyr-median_rec_yrs-1981+1-1)+0.001);
// 0.038 is the variance of the rec_dev from 1981 to 1999 - this restricts recruits in last years
//to mean log recruits if don't have survey data in current year.
//   rec_like2 =0.00000001*0.4*(5./(2*.038))*norm2(rec_dev(endyr-median_rec_yrs+1,endyr)-mean(rec_dev(81,endyr-median_rec_yrs)));

//  cout<<rec_like<<" , "<<rec_like2<<endl;
    f+=rec_like;
    f+=rec_like2;
    //cout<<"to loop"<<endl;
//   if(active(fish_slope_f))
//    {
     for(k=1;k<=2;k++)
     {
        for (i=1; i <= nobs_fish; i++)
        {
          ii=yrs_fish(i);
          //   cout<<pred_p_fish(ii)<<endl<<endl;
          //fishery length likelihood 
           for (j=1; j<=nlen; j++)
           {
             age_like(1)-=nsamples_fish(k,i)*(1e-5+obs_p_fish(k,i,j))*log(pred_p_fish(k,ii,j)+1e-5);
             //cout << age_like << " " << i << " "<<j<< endl;
            }
         }
      }
     //add the offset to the likelihood   
     age_like(1)-=offset(1);
     //survey
       for(k=1;k<=2;k++)
       {
         for (i=1; i <=nobs_srv1_length; i++)
         {
            ii=yrs_srv1_length(i);
            //survey likelihood 
             for (j=1; j<=nlen; j++)
             {
               age_like(2)-=nsamples_srv1_length(k,i)*(1e-3+obs_p_srv1_length(k,i,j))*log(pred_p_srv1_len(k,ii,j)+1e-3);
               //cout << age_like << " " << i << " "<<j<< endl;
             }
         }
       }
     age_like(2)-=offset(2);
  //bracket for active(fish_slope_f)
//   }   
  //survey ages
         for(k=1;k<=2;k++)
         {
           for (i=1; i <=nobs_srv1_age; i++)
           {
              ii=yrs_srv1_age(i);
              //survey likelihood 
              for (j=1; j<=nages; j++)
              {
                age_like(3)-=nsamples_srv1_age(k,i)*(1e-3+obs_p_srv1_age(k,i,j))*log(pred_p_srv1_age(k,ii,j)+1e-3);
 //               cout << age_like << " " << i << " "<<j<< endl;
              }
            }
         } 
      age_like(3)-=offset(3);
      f+=like_wght(1)*age_like(1);     
      f+=like_wght(2)*age_like(2);
      f+=like_wght(3)*age_like(3);
   //end of if(active (rec_dev))
    }
  //cout<<" f before survey = "<<f<<endl;
  // Fit to indices (lognormal) 
  //weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 

   surv_like = norm2(elem_div(log(obs_srv1+.000001)-log(pred_srv1(yrs_srv1)+.000001),sqrt(2)*sqrt(log(elem_prod(cv_srv1,cv_srv1)+1.0))));
 //this subtracts the log(sd) from the likelihood - is a constant so I'm not adding it.
 //   surv_like-= sum(log(sqrt(log(elem_prod(cv_srv1,cv_srv1)+1.0))));

   // surv_like = norm2(log(obs_srv1+.01)-log(pred_srv1(yrs_srv1)+.01));

    catch_like=norm2(log(catch_bio+.000001)-log(pred_catch+.000001));
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences 
  if(active(log_selcoffs_fish))
  {  
    sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1))));
    sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2))));
   for (j=1;j<nages;j++)
   {
    if(monot_sel>0.0001)
    { 
     if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
        sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));
      if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
        sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));
     }

   } 
    f+=1.*(sel_like(1)+sel_like(3));

    f+=1.*square(avgsel_fish(1));
    f+=1.*square(avgsel_fish(2));
  }
  if(active(log_selcoffs_srv1))
  {  
    sel_like(2)=wt_like(2)*norm2(first_difference(first_difference(log_sel_srv1(1))));
    sel_like(4)=wt_like(4)*norm2(first_difference(first_difference(log_sel_srv1(2)))); 
   for (j=1;j<nages;j++)
   {
    if(monot_sel_srv1>0.0001)
     {
      if (log_sel_srv1(1,j)>log_sel_srv1(1,j+1))
        sel_like(2)+=wt_like(7)*square(log_sel_srv1(1,j)-log_sel_srv1(1,j+1));
      if (log_sel_srv1(2,j)>log_sel_srv1(2,j+1))
        sel_like(4)+=wt_like(8)*square(log_sel_srv1(2,j)-log_sel_srv1(2,j+1));
     }
   } 
    f+=1.*(sel_like(2)+sel_like(4));

    f+=1.*square(avgsel_srv1(1));
    f+=1.*square(avgsel_srv1(2));
  }

//put in mult by more than other likehs to improve fit to survey biomass
    f+=like_wght(5)*surv_like;
  //cout<<" f after survey = "<<f<<endl;  
    f+=like_wght(4)*catch_like;
  //cout<<" f after catch = "<<f<<endl;
  // Phases less than 3, penalize High F's
    if (current_phase()<2)
     {
       //F's are low for arrowtooth changed the value to compare from .2 to .001
       //don't know if makes any difference since the penalty is reduced at the end
       fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-.001);
     }
    else
     {
       fpen=.01*norm2(mfexp(fmort_dev+log_avg_fmort)-.001);
     }

    if (active(fmort_dev))
      {
        fpen+=norm2(fmort_dev);
      }
   f+=fpen;
   f+=sprpen;
//bayesian on nat mortalities
//  if(active(Mm))
//   {  
//    like_natm   += .5 * square((Mm    - 0.35)    / (0.35*10.0));
//  like_natm   += .5 * square((M(1)    - 0.2)    / (0.2*.5));
//  like_natm   += .5 * square((M(2)    - 0.35)    / (0.35*.5));
//    f += like_natm;
//   }

//bayesian prior on q
//  if(active(q1))
//    {
//      like_q   += .5 * square((q1    - 1.34)    / (1.34*0.1));
//     f += like_q;
//      }

//   cout<<"q= "<<q1<<endl;
//   cout<<"pred bio end year = "<<pred_bio(endyr)<<endl;
   // cout<<" f after fpen = "<<f<<endl;
REPORT_SECTION
  report << "Estimated numbers of female fish: seq(3,15) " << endl;
  report << natage(1) << endl;
  report << "Estimated numbers of male fish: seq(3,15) " << endl;
  report << natage(2) << endl;
  report << "Estimated numbers of female fish by length: '25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
  report << natlength(1) << endl;
  report << "Estimated numbers of male fish by length: '25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
  report << natlength(2) << endl;
  report << "Estimated numbers of female fish by length MD: '24','30','40','52','61','70'" << endl;
  for(i=styr;i<=endyr;i++){
  report <<sum(natlength(1,i)(1,3))<<" "<<sum(natlength(1,i)(4,8))<<" "<<sum(natlength(1,i)(9,12))<<" "<<sum(natlength(1,i)(13,15))<<" "<<sum(natlength(1,i)(16,18))<<" "<<sum(natlength(1,i)(19,20))<<endl;
  }
  report << "Estimated numbers of male fish by length MD: '25','35','45','55','65','75'" << endl;
  for(i=styr;i<=endyr;i++){
  report <<sum(natlength(2,i)(1,3))<<" "<<sum(natlength(2,i)(4,8))<<" "<<sum(natlength(2,i)(9,12))<<" "<<sum(natlength(2,i)(13,15))<<" "<<sum(natlength(2,i)(16,18))<<" "<<sum(natlength(2,i)(19,20))<<endl;
  }

//  report << "Estimated F mortality: seq(3,15)" << endl;
//  report << F << endl;
  report << "Observed Survey 1:  '1961','1973','1984','1987','1990','1993','1996','1999','2001','2003','2005','2007','2009','2011','2013'" << endl;
  report << obs_srv1 << endl;
  report << "Survey Q: 'Q'" << endl;
  report << q1 << endl;
  report << "natural mortality females, males: 'FemM','MaleM'" << endl;
  report << M << endl;
  report << "Predicted Survey 1: seq(1961,"<<endyr+1900<<") " << endl;
  report << pred_srv1 << endl;
  report << "Predicted Biomass : seq(1961,"<<endyr+1900<<")" << endl;
  report << pred_bio << endl;
  report << "Predicted Exploitable Biomass : seq(1961,"<<endyr+1900<<") " << endl;
  report << explbiom << endl;
  report << "Female Spawning Biomass : seq(1961,"<<endyr+1900<<") " << endl;
  report << fspbio << endl;
  report << "Observed Prop fishery females:'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'"<< endl;
      for (i=1; i<=nobs_fish; i++)
        {
          report << yrs_fish(i)+1900 << " " << obs_p_fish(1,i) << endl;
         }
    report << "Predicted length prop fishery females : 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
   //report<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          report <<  ii+1900  <<  " "  <<  pred_p_fish(1,ii)  << endl;
         }
  report << "Observed length Prop fishery males:'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'"<< endl;
      for (i=1; i<=nobs_fish; i++)
        {
          report << yrs_fish(i)+1900 << " " << obs_p_fish(2,i) << endl;
         }
    report << "Predicted length prop fishery males : 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
   //report<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          report <<  ii +1900 <<  " "  <<  pred_p_fish(2,ii)  << endl;
         }
 report << "Observed Length Prop survey females: 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
              ii=yrs_srv1_length(i);
              report << ii +1900 <<" " <<obs_p_srv1_length(1,i) << endl;
           }
  report << "Predicted length prop survey females: 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             report << ii +1900<< " " << pred_p_srv1_len(1,ii) << endl;
          }
 report << "Observed Length Prop survey males: 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             report << yrs_srv1_length(i)+1900 <<" " <<obs_p_srv1_length(2,i) << endl;
           }
  report << "Predicted length prop survey males: 'year','25','27','29','31','33','35','37','39','41.5','44.5','47.5','50.5','53.5','56.5','59.5','62.5','65.5','68.5','72.5','85'" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             report << ii +1900<< " " << pred_p_srv1_len(2,ii) << endl;
          }
 report << "Observed Age Prop survey females: 'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
         for (i=1; i<=nobs_srv1_age; i++)
          {
             report << yrs_srv1_age(i)+1900 <<" " <<obs_p_srv1_age(1,i) << endl;
           }
  report << "Predicted age prop survey females: 'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
       for (i=1; i<=nobs_srv1_age; i++)
          {
             ii=yrs_srv1_age(i);  
             report << ii +1900<< " " << pred_p_srv1_age(1,ii) << endl;
          }
 report << "Observed age Prop survey males: 'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
         for (i=1; i<=nobs_srv1_age; i++)
          {
             report << yrs_srv1_age(i)+1900 <<" " <<obs_p_srv1_age(2,i) << endl;
           }
  report << "Predicted age prop survey males: 'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
       for (i=1; i<=nobs_srv1_age; i++)
          {
             ii=yrs_srv1_age(i);  
             report << ii +1900 << " " << pred_p_srv1_age(2,ii) << endl;
          }
  //report << pred_p_srv1<<endl;
  report << "observed catch biomass : seq(1961,"<<endyr+1900<<")" << endl;
  report << catch_bio << endl;
  report << "predicted catch biomass : seq(1961,"<<endyr+1900<<")" << endl;
  report << pred_catch << endl;
  report << "estimated annual fishing mortality : seq(1961,"<<endyr+1900<<")" << endl;
  report << mfexp(log_avg_fmort+fmort_dev) << endl;
  report << "initial number of recruitments age 3 : seq(1949,1960)" << endl;
      for(i=styr+1-nages; i<=styr-1; i++)
       {
          report << 2*mfexp(mean_log_rec+rec_dev(i))<<" ";
       } 
  report <<endl<< "estimated number of recruitments age 3(males+females) : seq(1961,"<<endyr+1900<<")" << endl;
      for(i=styr; i<=endyr; i++)
       {
          report << 2*mfexp(mean_log_rec+rec_dev(i))<<" ";
       }
   //   for(i=1;i<=median_rec_yrs;i++)report<<2*median_rec<<" ";
  report <<endl<< "selectivity fishery females: seq(3,15)" << endl;
  report << sel(1)/maxsel_fish << endl;
  report << "selectivity fishery males: seq(3,15)" << endl;
  report << sel(2)/maxsel_fish << endl;
  report << "selectivity survey females: seq(3,15)" << endl;  
  report << sel_srv1(1)/maxsel_srv1(1) << endl;
  report << "selectivity survey males: seq(3,15)" << endl;  
  report << sel_srv1(2)/maxsel_srv1(2) << endl;
  report << "obs_sexr: seq(1,"<<nobs_fish<<")" << endl;
  report << yrs_fish <<endl;
   report << obs_sexr<<endl;
  for (i=styr;i<=endyr;i++)
  {
  pred_sexr(i)=popn(1,i)/(popn(1,i)+popn(2,i));
  }
  report << "pred_sexr: seq(1961,"<<endyr+1900<<")" << endl;
  report << pred_sexr <<endl;
 for (i=styr;i<=endyr;i++)
  {
  preds_sexr(i)=totn_srv1(1,i)/(totn_srv1(1,i)+totn_srv1(2,i));
  }
  report << "obs_sexr_srv1 from lengths: seq(1,"<<nobs_srv1_length<<")" << endl;
  report << yrs_srv1_length <<endl;
   report << obs_sexr_srv1_l<<endl; 
  report << "obs_sexr_srv1 from ages: seq(1,"<<nobs_srv1_age<<")" << endl;
  report << yrs_srv1_age <<endl;
   report << obs_sexr_srv1<<endl;  
  report << "pred_sexr survey: seq(1961,"<<endyr+1900<<")" << endl;
  report << preds_sexr <<endl;
  report <<"likelihood : 'rec_like',  'len_like_fishery', 'len_like_surv', 'age_like_surv', 'fpen',  'catch_like',  'surv_like',' sel_like_fish',' sel_like_survey','total likelihood'"<<endl;
  report <<rec_like<<"  "<< age_like<< " "<<fpen<<" "<<catch_like<<" "<<surv_like<<" "<<sel_like(1)+sel_like(3)<<" "<<sel_like(2)+sel_like(4)<<" "<<f<<endl;
  report <<"projected biomass: seq("<<endyr+1900+1<<","<<endyr+1900+5<<")" << endl;
  report <<future_biomass<<endl;
  report <<"projected female spawning biomass: seq("<<endyr+1900+1<<","<<endyr+1900+5<<")" << endl;
  report <<fspbiom_fut<<endl;
  report <<"Future Yield F40 F35 : seq("<<endyr+1900+1<<","<<endyr+1900+5<<")"<<endl;
  report <<catch_future<<endl;
  report <<"Future Fishing mortalities: 'F40','F35'"<<endl;
  report <<F40<<" "<<F35<<endl;
  report <<"Biomass ref points: 'B0','B40','B35'"<<endl;
  report<<B0<<"  "<<B40<<"  "<<B35<<endl;
  report <<"mean recruitment: 'Ameanrec'"<<endl;
  report <<AMeanRec<<endl;
  report << "future Estimated numbers of female fish: seq(3,15) " << endl;
  report << nage_future(1) << endl;
  report << "future Estimated numbers of male fish: seq(3,15) " << endl;
  report << nage_future(2) << endl;

    if (last_phase())
  {
// r output file  
  R_out << "$Estimated.numbers.of.female.fish" << endl;
  R_out << natage(1) << endl;
  R_out << "$Estimated.numbers.of.male.fish" << endl;
  R_out << natage(2) << endl;
  R_out << "$Estimated.numbers.of.female.fish.by.length" << endl;
  R_out << natlength(1) << endl;
  R_out << "$Estimated.numbers.of.male.fish.by.length" << endl;
  R_out << natlength(2) << endl;
  R_out << "$Estimated.numbers.of.female.fish.by.length.MD" << endl;
  for(i=styr;i<=endyr;i++){
  R_out <<sum(natlength(1,i)(1,3))<<" "<<sum(natlength(1,i)(4,8))<<" "<<sum(natlength(1,i)(9,12))<<" "<<sum(natlength(1,i)(13,15))<<" "<<sum(natlength(1,i)(16,18))<<" "<<sum(natlength(1,i)(19,20))<<endl;
  }
  R_out << "$Estimated.numbers.of.male.fish.by.length.MD" << endl;
  for(i=styr;i<=endyr;i++){
  R_out <<sum(natlength(2,i)(1,3))<<" "<<sum(natlength(2,i)(4,8))<<" "<<sum(natlength(2,i)(9,12))<<" "<<sum(natlength(2,i)(13,15))<<" "<<sum(natlength(2,i)(16,18))<<" "<<sum(natlength(2,i)(19,20))<<endl;
  }

//  report << "Estimated F mortality: seq(3,15)" << endl;
//  report << F << endl;
  R_out << "$Observed.Survey" << endl;
  R_out << obs_srv1 << endl;
  R_out << "$Observed.Survey.sd" << endl;
  R_out << obs_srv1_sd << endl;
  R_out << "$Survey.Q" << endl;
  R_out << q1 << endl;
  R_out << "$natural.mortality.females.males" << endl;
  R_out << M << endl;
  R_out << "$Predicted.Survey" << endl;
  R_out << pred_srv1 << endl;
  R_out << "$Predicted.Biomass" << endl;
  R_out << pred_bio << endl;
  R_out << "$Predicted.Exploitable.Biomass" << endl;
  R_out << explbiom << endl;
  R_out << "$Female.Spawning.Biomass" << endl;
  R_out << fspbio << endl;
  R_out << "$Observed.Prop.fishery.females"<< endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i)+1900 << " " << obs_p_fish(1,i) << endl;
         }
    R_out << "$Predicted.length.prop.fishery.females" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii+1900  <<  " "  <<  pred_p_fish(1,ii)  << endl;
         }
  R_out << "$Observed.length.Prop.fishery.males"<< endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i)+1900 << " " << obs_p_fish(2,i) << endl;
         }
    R_out << "$Predicted.length.prop.fishery.males" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii +1900 <<  " "  <<  pred_p_fish(2,ii)  << endl;
         }
 R_out << "$Observed.Length.Prop.survey.females" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
              ii=yrs_srv1_length(i);
              R_out << ii +1900 <<" " <<obs_p_srv1_length(1,i) << endl;
           }
  R_out << "$Predicted.length.prop.survey.females" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii +1900<< " " << pred_p_srv1_len(1,ii) << endl;
          }
 R_out << "$Observed.Length.Prop.survey.males" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             R_out << yrs_srv1_length(i)+1900 <<" " <<obs_p_srv1_length(2,i) << endl;
           }
  R_out << "$Predicted.length.prop.survey.males" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii +1900<< " " << pred_p_srv1_len(2,ii) << endl;
          }
 R_out << "$Observed.Age.Prop.survey.females" << endl;
         for (i=1; i<=nobs_srv1_age; i++)
          {
             R_out << yrs_srv1_age(i)+1900 <<" " <<obs_p_srv1_age(1,i) << endl;
           }
  R_out << "$Predicted.age.prop.survey.females" << endl;
       for (i=1; i<=nobs_srv1_age; i++)
          {
             ii=yrs_srv1_age(i);  
             R_out << ii +1900<< " " << pred_p_srv1_age(1,ii) << endl;
          }
 R_out << "$Observed.age.Prop.survey.males" << endl;
         for (i=1; i<=nobs_srv1_age; i++)
          {
             R_out << yrs_srv1_age(i)+1900 <<" " <<obs_p_srv1_age(2,i) << endl;
           }
  R_out << "$Predicted.age.prop.survey.males" << endl;
       for (i=1; i<=nobs_srv1_age; i++)
          {
             ii=yrs_srv1_age(i);  
             R_out << ii +1900 << " " << pred_p_srv1_age(2,ii) << endl;
          }
  //R_out << pred_p_srv1<<endl;
  R_out << "$observed.catch.biomass" << endl;
  R_out << catch_bio << endl;
  R_out << "$predicted.catch.biomass" << endl;
  R_out << pred_catch << endl;
  R_out << "$estimated.annual.fishing.mortality" << endl;
  R_out << mfexp(log_avg_fmort+fmort_dev) << endl;
  R_out << "$initial.number.of.recruitments.age.3" << endl;
      for(i=styr+1-nages; i<=styr-1; i++)
       {
          R_out << 2*mfexp(mean_log_rec+rec_dev(i))<<" ";
       } 
  R_out <<endl<< "$estimated.number.of.recruitments.age.3" << endl;
      for(i=styr; i<=endyr; i++)
       {
          R_out << 2*mfexp(mean_log_rec+rec_dev(i))<<" ";
       }
   //   for(i=1;i<=median_rec_yrs;i++)R_out<<2*median_rec<<" ";
  R_out <<endl<< "$selectivity.fishery.females" << endl;
  R_out << sel(1)/maxsel_fish << endl;
  R_out << "$selectivity.fishery.males" << endl;
  R_out << sel(2)/maxsel_fish << endl;
  R_out << "$selectivity.survey.females" << endl;  
  R_out << sel_srv1(1)/maxsel_srv1(1) << endl;
  R_out << "$selectivity.survey.males" << endl;  
  R_out << sel_srv1(2)/maxsel_srv1(2) << endl;
  R_out << "$obs.sexr" << endl;
  R_out << yrs_fish <<endl;
   R_out << obs_sexr<<endl;
  for (i=styr;i<=endyr;i++)
  {
  pred_sexr(i)=popn(1,i)/(popn(1,i)+popn(2,i));
  }
  R_out << "$pred.sexr" << endl;
  R_out << pred_sexr <<endl;
 for (i=styr;i<=endyr;i++)
  {
  preds_sexr(i)=totn_srv1(1,i)/(totn_srv1(1,i)+totn_srv1(2,i));
  }
  R_out << "$obs.sexr.srv1.from.lengths" << endl;
  R_out << yrs_srv1_length <<endl;
   R_out << obs_sexr_srv1_l<<endl; 
  R_out << "$obs.sexr.srv1.from.ages" << endl;
  R_out << yrs_srv1_age <<endl;
   R_out << obs_sexr_srv1<<endl;  
  R_out << "$pred.sexr.survey" << endl;
  R_out << preds_sexr <<endl;
  R_out <<"$likelihood"<<endl;
  R_out <<rec_like<<"  "<< age_like<< " "<<fpen<<" "<<catch_like<<" "<<surv_like<<" "<<sel_like(1)+sel_like(3)<<" "<<sel_like(2)+sel_like(4)<<" "<<f<<endl;
  R_out <<"$projected.biomass" << endl;
  R_out <<future_biomass<<endl;
  R_out <<"$projected.female.spawning.biomass" << endl;
  R_out <<fspbiom_fut<<endl;
  R_out <<"$Future.Yield.F40.F35"<<endl;
  R_out <<catch_future<<endl;
  R_out <<"$Future.Fishing.mortalities"<<endl;
  R_out <<F40<<" "<<F35<<endl;
  R_out <<"$Biomass.ref.points"<<endl;
  R_out<<B0<<"  "<<B40<<"  "<<B35<<endl;
  R_out <<"$mean.recruitment"<<endl;
  R_out <<AMeanRec<<endl;
  R_out << "$future.Estimated.numbers.of.female.fish" << endl;
  R_out << nage_future(1) << endl;
  R_out << "$future.Estimated.numbers.of.male.fish" << endl;
  R_out << nage_future(2) << endl;  

  }
// ===============================================================================


  


  

  report << "SARA form for Angie Grieg" << endl;

  report << "ATF        # stock  " << endl;
  report << "BSAI       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
  report << "2013       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
  report << "3a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
  report << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
  report << "partial    # UPDATE (new benchmark full partial)" << endl;
  report << "2          # LIFE_HIST - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "2          # ASSES_FREQ - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "5          # ASSES_LEV - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "5          # CATCH_DAT - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "3          # ABUND_DAT - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "567000     # Minimum B  Lower 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "665000     # Maximum B  Upper 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "202138     # BMSY  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
  report << "ADMB       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
  report << "           # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
  report << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
  report << "1          # number of fisheries" << endl;
  report << "1          # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
  report << "1          # recruitment age used by model or size" << endl;
  report << "1          # age+ or mmCW+ used for biomass estimate" << endl;
  report << "Single age        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
  report << "Age model         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
  report << "Age of maximum F  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
  report << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
  report << "ALL" << endl; 

  report <<"#FISHERYYEAR - list years used in the model " << endl;
    for (i=styr;  i<=endyr; i++)
       report << i << "	";
       report<<endl;  

  report<<"#AGE - list of ages used in the model"<<endl;
    for (i=1; i<=21;i++)
       report << i << "	";
       report<<endl;    

  report <<"#RECRUITMENT - Number of recruits by year " << endl;
    for (i=styr;  i<=endyr;  i++)
	   report  << 2*natage(1,i,1) << "	";
	   report<<endl;     

  report <<"#SPAWNBIOMASS - Spawning biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << fspbio(i) << "	";
       report<<endl;  

  report <<"#TOTALBIOMASS - Total biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << natage(1,i)*wt(1) + natage(2,i)*wt(2) << "	";
       report<<endl;

  report <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
	for (i=styr;  i<=endyr;  i++)
	   report  << (F(1,i,13)+ F(2,i,13))/2<< "	";
	   report<<endl;

  report <<"#TOTALCATCH - Total catch by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << catch_bio(i) << "	";
       report<<endl;

 report <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
       for (i=1;  i<=13;  i++) 
       report  << maturity(i) <<"	";
       report<< endl; 

 report <<"#SPAWNWT - Average spawning weight (in kg) by age"<< endl; 
       report <<"0.019	0.041	0.113	0.224	0.376	0.566	0.784	1.028	1.292	1.569	1.855	2.142	2.417	2.667	2.881	3.057	3.198	3.308	3.393"<<endl;                              
       report<<endl;

 report <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
 for (i=1;  i<=13;  i++) 
 report  << 0.2 <<"	";
 report<< endl;   
 for (i=1;  i<=13;  i++) 
 report  << 0.35 <<"	";
 report<< endl;

 report << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
   for (i=styr; i<=endyr;i++)
     report <<natage(1,i)<< "	";
     report<<endl;

   for (i=styr; i<=endyr;i++)
     report <<natage(2,i)<< "	";
     report<<endl;

 report <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
    report <<wt(1)*1000  << "	";
    report<<endl; //1 is females        

    report <<wt(2)*1000  << "	";
    report<<endl; //2 is males


  report << "#SELECTIVITY - Estimated fishery selectivity for females (first) males (second) at age " << endl;
    for (j=1; j<=nages;j++)
      report <<" " <<sel(1,j)<< "	";
      report<<endl;

    for (j=1; j<=nages;j++)
      report <<" "  <<sel(2,j)<< "	";
      report<<endl;
 report << "#SURVEYDESC"<<endl;
 report<<"EBS_trawl_survey BS_slope_trawl_survey AI_trawl_survey"<<endl;

 report<<"SURVEYMULT"<<endl;
 report<<"1 1 1"<<endl;
 report << "#GOA_trawl_survey - Gulf of Alaska survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nobs_srv1;i++)
      report << yrs_srv1(i) << "	";
      report<<endl;
    for (i=1; i<=nobs_srv1;i++) 
      report<< obs_srv1(i)<< "	";
      report<< endl;

 report<<"#STOCKNOTES"<<endl;
 report<<"SAFE report indicates that this stock was not subjected to overfishing in 2012 and is neither overfished nor approaching a condition of being overfished in 2013."<<endl;
 
GLOBALS_SECTION
 #include <math.h>
 #include <admodel.h>
  #include <time.h>
  ofstream R_out;
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

RUNTIME_SECTION
  maximum_function_evaluations 4000
  convergence_criteria 1e-7

TOP_OF_MAIN_SECTION
  arrmblsize = 2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(300);
  time(&start);
//  CheckFile.open("Check.Out");
  R_out.open("R_input.txt");

FINAL_SECTION
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;       

