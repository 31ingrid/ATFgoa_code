 #include <math.h>
 #include <admodel.h>
  #include <time.h>
  ofstream R_out;
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <atfgoa2013.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("atfgoa2013.mcmc.out");;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_fut.allocate("styr_fut");
  endyr_fut.allocate("endyr_fut");
  M.allocate(1,2,"M");
  phase_F40.allocate("phase_F40");
  median_rec.allocate("median_rec");
  median_rec_yrs.allocate("median_rec_yrs");
  nages_read.allocate("nages_read");
  nages.allocate("nages");
  nselages.allocate("nselages");
  nselages_srv1.allocate("nselages_srv1");
  monot_sel.allocate("monot_sel");
  monot_sel_srv1.allocate("monot_sel_srv1");
  phase_logistic_sel.allocate("phase_logistic_sel");
  phase_selcoffs.allocate("phase_selcoffs");
  phase_logistic_sel_srv1.allocate("phase_logistic_sel_srv1");
  phase_selcoffs_srv1.allocate("phase_selcoffs_srv1");
  wt_like.allocate(1,8,"wt_like");
  nlen_r.allocate("nlen_r");
  nobs_fish.allocate("nobs_fish");
  yrs_fish.allocate(1,nobs_fish,"yrs_fish");
  nsamples_fish.allocate(1,2,1,nobs_fish,"nsamples_fish");
  nobs_srv1.allocate("nobs_srv1");
  yrs_srv1.allocate(1,nobs_srv1,"yrs_srv1");
  nobs_srv1_length.allocate("nobs_srv1_length");
  yrs_srv1_length.allocate(1,nobs_srv1_length,"yrs_srv1_length");
  nsamples_srv1_length.allocate(1,2,1,nobs_srv1_length,"nsamples_srv1_length");
  nobs_srv1_age.allocate("nobs_srv1_age");
  yrs_srv1_age.allocate(1,nobs_srv1_age,"yrs_srv1_age");
  like_wght.allocate(1,5,"like_wght");
  nsamples_srv1_age.allocate(1,2,1,nobs_srv1_age,"nsamples_srv1_age");
  obs_p_srv1_len_r.allocate(1,2,1,nobs_srv1_length,1,nlen_r,"obs_p_srv1_len_r");
  obs_p_srv1_age_read.allocate(1,2,1,nobs_srv1_age,1,nages_read,"obs_p_srv1_age_read");
  obs_p_fish_r.allocate(1,2,1,nobs_fish,1,nlen_r,"obs_p_fish_r");
  catch_bio.allocate(styr,endyr,"catch_bio");
  obs_srv1.allocate(1,nobs_srv1,"obs_srv1");
  obs_srv1_sd.allocate(1,nobs_srv1,"obs_srv1_sd");
  wt.allocate(1,2,1,nages,"wt");
  maturity.allocate(1,nages,"maturity");
  lenage.allocate(1,2,1,nages,1,nlen_r-1,"lenage");
  offset_const.allocate("offset_const");
  cv_srv1.allocate(1,nobs_srv1);
  nlen=nlen_r-1;
  styr_rec=styr-nages+1;
  if(nselages>nages) nselages=nages;
  if(nselages_srv1>nages) nselages_srv1=nages;
  cv_srv1=elem_div(obs_srv1_sd,obs_srv1);
  wt=wt*.001;
  obs_p_srv1_age.allocate(1,2,1,nobs_srv1_age,1,nages);
  obs_p_srv1_age_r.allocate(1,2,1,nobs_srv1_age,1,nages);
  obs_p_srv1_length.allocate(1,2,1,nobs_srv1_length,1,nlen);
  obs_p_fish.allocate(1,2,1,nobs_fish,1,nlen);
  obs_sexr.allocate(1,nobs_fish);
  obs_sexr_srv1.allocate(1,nobs_srv1_age);
  obs_sexr_srv1_l.allocate(1,nobs_srv1_length);
}

void model_parameters::initializationfunction(void)
{
  mean_log_rec.set_initial_value(18);
  log_avg_fmort.set_initial_value(-5.);
  q1.set_initial_value(1.);
  fmort_dev.set_initial_value(0.00001);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  q1.allocate(0.01,20.0,-8,"q1");
  mean_log_rec.allocate(1,"mean_log_rec");
  rec_dev.allocate(styr_rec,endyr,-15,15,2,"rec_dev");
  log_avg_fmort.allocate(1,"log_avg_fmort");
  fmort_dev.allocate(styr,endyr,-5.0,3.5,1,"fmort_dev");
  log_selcoffs_fish.allocate(1,2,1,nselages,phase_selcoffs,"log_selcoffs_fish");
  log_selcoffs_srv1.allocate(1,2,1,nselages_srv1,phase_selcoffs_srv1,"log_selcoffs_srv1");
  fish_slope_f.allocate(.1,5.,phase_logistic_sel,"fish_slope_f");
  fish_sel50_f.allocate(1.,8.,phase_logistic_sel,"fish_sel50_f");
  fish_slope_m.allocate(.05,2.0,phase_logistic_sel,"fish_slope_m");
  fish_sel50_m.allocate(1.,25.,phase_logistic_sel,"fish_sel50_m");
  srv1_slope_f.allocate(.1,5.,phase_logistic_sel_srv1,"srv1_slope_f");
  srv1_sel50_f.allocate(1.,10.,phase_logistic_sel_srv1,"srv1_sel50_f");
  srv1_slope_m.allocate(.01,5.,phase_logistic_sel_srv1,"srv1_slope_m");
  srv1_sel50_m.allocate(1.,10.,phase_logistic_sel_srv1,"srv1_sel50_m");
  sexr_param_fish.allocate(1.0,1.0,-5,"sexr_param_fish");
  sexr_param_srv.allocate(.25,1.0,phase_selcoffs_srv1,"sexr_param_srv");
  F40.allocate(0.05,.5,phase_F40,"F40");
  F35.allocate(0.05,.5,phase_F40,"F35");
  log_sel_fish.allocate(1,2,1,nages,"log_sel_fish");
  #ifndef NO_AD_INITIALIZE
    log_sel_fish.initialize();
  #endif
  log_sel_srv1.allocate(1,2,1,nages,"log_sel_srv1");
  #ifndef NO_AD_INITIALIZE
    log_sel_srv1.initialize();
  #endif
  sel.allocate(1,2,1,nages,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  sel_srv1.allocate(1,2,1,nages,"sel_srv1");
  #ifndef NO_AD_INITIALIZE
    sel_srv1.initialize();
  #endif
  avgsel_fish.allocate(1,2,"avgsel_fish");
  #ifndef NO_AD_INITIALIZE
    avgsel_fish.initialize();
  #endif
  avgsel_srv1.allocate(1,2,"avgsel_srv1");
  #ifndef NO_AD_INITIALIZE
    avgsel_srv1.initialize();
  #endif
  popn.allocate(1,2,styr,endyr,"popn");
  #ifndef NO_AD_INITIALIZE
    popn.initialize();
  #endif
  totn_srv1.allocate(1,2,styr,endyr,"totn_srv1");
  #ifndef NO_AD_INITIALIZE
    totn_srv1.initialize();
  #endif
  M.allocate(1,2,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  explbiom.allocate(styr,endyr,"explbiom");
  #ifndef NO_AD_INITIALIZE
    explbiom.initialize();
  #endif
  pred_bio.allocate(styr,endyr,"pred_bio");
  #ifndef NO_AD_INITIALIZE
    pred_bio.initialize();
  #endif
  fspbio.allocate(styr,endyr,"fspbio");
  #ifndef NO_AD_INITIALIZE
    fspbio.initialize();
  #endif
  pred_srv1.allocate(styr,endyr,"pred_srv1");
  #ifndef NO_AD_INITIALIZE
    pred_srv1.initialize();
  #endif
  pred_p_fish.allocate(1,2,styr,endyr,1,nlen,"pred_p_fish");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish.initialize();
  #endif
  pred_p_srv1_age.allocate(1,2,styr,endyr,1,nages,"pred_p_srv1_age");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv1_age.initialize();
  #endif
  pred_p_srv1_len.allocate(1,2,styr,endyr,1,nlen,"pred_p_srv1_len");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv1_len.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  natage.allocate(1,2,styr,endyr,1,nages,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  catage.allocate(1,2,styr,endyr,1,nages,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  natlength.allocate(1,2,styr,endyr,1,nlen,"natlength");
  #ifndef NO_AD_INITIALIZE
    natlength.initialize();
  #endif
  Z.allocate(1,2,styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,2,styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,2,styr,endyr,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  fmort.allocate(styr,endyr,"fmort");
  #ifndef NO_AD_INITIALIZE
    fmort.initialize();
  #endif
  rbar.allocate("rbar");
  #ifndef NO_AD_INITIALIZE
  rbar.initialize();
  #endif
  surv.allocate(1,2,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  offset.allocate(1,3,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  catch_like.allocate("catch_like");
  #ifndef NO_AD_INITIALIZE
  catch_like.initialize();
  #endif
  age_like.allocate(1,3,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  sel_like.allocate(1,4,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  surv_like.allocate("surv_like");
  #ifndef NO_AD_INITIALIZE
  surv_like.initialize();
  #endif
  recruits.allocate(styr,endyr,"recruits");
  biomassrep.allocate(styr,endyr,"biomassrep");
  fspbiorep.allocate(styr,endyr,"fspbiorep");
  endbiom.allocate("endbiom");
  depletion.allocate("depletion");
  f.allocate("f");
  tmp.allocate("tmp");
  #ifndef NO_AD_INITIALIZE
  tmp.initialize();
  #endif
  pred_sexr.allocate(styr,endyr,"pred_sexr");
  #ifndef NO_AD_INITIALIZE
    pred_sexr.initialize();
  #endif
  preds_sexr.allocate(styr,endyr,"preds_sexr");
  #ifndef NO_AD_INITIALIZE
    preds_sexr.initialize();
  #endif
  sigmar.allocate("sigmar");
  #ifndef NO_AD_INITIALIZE
  sigmar.initialize();
  #endif
  ftmp.allocate("ftmp");
  #ifndef NO_AD_INITIALIZE
  ftmp.initialize();
  #endif
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,3,1,nages,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  nage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"nage_future");
  #ifndef NO_AD_INITIALIZE
    nage_future.initialize();
  #endif
  F_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"F_future");
  #ifndef NO_AD_INITIALIZE
    F_future.initialize();
  #endif
  Z_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"Z_future");
  #ifndef NO_AD_INITIALIZE
    Z_future.initialize();
  #endif
  S_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"S_future");
  #ifndef NO_AD_INITIALIZE
    S_future.initialize();
  #endif
  catage_future.allocate(1,2,styr_fut,endyr_fut,1,nages,"catage_future");
  #ifndef NO_AD_INITIALIZE
    catage_future.initialize();
  #endif
  avg_rec_dev_future.allocate("avg_rec_dev_future");
  #ifndef NO_AD_INITIALIZE
  avg_rec_dev_future.initialize();
  #endif
  avg_F_future.allocate(1,3,"avg_F_future");
  #ifndef NO_AD_INITIALIZE
    avg_F_future.initialize();
  #endif
  catch_future.allocate(1,5,styr_fut,endyr_fut,"catch_future");
  #ifndef NO_AD_INITIALIZE
    catch_future.initialize();
  #endif
  fspbiom_fut.allocate(1,5,styr_fut,endyr_fut,"fspbiom_fut");
  future_biomass.allocate(1,5,styr_fut,endyr_fut,"future_biomass");
  explbiom_fut.allocate(styr_fut,endyr_fut,"explbiom_fut");
  #ifndef NO_AD_INITIALIZE
    explbiom_fut.initialize();
  #endif
  maxsel_fish.allocate("maxsel_fish");
  #ifndef NO_AD_INITIALIZE
  maxsel_fish.initialize();
  #endif
  maxsel_srv1.allocate(1,2,"maxsel_srv1");
  #ifndef NO_AD_INITIALIZE
    maxsel_srv1.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  B40.allocate("B40");
  #ifndef NO_AD_INITIALIZE
  B40.initialize();
  #endif
  B35.allocate("B35");
  #ifndef NO_AD_INITIALIZE
  B35.initialize();
  #endif
  AMeanRec.allocate("AMeanRec");
  #ifndef NO_AD_INITIALIZE
  AMeanRec.initialize();
  #endif
  like_natm.allocate("like_natm");
  #ifndef NO_AD_INITIALIZE
  like_natm.initialize();
  #endif
  like_q.allocate("like_q");
  #ifndef NO_AD_INITIALIZE
  like_q.initialize();
  #endif
 cout<<"wt_like"<<endl;
 cout<< wt_like <<endl;
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
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
      for(i=1; i<=nobs_fish;i++)
        {
         obs_sexr(i)=(sum(obs_p_fish_r(1,i)(2,nlen+1)))/(sum(obs_p_fish_r(1,i)(2,nlen+1))+sum(obs_p_fish_r(2,i)(2,nlen+1)));
        }
      for(i=1; i<=nobs_srv1_age;i++)
      {
       obs_sexr_srv1(i)=(sum(obs_p_srv1_age_r(1,i)))/(sum(obs_p_srv1_age_r(1,i))+sum(obs_p_srv1_age_r(2,i)));
      }
      for(i=1; i<=nobs_srv1_length;i++)
      {
       obs_sexr_srv1_l(i)=(sum(obs_p_srv1_len_r(1,i)(2,nlen+1)))/(sum(obs_p_srv1_len_r(1,i)(2,nlen+1))+sum(obs_p_srv1_len_r(2,i)(2,nlen+1)));
      }
 offset=0.; 
   for(k=1; k<=2; k++)
    {
       for (i=1; i <= nobs_fish; i++)
       {
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
         for (j=1; j<=nages; j++)
         {
            if (obs_p_srv1_age(k,i,j)>0.0)
             {
               offset(3)-=nsamples_srv1_age(k,i)*obs_p_srv1_age(k,i,j)*log(obs_p_srv1_age(k,i,j));
             }
         }
       }
    }
  M(1)=0.20;
  M(2)=0.35;
}

void model_parameters::userfunction(void)
{
  ofstream& evalout= *pad_evalout;
   get_selectivity();
   get_mortality();
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));
   get_numbers_at_age();
   get_catch_at_age();
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
}

void model_parameters::get_selectivity(void)
{
  ofstream& evalout= *pad_evalout;
  if(active(log_selcoffs_fish))
  {
   phase_logistic_sel=-2;
    for(k=1;k<=2;k++)
    {
      for (j=1;j<=nselages;j++)
      {
        log_sel_fish(k,j)=log_selcoffs_fish(k,j);
      }
    }
   //sets selectivity of ages older than nselages to selectivity at nselages  
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
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
      if(k==2)
       {
         sel(k)=sel(k)*sexr_param_fish;
       }
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
  if(active(log_selcoffs_srv1))
  {
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
}

void model_parameters::get_mortality(void)
{
  ofstream& evalout= *pad_evalout;
  maxsel_fish=max(sel(1));
  if(maxsel_fish<max(sel(2)))
      maxsel_fish=max(sel(2));
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
}

void model_parameters::get_numbers_at_age(void)
{
  ofstream& evalout= *pad_evalout;
  if(maxsel_srv1(1)<max(sel_srv1(2)))
     {
      maxsel_srv1(1)=max(sel_srv1(2)); 
      maxsel_srv1(2)=maxsel_srv1(1);
     }
  int itmp;
 //calc initial population  
  for (j=1;j<nages;j++)
    {
      itmp=styr+1-j;
      natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
      natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));
    }
    itmp=styr+1-nages;
    natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
    natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));
 // Now do for next several years----------------------------------
  for (i=styr+1;i<=endyr;i++)
    {
      natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
      natage(2,i,1)=natage(1,i,1);
    }
 //numbers at age
  for(k=1;k<=2;k++)
    {
      for (i=styr;i< endyr;i++)
       {
         //subvector - avoids writing a j loop  =++ increments the right side 
         //(1,nages-1) to 1+1 to nages-1+1 then does the assignment x(i)(1,n) 
         //takes the ith row of x the columns 1 to n
        for(j=1;j<nages;j++)
          {
            natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
          }
         //accumulates oldest ages
         natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
         popn(k,i)= natage(k,i)*sel(k);
     }
     popn(k,endyr)=natage(k,endyr)*sel(k);
    }
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
       explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
       pred_bio(i)+=natage(k,i)*wt(k);
       pred_p_srv1_len(k,i)=elem_prod(sel_srv1(k),natage(k,i))*lenage(k)/(totn_srv1(1,i)+totn_srv1(2,i));
       pred_p_srv1_age(k,i)=elem_prod(sel_srv1(k),natage(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      }
    }
    biomassrep=pred_bio;
    fspbiorep=fspbio;
    for(i=styr;i<=endyr;i++)
    {
       recruits(i)=mfexp(mean_log_rec+rec_dev(i));
    } 
    depletion=pred_bio(endyr)/pred_bio(styr);
    endbiom=pred_bio(endyr);
}

void model_parameters::get_catch_at_age(void)
{
  ofstream& evalout= *pad_evalout;
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
      }
   }
}

void model_parameters::Future_projections(void)
{
  ofstream& evalout= *pad_evalout;
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
        case 5:
          ftmp=mean(mfexp(log_avg_fmort+fmort_dev(endyr-4,endyr)));
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
}

void model_parameters::compute_spr_rates(void)
{
  ofstream& evalout= *pad_evalout;
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0=0.;
  SBF40=0.;
  SBF35=0.;
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
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*mfexp(-1.*M(1))/(1.-mfexp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/(1.-mfexp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/(1.-mfexp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  for (j=1;j<=nages;j++)
  {
   // Kill them off till- use fraction of the year e.g. april=0.25 - not for atf they spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*mfexp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
   }
  sprpen    = 300.*square((SBF40/SB0)-0.4);
  sprpen   += 300.*square((SBF35/SB0)-0.35);
  B0  = AMeanRec*SB0  ;
  B40 = AMeanRec*SBF40  ;
  B35 = AMeanRec*SBF35  ;
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& evalout= *pad_evalout;
  age_like=0.; //just having the . as part of the number identifies it as a floating point number (a double)
  sel_like=0.;
  fpen=0.;
  rec_like=0.;
  surv_like=0.;
  catch_like=0.;
  f=0.;
  like_natm=0.;
  like_q=0.;
 if (active(rec_dev))
   {
    age_like=0.;
    int ii;
    //recruitment likelihood - norm2 is sum of square values   
    rec_like=norm2(rec_dev(styr_rec,endyr));
    f+=rec_like;
     for(k=1;k<=2;k++)
     {
        for (i=1; i <= nobs_fish; i++)
        {
          ii=yrs_fish(i);
          //fishery length likelihood 
           for (j=1; j<=nlen; j++)
           {
             age_like(1)-=nsamples_fish(k,i)*(offset_const+obs_p_fish(k,i,j))*log(pred_p_fish(k,ii,j)+offset_const);
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
               age_like(2)-=nsamples_srv1_length(k,i)*(offset_const+obs_p_srv1_length(k,i,j))*log(pred_p_srv1_len(k,ii,j)+offset_const);
             }
         }
       }
     age_like(2)-=offset(2);
  //bracket for active(fish_slope_f) 
  //survey ages
         for(k=1;k<=2;k++)
         {
           for (i=1; i <=nobs_srv1_age; i++)
           {
              ii=yrs_srv1_age(i);
              //survey likelihood 
              for (j=1; j<=nages; j++)
              {
                age_like(3)-=nsamples_srv1_age(k,i)*(offset_const+obs_p_srv1_age(k,i,j))*log(pred_p_srv1_age(k,ii,j)+offset_const);
              }
            }
         } 
      age_like(3)-=offset(3);
      f+=like_wght(1)*age_like(1);     
      f+=like_wght(2)*age_like(2);
      f+=like_wght(3)*age_like(3);
   //end of if(active (rec_dev))
    }
  // Fit to indices (lognormal) 
  //weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 
   surv_like = norm2(elem_div(log(obs_srv1+offset_const)-log(pred_srv1(yrs_srv1)+offset_const),sqrt(2)*sqrt(log(elem_prod(cv_srv1,cv_srv1)+1.0))));
 //this subtracts the log(sd) from the likelihood - is a constant so I'm not adding it.
    catch_like=norm2(log(catch_bio+offset_const)-log(pred_catch+offset_const));
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
    f+=like_wght(5)*surv_like;  
    f+=like_wght(4)*catch_like;
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{4000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::final_calcs()
{
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;       
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(300);
  time(&start);
  R_out.open("R_input.txt");
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
