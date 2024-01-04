% This code reproduce the granger causality analysis in the direction
% amygdala to hippocampus and hippocampus to amygdala. The code reproduces
% Fig 2

clear all
close all
clc

restoredefaultpath
addpath '/Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox/fieldtrip-20210212/'
ft_defaults

mlist={'s6 right', 's13', 's15', 's16', 's16 right', 's25', 's32', 's33', 'z6', 'z6 right', 'z10', 'z1 right', 'z2','s36','s38','sz8'};%;
subjects = [60 13 15 16 160 25 32 33 6 600 10 100 2 36 38 8]
behavlist = {'Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient32','PatientZ_10','Patient33','Patient36','Patient38','PatientZ8'};

cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/0000_IAPS_Rec/DataRipples_Connectivity
load 'DataCleanAH_recognitionN16.mat'

%%
for subj = 1:length(dataeHit_AH);
    
  if subj == 8
  cfg= [];
  dataemoold{subj} = ft_appenddata(cfg, dataeHit_AH{subj},  dataeHitK_AH{subj}, dataeMiss_AH{subj});
  dataneuold{subj} = ft_appenddata(cfg, datanHit_AH{subj},  datanHitK_AH{subj}, datanMiss_AH{subj});
  dataemonew{subj} = ft_appenddata(cfg,dataeFAlR_AH{subj}, dataeFAlK_AH{subj}, dataeCrj_AH{subj});
  dataneunew{subj} = ft_appenddata(cfg,datanFAlK_AH{subj}, datanCrj_AH{subj});
  dataemo{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeFAlR_AH{subj}, dataeFAlK_AH{subj}, dataeCrj_AH{subj}, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  dataneu{subj} = ft_appenddata(cfg, datanHit_AH{subj}, datanFAlK_AH{subj}, datanCrj_AH{subj}, datanHitK_AH{subj}, datanMiss_AH{subj});
  
  data_eoldhit{subj}=dataeHit_AH{subj};
  data_noldhit{subj}=datanHit_AH{subj};
  cfg=[]
  data_ehitkmiss{subj}= ft_appenddata(cfg, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  data_nhitkmiss{subj}= ft_appenddata(cfg, datanHitK_AH{subj}, datanMiss_AH{subj});
  
  elseif subj == 11
  cfg= [];
  dataemoold{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  dataneuold{subj} = ft_appenddata(cfg, datanHit_AH{subj},datanMiss_AH{subj});
  dataemonew{subj} = ft_appenddata(cfg,dataeFAlR_AH{subj},dataeCrj_AH{subj});
  dataneunew{subj} = ft_appenddata(cfg,datanFAlK_AH{subj}, datanCrj_AH{subj});
  dataemo{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeHitK_AH{subj}, dataeMiss_AH{subj},dataeFAlR_AH{subj},dataeCrj_AH{subj});
  dataneu{subj} = ft_appenddata(cfg,  datanFAlK_AH{subj}, datanCrj_AH{subj},datanHit_AH{subj},datanMiss_AH{subj});
  
  data_eoldhit{subj}= dataeHit_AH{subj};
  data_noldhit{subj}= datanHit_AH{subj};
  cfg=[]
  data_ehitkmiss{subj}= ft_appenddata(cfg, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  data_nhitkmiss{subj}= datanMiss_AH{subj};
    
  elseif subj == 12
  cfg= [];
  dataemoold{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeHitK_AH{subj});
  dataneuold{subj} = ft_appenddata(cfg, datanHit_AH{subj}, datanHitK_AH{subj},datanMiss_AH{subj});
  dataemonew{subj} = ft_appenddata(cfg, dataeFAlR_AH{subj}, dataeCrj_AH{subj},dataeFAlK_AH{subj});
  dataneunew{subj} = ft_appenddata(cfg, datanFAlR_AH{subj}, datanFAlK_AH{subj}, datanCrj_AH{subj});
  dataemo{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeFAlR_AH{subj}, dataeHitK_AH{subj}, dataeCrj_AH{subj},dataeFAlK_AH{subj});
  dataneu{subj} = ft_appenddata(cfg, datanHit_AH{subj}, datanFAlR_AH{subj}, datanFAlK_AH{subj}, datanCrj_AH{subj},datanHitK_AH{subj},datanMiss_AH{subj});
  
  data_eoldhit{subj}= dataeHit_AH{subj};
  data_noldhit{subj}= datanHit_AH{subj};
  cfg=[]
  data_ehitkmiss{subj}= dataeHitK_AH{subj}
  data_nhitkmiss{subj}= ft_appenddata(cfg, datanHitK_AH{subj}, datanMiss_AH{subj});
    
  elseif subj == 13
  cfg= [];
  dataemoold{subj} = dataeHit_AH{subj};
  dataneuold{subj} = ft_appenddata(cfg, datanHit_AH{subj},datanMiss_AH{subj})
  dataemonew{subj} = dataeCrj_AH{subj}
  dataneunew{subj} = ft_appenddata(cfg, datanFAlR_AH{subj},datanCrj_AH{subj})
  dataemo{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeCrj_AH{subj});
  dataneu{subj} = ft_appenddata(cfg, datanHit_AH{subj},datanMiss_AH{subj},datanFAlR_AH{subj},datanCrj_AH{subj})
  
  % int cannot be compute
  
  else
  cfg= [];
  dataemo{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeHitK_AH{subj}, dataeMiss_AH{subj},dataeFAlR_AH{subj},dataeFAlK_AH{subj},dataeCrj_AH{subj});
  dataneu{subj} = ft_appenddata(cfg,  datanHit_AH{subj},datanHitK_AH{subj},datanMiss_AH{subj}, datanFAlR_AH{subj},datanFAlK_AH{subj}, datanCrj_AH{subj});
  dataemoold{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  dataneuold{subj} = ft_appenddata(cfg, datanHit_AH{subj},datanHitK_AH{subj},datanMiss_AH{subj});
  dataemonew{subj} = ft_appenddata(cfg, dataeFAlR_AH{subj},dataeFAlK_AH{subj},dataeCrj_AH{subj});
  dataneunew{subj} = ft_appenddata(cfg, datanFAlR_AH{subj},datanFAlK_AH{subj},datanCrj_AH{subj});
  
  data_eoldhit{subj}=dataeHit_AH{subj};
  data_noldhit{subj}=datanHit_AH{subj};
  cfg=[]
  data_ehitkmiss{subj}= ft_appenddata(cfg, dataeHitK_AH{subj}, dataeMiss_AH{subj});
  data_nhitkmiss{subj}= ft_appenddata(cfg, datanHitK_AH{subj}, datanMiss_AH{subj});
     
  end
end

%%
for subj=1:length(dataeHit_AH);
    if subj == 12;
    cfg= [];
    datahit{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, datanHit_AH{subj});
    datacrj{subj} = ft_appenddata(cfg, dataeCrj_AH{subj}, datanCrj_AH{subj});
    datamiss{subj} = datanMiss_AH{subj};
    else
    cfg= [];
    datahit{subj} = ft_appenddata(cfg, dataeHit_AH{subj}, datanHit_AH{subj});
    datacrj{subj} = ft_appenddata(cfg, dataeCrj_AH{subj}, datanCrj_AH{subj});
    datamiss{subj} = ft_appenddata(cfg, dataeMiss_AH{subj}, datanMiss_AH{subj});
end
end

%%
for subj=1:length(dataeHit_AH);
    if subj==11
    cfg= [];
    datahitkmiss{subj} = ft_appenddata(cfg, dataeHitK_AH{subj}, dataeMiss_AH{subj}, datanMiss_AH{subj});
    elseif subj==12
    cfg= [];
    datahitkmiss{subj} = ft_appenddata(cfg, dataeHitK_AH{subj}, datanHitK_AH{subj}, datanMiss_AH{subj});
    elseif subj==13
    cfg= [];
    datahitkmiss{subj} = datamiss{subj}
    else
    cfg= [];
    datahitkmiss{subj} = ft_appenddata(cfg, dataeHitK_AH{subj}, datanHitK_AH{subj}, dataeMiss_AH{subj}, datanMiss_AH{subj});
end
end
%% select the most lateral channels
%s60
Amy_list{1}=[2];
HC_list{1}=[3];
%s13
Amy_list{2}=[2];
HC_list{2}=[4];
%s15
Amy_list{3}=[2];
HC_list{3}=[3];
%s16
Amy_list{4}=[1];
HC_list{4}=[3];
%s160
Amy_list{5}=[1];
HC_list{5}=[3];
%s25
Amy_list{6}=[1];
HC_list{6}=[2];
%s32
Amy_list{7}=[1];
HC_list{7}=[2];
%s33
Amy_list{8}=[3];
HC_list{8}=[4];
%sz6
Amy_list{9}=[1];
HC_list{9}=[2];
%sz600
Amy_list{10}=[1];
HC_list{10}=[2];
%sz10
Amy_list{11}=[2];
HC_list{11}=[3];
%sz1
Amy_list{12}=[1];
HC_list{12}=[2];
%sz2
Amy_list{13}=[1];
HC_list{13}=[2];
%s36
Amy_list{14}=[1];
HC_list{14}=[3];
%s38
Amy_list{15}=[2];
HC_list{15}=[4];
%sz8
Amy_list{16}=[2];
HC_list{16}=[3];
%%
for subj = 1:16;
     
    time_beginning=1;
    time_end=size(dataemo{subj}.time{1},2);

    
    for trial=1:size(dataemo{subj}.trial,2)
        dataemo_alltrials{subj}(trial,:,:)=dataemo{subj}.trial{trial};
    end
    
    for trial=1:size(dataneu{subj}.trial,2)
        dataneu_alltrials{subj}(trial,:,:)=dataneu{subj}.trial{trial};
    end
    
    for trial=1:size(dataemoold{subj}.trial,2)
        dataemoold_alltrials{subj}(trial,:,:)=dataemoold{subj}.trial{trial};
    end
    
    for trial=1:size(dataneuold{subj}.trial,2)
        dataneuold_alltrials{subj}(trial,:,:)=dataneuold{subj}.trial{trial};
    end
    
    for trial=1:size(dataemonew{subj}.trial,2)
        dataemonew_alltrials{subj}(trial,:,:)=dataemonew{subj}.trial{trial};
    end
    
    for trial=1:size(dataneunew{subj}.trial,2)
        dataneunew_alltrials{subj}(trial,:,:)=dataneunew{subj}.trial{trial};
    end
    
    for trial=1:size(datahit{subj}.trial,2)
        datahit_alltrials{subj}(trial,:,:)=datahit{subj}.trial{trial};
    end
    
    for trial=1:size(datahitkmiss{subj}.trial,2)
        datahitkmiss_alltrials{subj}(trial,:,:)=datahitkmiss{subj}.trial{trial};
    end
  %%  
    
    for trial=1:size(dataemo{subj}.trial,2)
        for chan=1:size(dataemo{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataemo{subj}.trial{trial}(chan,timepoint)=dataemo{subj}.trial{trial}(chan,timepoint)-mean(dataemo_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    for trial=1:size(dataneu{subj}.trial,2)
        for chan=1:size(dataneu{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataneu{subj}.trial{trial}(chan,timepoint)=dataneu{subj}.trial{trial}(chan,timepoint)-mean(dataneu_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(dataemoold{subj}.trial,2)
        for chan=1:size(dataemoold{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataemoold{subj}.trial{trial}(chan,timepoint)=dataemoold{subj}.trial{trial}(chan,timepoint)-mean(dataemoold_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    for trial=1:size(dataneuold{subj}.trial,2)
        for chan=1:size(dataneuold{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataneuold{subj}.trial{trial}(chan,timepoint)=dataneuold{subj}.trial{trial}(chan,timepoint)-mean(dataneuold_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(dataemonew{subj}.trial,2)
        for chan=1:size(dataemonew{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataemonew{subj}.trial{trial}(chan,timepoint)=dataemonew{subj}.trial{trial}(chan,timepoint)-mean(dataemonew_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    for trial=1:size(dataneunew{subj}.trial,2)
        for chan=1:size(dataneunew{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                dataneunew{subj}.trial{trial}(chan,timepoint)=dataneunew{subj}.trial{trial}(chan,timepoint)-mean(dataneunew_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(datahit{subj}.trial,2)
        for chan=1:size(datahit{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                datahit{subj}.trial{trial}(chan,timepoint)=datahit{subj}.trial{trial}(chan,timepoint)-mean(datahit_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    
    for trial=1:size(datahitkmiss{subj}.trial,2)
        for chan=1:size(datahitkmiss{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                datahitkmiss{subj}.trial{trial}(chan,timepoint)=datahitkmiss{subj}.trial{trial}(chan,timepoint)-mean(datahitkmiss_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
%%
    
    %II.III Preprocessing: 
    %-----------------------------------------------------------------------
    cfg = [];
    cfg.dftfilter ='yes';
    cfg.lpfilter = 'yes'
    cfg.lpfreq = 85
    preprocessing_emo{subj} = ft_preprocessing(cfg,dataemo{subj});
    preprocessing_neu{subj} = ft_preprocessing(cfg,dataneu{subj});
    preprocessing_emoold{subj} = ft_preprocessing(cfg,dataemoold{subj});
    preprocessing_neuold{subj} = ft_preprocessing(cfg,dataneuold{subj});
    preprocessing_emonew{subj} = ft_preprocessing(cfg,dataemonew{subj});
    preprocessing_neunew{subj} = ft_preprocessing(cfg,dataneunew{subj});
    
    preprocessing_hit{subj} = ft_preprocessing(cfg,datahit{subj});
    preprocessing_hitkmiss{subj} = ft_preprocessing(cfg,datahitkmiss{subj});

    

    
    % resample
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'no';
    preprocessing_emo_ds{subj}  = ft_resampledata(cfg,preprocessing_emo{subj});
    preprocessing_neu_ds{subj}  = ft_resampledata(cfg,preprocessing_neu{subj});
    preprocessing_emoold_ds{subj}  = ft_resampledata(cfg,preprocessing_emoold{subj});
    preprocessing_neuold_ds{subj}  = ft_resampledata(cfg,preprocessing_neuold{subj});
    preprocessing_emonew_ds{subj}  = ft_resampledata(cfg,preprocessing_emonew{subj});
    preprocessing_neunew_ds{subj}  = ft_resampledata(cfg,preprocessing_neunew{subj});
    
    preprocessing_hit_ds{subj}  = ft_resampledata(cfg,preprocessing_hit{subj});
    preprocessing_hitkmiss_ds{subj}  = ft_resampledata(cfg,preprocessing_hitkmiss{subj});

    
    cfg         = [];
    cfg.order   = 9;
    cfg.toolbox = 'bsmart';
    cfg.toi       = -0.5:0.01:1.5;
    cfg.t_ftimwin = 0.4;
    mdata_emo{subj} = ft_mvaranalysis(cfg, preprocessing_emo_ds{subj});
    mdata_neu{subj} = ft_mvaranalysis(cfg, preprocessing_neu_ds{subj});
    mdata_emoold{subj} = ft_mvaranalysis(cfg, preprocessing_emoold_ds{subj});
    mdata_neuold{subj} = ft_mvaranalysis(cfg, preprocessing_neuold_ds{subj});
    mdata_emonew{subj} = ft_mvaranalysis(cfg, preprocessing_emonew_ds{subj});
    mdata_neunew{subj} = ft_mvaranalysis(cfg, preprocessing_neunew_ds{subj});
    
    mdata_hit{subj} = ft_mvaranalysis(cfg, preprocessing_hit_ds{subj});
    mdata_hitkmiss{subj} = ft_mvaranalysis(cfg, preprocessing_hitkmiss_ds{subj});

   
        
    %     I.III Analyse mit Wavelet Convolution Methode
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    freq_choice=2:80;
    
    
    cfg = [];
    cfg.method = 'mvar';
    cfg.foi        = freq_choice;
    Freq_emo{subj} =  ft_freqanalysis_mvar(cfg,mdata_emo{subj});
    Freq_neu{subj} =  ft_freqanalysis_mvar(cfg,mdata_neu{subj});
    Freq_emoold{subj} =  ft_freqanalysis_mvar(cfg,mdata_emoold{subj});
    Freq_neuold{subj} =  ft_freqanalysis_mvar(cfg,mdata_neuold{subj});
    Freq_emonew{subj} =  ft_freqanalysis_mvar(cfg,mdata_emonew{subj});
    Freq_neunew{subj} =  ft_freqanalysis_mvar(cfg,mdata_neunew{subj});
    
    Freq_hit{subj} =  ft_freqanalysis_mvar(cfg,mdata_hit{subj});
    Freq_hitkmiss{subj} =  ft_freqanalysis_mvar(cfg,mdata_hitkmiss{subj});

    
    cfg           = [];
    cfg.method    = 'granger';
    emo_granger{subj}  = ft_connectivityanalysis(cfg, Freq_emo{subj});
    neu_granger{subj}  = ft_connectivityanalysis(cfg, Freq_neu{subj});
    emoold_granger{subj}  = ft_connectivityanalysis(cfg, Freq_emoold{subj});
    neuold_granger{subj}  = ft_connectivityanalysis(cfg, Freq_neuold{subj});
    emonew_granger{subj}  = ft_connectivityanalysis(cfg, Freq_emonew{subj});
    neunew_granger{subj}  = ft_connectivityanalysis(cfg, Freq_neunew{subj});
    
    hit_granger{subj}  = ft_connectivityanalysis(cfg, Freq_hit{subj});
    hitkmiss_granger{subj}  = ft_connectivityanalysis(cfg, Freq_hitkmiss{subj});


    
    % interaction
    GA_emo_HCtoAmyg(subj,:,:)= emo_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emo_AmygtoHC(subj,:,:)= emo_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_neu_HCtoAmyg(subj,:,:)= neu_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neu_AmygtoHC(subj,:,:)= neu_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_emoold_HCtoAmyg(subj,:,:)= emoold_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emoold_AmygtoHC(subj,:,:)= emoold_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_neuold_HCtoAmyg(subj,:,:)= neuold_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neuold_AmygtoHC(subj,:,:)= neuold_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_emonew_HCtoAmyg(subj,:,:)= emonew_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emonew_AmygtoHC(subj,:,:)= emonew_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_neunew_HCtoAmyg(subj,:,:)= neunew_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neunew_AmygtoHC(subj,:,:)= neunew_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_hit_HCtoAmyg(subj,:,:)= hit_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_hit_AmygtoHC(subj,:,:)= hit_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_hitkmiss_HCtoAmyg(subj,:,:)= hitkmiss_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_hitkmiss_AmygtoHC(subj,:,:)= hitkmiss_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);

end
%% now compute the interaction
subjects_renam=[60 13 15 16 160 25 32 33 6 600 10 11 2000 36 38 8]
method = 'oneel'
subjects_ers1 = [60 13 15 160 25 32 33 600 2000 36 38 8];
posint = find(ismember(subjects_renam,subjects_ers1))
%% subjects included in the int [60 13 15 160 25 32 33 600 2000 36 38 8]
%%
for subj= [1:3,5:8,10,14:16];
        
    %% substract ensemble mean
      
    time_beginning=1;
    time_end=size(data_eoldhit{subj}.time{1},2);
  
    for trial=1:size(data_eoldhit{subj}.trial,2)
        dataeoldhit_alltrials{subj}(trial,:,:)=data_eoldhit{subj}.trial{trial};
    end
    
    for trial=1:size(data_ehitkmiss{subj}.trial,2)
        dataehitkmiss_alltrials{subj}(trial,:,:)=data_ehitkmiss{subj}.trial{trial};
    end
    %%
    
    for trial=1:size(data_noldhit{subj}.trial,2)
        datanoldhit_alltrials{subj}(trial,:,:)=data_noldhit{subj}.trial{trial};
    end
    
    for trial=1:size(data_nhitkmiss{subj}.trial,2)
        datanhitkmiss_alltrials{subj}(trial,:,:)=data_nhitkmiss{subj}.trial{trial};
    end
 %%
    for trial=1:size(data_eoldhit{subj}.trial,2)
        for chan=1:size(data_eoldhit{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                if size(data_eoldhit{subj}.trial,2) ~= 1
                data_eoldhit{subj}.trial{trial}(chan,timepoint)=data_eoldhit{subj}.trial{trial}(chan,timepoint)-mean(dataeoldhit_alltrials{subj}(:,chan,timepoint),1);
                else
                data_eoldhit{subj}.trial{trial}(chan,timepoint)=data_eoldhit{subj}.trial{trial}(chan,timepoint);
                end
            end
        end
    end
    
    
    for trial=1:size(data_ehitkmiss{subj}.trial,2)
        for chan=1:size(data_ehitkmiss{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                if size(data_ehitkmiss{subj}.trial,2) ~= 1
                    data_ehitkmiss{subj}.trial{trial}(chan,timepoint)=data_ehitkmiss{subj}.trial{trial}(chan,timepoint)-mean(dataehitkmiss_alltrials{subj}(:,chan,timepoint),1);
                else
                    data_ehitkmiss{subj}.trial{trial}(chan,timepoint)=data_ehitkmiss{subj}.trial{trial}(chan,timepoint);
                end
            end
        end
    end
    %%
    
    for trial=1:size(data_noldhit{subj}.trial,2)
        for chan=1:size(data_noldhit{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                if size(data_noldhit{subj}.trial,2) ~= 1
                    data_noldhit{subj}.trial{trial}(chan,timepoint)=data_noldhit{subj}.trial{trial}(chan,timepoint)-mean(datanoldhit_alltrials{subj}(:,chan,timepoint),1);
                else
                    data_noldhit{subj}.trial{trial}(chan,timepoint)=data_noldhit{subj}.trial{trial}(chan,timepoint)
                end
            end
        end
    end
    
    
    for trial=1:size(data_nhitkmiss{subj}.trial,2)
        for chan=1:size(data_nhitkmiss{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                if size(data_nhitkmiss{subj}.trial,2) ~= 1
                    data_nhitkmiss{subj}.trial{trial}(chan,timepoint)=data_nhitkmiss{subj}.trial{trial}(chan,timepoint)-mean(datanhitkmiss_alltrials{subj}(:,chan,timepoint),1);
                else
                    data_nhitkmiss{subj}.trial{trial}(chan,timepoint)=data_nhitkmiss{subj}.trial{trial}(chan,timepoint);
                end
            end
        end
    end
%%
    
    %II.III Preprocessing: 
    %-----------------------------------------------------------------------
    cfg = [];
    cfg.dftfilter ='yes';
    cfg.lpfilter = 'yes'
    cfg.lpfreq = 85
    
    preprocessing_eoldhit{subj} = ft_preprocessing(cfg,data_eoldhit{subj});
    preprocessing_ehitkmiss{subj} = ft_preprocessing(cfg,data_ehitkmiss{subj});
    preprocessing_noldhit{subj} = ft_preprocessing(cfg,data_noldhit{subj});
    preprocessing_nhitkmiss{subj} = ft_preprocessing(cfg,data_nhitkmiss{subj});

    
    % resample
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'no';
    
    preprocessing_eoldhit_ds{subj}  = ft_resampledata(cfg,preprocessing_eoldhit{subj});
    preprocessing_ehitkmiss_ds{subj}  = ft_resampledata(cfg,preprocessing_ehitkmiss{subj});
    preprocessing_noldhit_ds{subj}  = ft_resampledata(cfg,preprocessing_noldhit{subj});
    preprocessing_nhitkmiss_ds{subj}  = ft_resampledata(cfg,preprocessing_nhitkmiss{subj});

 
    cfg         = [];
    cfg.order   = 9;
    cfg.toolbox = 'bsmart';
    cfg.toi       = -0.5:0.01:1.5;
    cfg.t_ftimwin = 0.4;
    %interaction
    
    mdata_eoldhit{subj} = ft_mvaranalysis(cfg, preprocessing_eoldhit_ds{subj});
    mdata_ehitkmiss{subj} = ft_mvaranalysis(cfg, preprocessing_ehitkmiss_ds{subj});
    mdata_noldhit{subj} = ft_mvaranalysis(cfg, preprocessing_noldhit_ds{subj});
    mdata_nhitkmiss{subj} = ft_mvaranalysis(cfg, preprocessing_nhitkmiss_ds{subj});

        
    %     I.III Analyse mit Wavelet Convolution Methode
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    freq_choice=2:80;
    
    
    cfg = [];
    cfg.method = 'mvar';
    cfg.foi        = freq_choice;
    %interaction

    Freq_eoldhit{subj} =  ft_freqanalysis_mvar(cfg,mdata_eoldhit{subj});
    Freq_ehitkmiss{subj} =  ft_freqanalysis_mvar(cfg,mdata_ehitkmiss{subj});
    
    Freq_noldhit{subj} =  ft_freqanalysis_mvar(cfg,mdata_noldhit{subj});
    Freq_nhitkmiss{subj} =  ft_freqanalysis_mvar(cfg,mdata_nhitkmiss{subj});

    
    cfg           = [];
    cfg.method    = 'granger';
    % interaction

    eoldhit_granger{subj}  = ft_connectivityanalysis(cfg, Freq_eoldhit{subj});
    noldhit_granger{subj}  = ft_connectivityanalysis(cfg, Freq_noldhit{subj});
    ehitkmiss_granger{subj}  = ft_connectivityanalysis(cfg, Freq_ehitkmiss{subj});
    nhitkmiss_granger{subj}  = ft_connectivityanalysis(cfg, Freq_nhitkmiss{subj});
       
    % interaction
    
    GA_eoldhit_HCtoAmyg(subj,:,:)= eoldhit_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_noldhit_HCtoAmyg(subj,:,:)= noldhit_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_ehitkmiss_HCtoAmyg(subj,:,:)= ehitkmiss_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_nhitkmiss_HCtoAmyg(subj,:,:)= nhitkmiss_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);

    
    GA_eoldhit_AmygtoHC(subj,:,:)= eoldhit_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_noldhit_AmygtoHC(subj,:,:)= noldhit_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_ehitkmiss_AmygtoHC(subj,:,:)= ehitkmiss_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_nhitkmiss_AmygtoHC(subj,:,:)= nhitkmiss_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
end
%% compute the GA
subjects_renam=[60 13 15 16 160 25 32 33 6 600 10 11 2000 36 38 8]
subjects_ers = [60 13 15 16 160 25 32 11 33 61 600 10 2000 36 38 8];

method = 'oneel'
if strcmp(method,'allel')
pos=find(ismember(subjects_renam,subjects_ers))
else
subjects_ers1 = [60 13 15 160 25 32 11 33 600 10 2000 36 38 8];
pos = find(ismember(subjects_renam,subjects_ers1))
end

posint=[1:3,5:8,10,14:16];
%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load Gatemp.mat

GTFs_Multitaper_baseline_Unpl.time = []
GTFs_Multitaper_baseline_Unpl.time = Freq_emo{1}.time
GTFs_Multitaper_baseline_Unpl.freq = []
GTFs_Multitaper_baseline_Unpl.freq = Freq_emo{1}.freq
GTFs_Multitaper_baseline_Unpl.label = []
GTFs_Multitaper_baseline_Unpl.label = {'AH'}


GA_stat_emoold_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_emoold_HCtoAmyg.powspctrm=[];
GA_stat_emoold_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_emoold_HCtoAmyg(pos,:,:,:),4));

GA_stat_neuold_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_neuold_HCtoAmyg.powspctrm=[];
GA_stat_neuold_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_neuold_HCtoAmyg(pos,:,:,:),4));

GA_stat_hit_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_hit_HCtoAmyg.powspctrm=[];
GA_stat_hit_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_hit_HCtoAmyg(pos,:,:,:),4));


GA_stat_hitkmiss_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_hitkmiss_HCtoAmyg.powspctrm=[];
GA_stat_hitkmiss_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_hitkmiss_HCtoAmyg(pos,:,:,:),4));
%%

GA_stat_emoold_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_emoold_AmygtoHC.powspctrm=[];
GA_stat_emoold_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_emoold_AmygtoHC(pos,:,:,:),4));

GA_stat_neuold_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_neuold_AmygtoHC.powspctrm=[];
GA_stat_neuold_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_neuold_AmygtoHC(pos,:,:,:),4));

GA_stat_hit_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_hit_AmygtoHC.powspctrm=[];
GA_stat_hit_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_hit_AmygtoHC(pos,:,:,:),4));

GA_stat_hitkmiss_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_hitkmiss_AmygtoHC.powspctrm=[];
GA_stat_hitkmiss_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_hitkmiss_AmygtoHC(pos,:,:,:),4));

%%

GA_stat_emoolda2h = GA_stat_emoold_AmygtoHC
GA_stat_emoolda2h.powspctrm = []
GA_stat_emoolda2h.powspctrm = GA_stat_emoold_AmygtoHC.powspctrm - GA_stat_neuold_AmygtoHC.powspctrm

GA_stat_emooldh2a = GA_stat_emoold_AmygtoHC
GA_stat_emooldh2a.powspctrm = []
GA_stat_emooldh2a.powspctrm = GA_stat_emoold_HCtoAmyg.powspctrm - GA_stat_neuold_HCtoAmyg.powspctrm

%%
GA_stat_memhitkmissa2h = GA_stat_emoold_AmygtoHC
GA_stat_memhitkmissa2h.powspctrm = []
GA_stat_memhitkmissa2h.powspctrm = GA_stat_hit_AmygtoHC.powspctrm - GA_stat_hitkmiss_AmygtoHC.powspctrm

GA_stat_memhitkmissh2a = GA_stat_emoold_AmygtoHC
GA_stat_memhitkmissh2a.powspctrm = []
GA_stat_memhitkmissh2a.powspctrm = GA_stat_hit_HCtoAmyg.powspctrm - GA_stat_hitkmiss_HCtoAmyg.powspctrm
%% cluter stat

cfg=[];
cfg.method = 'montecarlo'%
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=2;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 'all';

cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GA_stat_emoold_AmygtoHC.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

emoneuold_A2H = ft_freqstatistics(cfg, GA_stat_emoold_AmygtoHC, GA_stat_neuold_AmygtoHC);
emoneuold_H2A = ft_freqstatistics(cfg, GA_stat_emoold_HCtoAmyg, GA_stat_neuold_HCtoAmyg);

memhitkmiss_A2H = ft_freqstatistics(cfg, GA_stat_hit_AmygtoHC, GA_stat_hitkmiss_AmygtoHC);
memhitkmiss_H2A = ft_freqstatistics(cfg, GA_stat_hit_HCtoAmyg, GA_stat_hitkmiss_HCtoAmyg);

%%
GA_stat_eoldhit_a2h=GTFs_Multitaper_baseline_Unpl;
GA_stat_eoldhit_a2h.powspctrm=[];
GA_stat_eoldhit_a2h.powspctrm(:,1,:,:)=squeeze(mean(GA_eoldhit_AmygtoHC(posint,:,:,:),4));

GA_stat_noldhit_a2h=GTFs_Multitaper_baseline_Unpl;
GA_stat_noldhit_a2h.powspctrm=[];
GA_stat_noldhit_a2h.powspctrm(:,1,:,:)=squeeze(mean(GA_noldhit_AmygtoHC(posint,:,:,:),4));

GA_stat_ehitkmiss_a2h=GTFs_Multitaper_baseline_Unpl;
GA_stat_ehitkmiss_a2h.powspctrm=[];
GA_stat_ehitkmiss_a2h.powspctrm(:,1,:,:)=squeeze(mean(GA_ehitkmiss_AmygtoHC(posint,:,:,:),4));

GA_stat_nhitkmiss_a2h=GTFs_Multitaper_baseline_Unpl;
GA_stat_nhitkmiss_a2h.powspctrm=[];
GA_stat_nhitkmiss_a2h.powspctrm(:,1,:,:)=squeeze(mean(GA_nhitkmiss_AmygtoHC(posint,:,:,:),4));
%%
GA_stat_eoldhit_h2a=GTFs_Multitaper_baseline_Unpl;
GA_stat_eoldhit_h2a.powspctrm=[];
GA_stat_eoldhit_h2a.powspctrm(:,1,:,:)=squeeze(mean(GA_eoldhit_HCtoAmyg(posint,:,:,:),4));

GA_stat_noldhit_h2a=GTFs_Multitaper_baseline_Unpl;
GA_stat_noldhit_h2a.powspctrm=[];
GA_stat_noldhit_h2a.powspctrm(:,1,:,:)=squeeze(mean(GA_noldhit_HCtoAmyg(posint,:,:,:),4));

GA_stat_ehitkmiss_h2a=GTFs_Multitaper_baseline_Unpl;
GA_stat_ehitkmiss_h2a.powspctrm=[];
GA_stat_ehitkmiss_h2a.powspctrm(:,1,:,:)=squeeze(mean(GA_ehitkmiss_HCtoAmyg(posint,:,:,:),4));

GA_stat_nhitkmiss_h2a=GTFs_Multitaper_baseline_Unpl;
GA_stat_nhitkmiss_h2a.powspctrm=[];
GA_stat_nhitkmiss_h2a.powspctrm(:,1,:,:)=squeeze(mean(GA_nhitkmiss_HCtoAmyg(posint,:,:,:),4));
%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load Grangerdata.mat
load Grangerstat.mat
%% Reproduce Figure 2

h = figure;set(h,'position', [0 0 1000 400]) %this is in pixel
cfg = [];
cfg.figure = 'gcf';
cfg.xlim =[0 1.5]
cfg.ylim = [2 34]
cfg.zlim = [-4 4];
subplot(2,2,1)
logRelative1 = emoneuold_H2A
logRelative1.mask = emoneuold_H2A.mask;
logRelative1.powspctrm = emoneuold_H2A.stat% You shoud add a field mask with the mask from the statistics
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('emo vs neu H2A'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
%title(['emo vs neu h2a  p=' num2str([emoneuold_H2A.posclusters(1).prob])])
title('emotional vs neutral h2a');% p=' num2str([emoneuold_H2A.posclusters(1).prob])])

hold on
subplot(2,2,2)
E = emoneuold_A2H
%emoneu_A2H.mask=emoneu_A2H.posclusterslabelmat==1;
E.powspctrm = emoneuold_A2H.stat;%.*emoneuold_A2H.mask;
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('emo vs neu A2H'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
%title(['emo vs neu a2h  p=' num2str([emoneuold_A2H.posclusters(1).prob])])
title('emotional vs neutral a2h');% p=' num2str([emoneuold_H2A.posclusters(1).prob])])

subplot(2,2,3)
logRelative1 = memhitkmiss_H2A
logRelative1.mask = memhitkmiss_H2A.mask;
logRelative1.powspctrm = memhitkmiss_H2A.stat% You shoud add a field mask with the mask from the statistics
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('hit vs hitkmiss H2A'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
%title(['h2a  p=' num2str([memhitkmiss_H2A.posclusters(1).prob])])
title('hit vs hitkmiss h2a')

hold on
subplot(2,2,4)
logRelative1 = memhitkmiss_A2H
logRelative1.mask = memhitkmiss_A2H.mask;
logRelative1.powspctrm = memhitkmiss_A2H.stat% You shoud add a field mask with the mask from the statistics
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('hit vs hitkmiss A2H'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
%title(['ah2  p=' num2str([memhitkmiss_H2A.posclusters(1).prob])])
title('hit vs hitkmiss ah2')

%%
GA_stat_ediff_h2a= GA_stat_eoldhit_h2a
GA_stat_ediff_h2a.powspctrm = GA_stat_eoldhit_h2a.powspctrm - GA_stat_ehitkmiss_h2a.powspctrm

GA_stat_ndiff_h2a= GA_stat_noldhit_h2a
GA_stat_ndiff_h2a.powspctrm = GA_stat_noldhit_h2a.powspctrm - GA_stat_nhitkmiss_h2a.powspctrm

GA_stat_ediff_a2h= GA_stat_eoldhit_a2h
GA_stat_ediff_a2h.powspctrm = GA_stat_eoldhit_a2h.powspctrm - GA_stat_ehitkmiss_a2h.powspctrm

GA_stat_ndiff_a2h= GA_stat_noldhit_a2h
GA_stat_ndiff_a2h.powspctrm = GA_stat_noldhit_a2h.powspctrm - GA_stat_nhitkmiss_a2h.powspctrm
%%
% save ('Grangerdata.mat','GA_stat_emoold_AmygtoHC','GA_stat_neuold_AmygtoHC','GA_stat_emoold_HCtoAmyg','GA_stat_neuold_HCtoAmyg',...
%     'GA_stat_hit_AmygtoHC','GA_stat_hitkmiss_AmygtoHC','GA_stat_hit_HCtoAmyg','GA_stat_hitkmiss_HCtoAmyg',...
%     'GA_stat_eoldhit_h2a','GA_stat_ehitkmiss_h2a','GA_stat_noldhit_h2a','GA_stat_nhitkmiss_h2a',...
%     'GA_stat_eoldhit_a2h','GA_stat_ehitkmiss_a2h','GA_stat_noldhit_a2h','GA_stat_nhitkmiss_a2h')
%% cluter stat

cfg=[];
cfg.method = 'montecarlo'%
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=2;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 'all';

cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GA_stat_ediff_h2a.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

enhithitkmiss_H2A = ft_freqstatistics(cfg, GA_stat_ediff_h2a, GA_stat_ndiff_h2a);
enhithitkmiss_A2H = ft_freqstatistics(cfg, GA_stat_ediff_a2h, GA_stat_ndiff_a2h);
%%
save ('Grangerstat.mat','emoneuold_A2H','emoneuold_H2A','memhitkmiss_A2H','memhitkmiss_H2A');%,'enhithitkmiss_H2A','enhithitkmiss_A2H')