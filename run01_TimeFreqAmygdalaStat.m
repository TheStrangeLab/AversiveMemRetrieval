
clear all
close all
clc

restoredefaultpath
addpath '/Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox/fieldtrip-20210212/'
ft_defaults

cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load 'Ga_Amygdala_Timefreq.mat'

%% compute the diff for the interaction
GTFs_diff_eHitnHit = Ga_eHit
GTFs_diff_eHitnHit.powspctrm = Ga_eHit.powspctrm - Ga_nHit.powspctrm

GTFs_diff_eKMissnKMiss = Ga_eKMiss
GTFs_diff_eKMissnKMiss.powspctrm = Ga_eKMiss.powspctrm - Ga_nKMiss.powspctrm

%%
cfg=[];
cfg.method = 'montecarlo'%; %'analytic
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=35;
f2=150;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

%cfg.minnbchan      = 2;
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GTFs_diff_eKMissnKMiss.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Int_hitKmiss = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eKMissnKMiss);
%%
ns=size(Ga_Unpl.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Emo = ft_freqstatistics(cfg,Ga_Unpl, Ga_Neu);
%%
ns=size(Ga_KMiss.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
HitKMiss = ft_freqstatistics(cfg,Ga_Hit, Ga_KMiss);
%%
%cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
%save ('statTFAmygdala.mat','Emo','Int_hitKmiss','HitKMiss')
%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load statTFAmygdala.mat
%% plot data for Fig 1
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];%
cfg.ylim = [35 150];%'maxmin';
cfg.zlim = [-0.3 0.3];%[0 0.2];%'maxmin';
h = figure;%set(h,'position', [0 0 400 188])
subplot(2,2,1); ft_singleplotTFR(cfg,Ga_eHit); 
tit=title('eRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nHit);
tit=title('nRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3); ft_singleplotTFR(cfg,Ga_eKMiss);
tit=title('eKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKMiss);
tit=title('nKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on

%% plot unpleasant vs neutral stat for Fig 1 

h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg.zlim = [-4 4];
cfg.xlim = [0 1.5]
M = Emo;
M.powspctrm = Emo.stat.*Emo.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('emo'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16); 
%%
h = figure;%set(h,'position', [0 0 250 188]) %this is in pixel
logRelative1 = Emo
logRelative1.mask = Emo.mask;
logRelative1.powspctrm = Emo.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('emotional effect'); xlabel('time (sec)'); ylabel('frequency (Hz)');

%%
maskt= sum(squeeze(Emo.mask),1);
tvec = Emo.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(Emo.mask),2);
fvec = Emo.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t=toi
f=foi


pt1 = nearest(Ga_Unpl.time,t(1));
pt2 = nearest(Ga_Unpl.time,t(2));
pf1 = nearest(Ga_Unpl.freq,f(1));
pf2 = nearest(Ga_Unpl.freq,f(2));

clear mat

mat(:,1) = squeeze(mean(mean(Ga_Unpl.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_Neu.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
%%
addpath /Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox_plot/beeswarm-master
addpath /Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox_plot/stdshade.m
h = figure;set(h,'position', [0 0 250 188]) %this is in pixel
nsj=size(mat,1)
x = [ones(nsj,1) ones(nsj,1)*2];
y = [mat(:,1) mat(:,2)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0 0; 0 0 1])
ylim([-0.4 0.4]);
xticklabels({[],'Aversive',[],'Neutral'})
ylabel('mean gamma power')

%% Supplementary Figure 2
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];
cfg.ylim = [0 34];
cfg.zlim = [-1 1];
h = figure;%set(h,'position', [0 0 400 188])
subplot(2,2,1); ft_singleplotTFR(cfg,Ga_eHit); 
tit=title('eRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nHit);
tit=title('nRHit');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3); ft_singleplotTFR(cfg,Ga_eKMiss);
tit=title('eKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKMiss);
tit=title('nKHit&Miss');set(findobj(tit,'type','text'),'FontSize',36);
hold on
%%
cfg=[];
cfg.method = 'montecarlo'%; %'analytic
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=0;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

%cfg.minnbchan      = 2;
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(GTFs_diff_eKMissnKMiss.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Int_hitKmisslow = ft_freqstatistics(cfg,GTFs_diff_eHitnHit, GTFs_diff_eKMissnKMiss);
%%
ns=size(Ga_Unpl.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
Emolow = ft_freqstatistics(cfg,Ga_Unpl, Ga_Neu);
%%
ns=size(Ga_KMiss.powspctrm,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];
HitKMisslow = ft_freqstatistics(cfg,Ga_Hit, Ga_KMiss);
%%
% cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
% save('statTFAmygdalaLF.mat','Int_hitKmisslow','Emolow','HitKMisslow')

cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load statTFAmygdalaLF.mat
%% Supplementary Fig. 2
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg.zlim = [-4 4];
M = Int_hitKmisslow;
M.powspctrm = Int_hitKmisslow.stat.*Int_hitKmisslow.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('int en hit Kmiss'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16);

%%
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative1 = Int_hitKmisslow
logRelative1.mask = Int_hitKmisslow.mask;
logRelative1.powspctrm = Int_hitKmisslow.stat% You shoud add a field mask with the mask from the statistics
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);

%%
clear mat
maskt= sum(squeeze(Int_hitKmisslow.mask),1);
tvec = Int_hitKmisslow.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(Int_hitKmisslow.mask),2);
fvec = Int_hitKmisslow.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t=toi
f=foi

pt1 = nearest(Ga_Unpl.time,t(1));
pt2 = nearest(Ga_Unpl.time,t(2));
pf1 = nearest(Ga_Unpl.freq,f(1));
pf2 = nearest(Ga_Unpl.freq,f(2));

clear mat 
mat(:,1) = squeeze(mean(mean(Ga_eHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_eKMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(Ga_nHit.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(Ga_nKMiss.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

addpath /Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox_plot/beeswarm-master
addpath /Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox_plot/stdshade.m

h = figure;set(h,'position', [0 0 250 188])
ns=size(mat,1)
x = [ones(ns,1) ones(ns,1)*2 ones(ns,1)*3 ones(ns,1)*4];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4)];
beeswarm(x(:),y(:),'sort_style','up','dot_size',2,'overlay_style','sd','colormap',[1 0 0; 1 0 1;0 0 1; 0.5 0.5 0.5])
ylabel('mean theta power','FontSize',10)
xticklabels({'eRHit','eKHit&Miss','nRHit','nKHit&Miss'})
