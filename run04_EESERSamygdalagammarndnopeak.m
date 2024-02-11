clear all
close all
clc

%%
restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')
addpath('~/utils')

%%
cd ~/data2github
load TF_data_16sj.mat
load('data_16sj.mat')

subjects = [60 13 15 16 160 25 32 10 33 6 600 10 2 36 38 8];
%%
Amylist=[2 2 2 1 1 1 1 1 2 1 1 2 1 1 2 2];
Hclist=[1 2 1 2 2 1 1 1 1 1 1 1 1 2 2 1];
freq = [23:47]
toi=[0 1.5]
%% find amygdala peaks during encoding and prepare the data
parfor s=1:16;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
end
%%
load Gatemp.mat

subjects_renamed = [60 13 15 16 160 25 32 11 33 6 600 10 2000 36 38 8];
method = 'oneel'
subjects_ers = [60 13 15 160 25 32 33 600 2000 36 38 8];
subjects_ersme = [60 13 15 160 25 32 11 33 600 10 2000 36 38 8]; 
pos = find(ismember(subjects_renamed,subjects_ers))
posme = find(ismember(subjects_renamed,subjects_ersme))

%%
GAs.time = []
GAs.time = TF_erh(1).values.time(151:221)
GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:221)
GAs.individual = []
%%
Ga_erah_rndnop = GAs
Ga_nrah_rndnop = GAs
Ga_efah_rndnop = GAs
Ga_nfah_rndnop = GAs
%%
eR_amysel = eR_amy(1,pos);
nR_amysel = nR_amy(1,pos);
eKF_amysel = eKF_amy(1,pos);
nKF_amysel = nKF_amy(1,pos);

TF_erasel = TF_era(1,pos);
TF_nrasel = TF_nra(1,pos);
TF_efasel = TF_efa(1,pos);
TF_nfasel = TF_nfa(1,pos);

TF_erhsel = TF_erh(1,pos);
TF_nrhsel = TF_nrh(1,pos);
TF_efhsel = TF_efh(1,pos);
TF_nfhsel = TF_nfh(1,pos);

TF_ehithsel = TF_ehith(1,pos);
TF_nhithsel = TF_nhith(1,pos);
TF_emisshsel = TF_emissh(1,pos);
TF_nmisshsel = TF_nmissh(1,pos);

Amylistsel= Amylist(1,pos);
Hclistsel= Hclist(1,pos);
%% compute the rsa using random no peak selection
tic
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

nbsmp= size(tsamps,1)
Ga_erah_rndnop = zeros(12, 1, nbsmp);
Ga_nrah_rndnop = zeros(12, 1, nbsmp);
Ga_efah_rndnop = zeros(12, 1, nbsmp);
Ga_nfah_rndnop = zeros(12, 1, nbsmp);


for rnd=1:1000; % repeat the analysis 1000 times and get a p-value
parfor s=1:12;
[TFpeakall_erarndnop{s}, TFpeakall_erhrndnop{s}, TFpeakall_ehithrndnop{s}] = findamypeaksEESERS_rndnop(eR_amysel{s},TF_erasel(s),Amylistsel(s), TF_erhsel(s),TF_ehithsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_nrarndnop{s}, TFpeakall_nrhrndnop{s}, TFpeakall_nhithrndnop{s}] = findamypeaksEESERS_rndnop(nR_amysel{s},TF_nrasel(s),Amylistsel(s), TF_nrhsel(s),TF_nhithsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_efarndnop{s}, TFpeakall_efhrndnop{s}, TFpeakall_emisshrndnop{s}] = findamypeaksEESERS_rndnop(eKF_amysel{s},TF_efasel(s),Amylistsel(s), TF_efhsel(s),TF_emisshsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_nfarndnop{s}, TFpeakall_nfhrndnop{s}, TFpeakall_nmisshrndnop{s}] = findamypeaksEESERS_rndnop(nKF_amysel{s},TF_nfasel(s),Amylistsel(s), TF_nfhsel(s),TF_nmisshsel(s),Hclistsel(s),50,freq,toi);

fsh=[1:25]%this is 90-150 hz;%

[alleraht_rnd_nop{s}] = rsa_peaks(TFpeakall_erarndnop{s},TFpeakall_erhrndnop{s},tsamps,fsh);
[allnraht_rnd_nop{s}] = rsa_peaks(TFpeakall_nrarndnop{s},TFpeakall_nrhrndnop{s},tsamps,fsh);
[allefaht_rnd_nop{s}] = rsa_peaks(TFpeakall_efarndnop{s},TFpeakall_efhrndnop{s},tsamps,fsh);
[allnfaht_rnd_nop{s}] = rsa_peaks(TFpeakall_nfarndnop{s},TFpeakall_nfhrndnop{s},tsamps,fsh);

Ga_erah_rndnop(s,1,:)= mean(alleraht_rnd_nop{s},1);
Ga_nrah_rndnop(s,1,:)= mean(allnraht_rnd_nop{s},1);
Ga_efah_rndnop(s,1,:)= mean(allefaht_rnd_nop{s},1);
Ga_nfah_rndnop(s,1,:)= mean(allnfaht_rnd_nop{s},1);

end

URUFah_rndnop = GAs;
URUFah_rndnop.individual = Ga_erah_rndnop - Ga_efah_rndnop;
NRNFah_rndnop = GAs;
NRNFah_rndnop.individual = Ga_nrah_rndnop - Ga_nfah_rndnop;

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 1000;
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=12;
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

intrndnop{rnd} = ft_timelockstatistics(cfg,URUFah_rndnop, NRNFah_rndnop);
end
toc

%%
for rnd=1:1000
if  ~isfield(intrndnop{rnd}, 'posclusters') %si no existe posclusters = 0, si existe = 1
nocluster{rnd} = intrndnop{rnd}.prob(1,1);   
elseif isfield(intrndnop{rnd}, 'posclusters')    
positive_cluster_value(rnd) = intrndnop{rnd}.posclusters.prob;
end
end
%%
figure
histogram(positive_cluster_value)
xlabel 'p-values'
ylabel 'Number of observations'
%% Do the same for the ERS
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]

nbsmp= size(tsamps,1)
Ga_erah_rndnopers = zeros(12, 1, nbsmp);
Ga_nrah_rndnopers = zeros(12, 1, nbsmp);
Ga_efah_rndnopers = zeros(12, 1, nbsmp);
Ga_nfah_rndnopers = zeros(12, 1, nbsmp);


for rnd=1:1000;
parfor s=1:12;
[TFpeakall_erarndnopers{s}, TFpeakall_ehithrndnopers{s}] = findamypeaksERS_rndnop(eR_amysel{s},TF_erasel(s),Amylistsel(s), TF_erhsel(s),TF_ehithsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_nrarndnopers{s}, TFpeakall_nhithrndnopers{s}] = findamypeaksERS_rndnop(nR_amysel{s},TF_nrasel(s),Amylistsel(s), TF_nrhsel(s),TF_nhithsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_efarndnopers{s}, TFpeakall_emisshrndnopers{s}] = findamypeaksERS_rndnop(eKF_amysel{s},TF_efasel(s),Amylistsel(s), TF_efhsel(s),TF_emisshsel(s),Hclistsel(s),50,freq,toi);
[TFpeakall_nfarndnopers{s}, TFpeakall_nmisshrndnopers{s}] = findamypeaksERS_rndnop(nKF_amysel{s},TF_nfasel(s),Amylistsel(s), TF_nfhsel(s),TF_nmisshsel(s),Hclistsel(s),50,freq,toi);

fsh=[1:25]%this is 90-150 hz;%

[alleraht_rnd_nop_ers{s}] = rsa_peaks(TFpeakall_erarndnopers{s},TFpeakall_ehithrndnopers{s},tsamps,fsh);
[allnraht_rnd_nop_ers{s}] = rsa_peaks(TFpeakall_nrarndnopers{s},TFpeakall_nhithrndnopers{s},tsamps,fsh);
[allefaht_rnd_nop_ers{s}] = rsa_peaks(TFpeakall_efarndnopers{s},TFpeakall_emisshrndnopers{s},tsamps,fsh);
[allnfaht_rnd_nop_ers{s}] = rsa_peaks(TFpeakall_nfarndnopers{s},TFpeakall_nmisshrndnopers{s},tsamps,fsh);


Ga_erah_rndnopers(s,1,:)= mean(alleraht_rnd_nop_ers{s},1);
Ga_nrah_rndnopers(s,1,:)= mean(allnraht_rnd_nop_ers{s},1);
Ga_efah_rndnopers(s,1,:)= mean(allefaht_rnd_nop_ers{s},1);
Ga_nfah_rndnopers(s,1,:)= mean(allnfaht_rnd_nop_ers{s},1);

end

URUFah_rndnopers = GAs;
URUFah_rndnopers.individual = Ga_erah_rndnopers - Ga_efah_rndnopers;
NRNFah_rndnopers = GAs;
NRNFah_rndnopers.individual = Ga_nrah_rndnopers - Ga_nfah_rndnopers;

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 1000;%'all';
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=12;
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

intrndnop_ers{rnd} = ft_timelockstatistics(cfg,URUFah_rndnopers, NRNFah_rndnopers);
end
toc

%%
for rnd=1:1000
if  ~isfield(intrndnop_ers{rnd}, 'posclusters') %si no existe posclusters = 0, si existe = 1
nocluster_ers{rnd} = intrndnop_ers{rnd}.prob(1,1);   
elseif isfield(intrndnop_ers{rnd}, 'posclusters')    
positive_cluster_value_ers(rnd) = intrndnop_ers{rnd}.posclusters.prob
end
end
%% to plot Supplementary Fig. 4 starts from here
clear all
close all
clc

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')
addpath('~/utils')

%%
cd ~/data2github
load 'statees_rndnopeaks1000.mat'
load 'stat_rndnopeaks_ers1000.mat'

%%
for rnd=1:1000
if  ~isfield(intrndnop{rnd}, 'posclusters') %si no existe posclusters = 0, si existe = 1
nocluster{rnd} = intrndnop{rnd}.prob(1,1);   
elseif isfield(intrndnop{rnd}, 'posclusters')    
positive_cluster_value(rnd) = intrndnop{rnd}.posclusters.prob;
end
end
%%
for rnd=1:1000
if  ~isfield(intrndnop_ers{rnd}, 'posclusters') %si no existe posclusters = 0, si existe = 1
nocluster_ers{rnd} = intrndnop_ers{rnd}.prob(1,1);   
elseif isfield(intrndnop_ers{rnd}, 'posclusters')    
positive_cluster_value_ers(rnd) = intrndnop_ers{rnd}.posclusters.prob;
end
end
%% Supplementary Fig. 4
figure
subplot(1,2,1)
edges = linspace(0, 0.15, 41);
histogram(positive_cluster_value,'BinEdges',edges, 'FaceColor',[0 0 0])
xlim([0, 0.15]);
ylim([0 180]);
xlabel 'p-values'
ylabel 'Number of observations'
hold all
line([0.025,0.025], ylim, 'LineWidth', 2, 'Color', 'r');

hold on
subplot(1,2,2)
histogram(positive_cluster_value_ers,'BinEdges',edges,'FaceColor',[0 0 0])
xlim([0, 0.15]);
ylim([0 250]);
xlabel 'p-values'
ylabel 'Number of observations'
hold all
line([0.025,0.025], ylim, 'LineWidth', 2, 'Color', 'r');

