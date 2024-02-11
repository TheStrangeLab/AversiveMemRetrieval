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

%% find amygdala peaks during encoding and prepare the data for the EES and ERS analysis
parfor s=1:16;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
end

%% control analysis 300 ipd Supplementary Fig. 5a & 6a
parfor s=1:16;
[TFpeakall_erac1{s}, TFpeakall_erhc1{s}, TFpeakall_ehithc1{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),150,freq,toi);
[TFpeakall_nrac1{s}, TFpeakall_nrhc1{s}, TFpeakall_nhithc1{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),150,freq,toi);
[TFpeakall_efac1{s}, TFpeakall_efhc1{s}, TFpeakall_emisshc1{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),150,freq,toi);
[TFpeakall_nfac1{s}, TFpeakall_nfhc1{s}, TFpeakall_nmisshc1{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),150,freq,toi);
end
%% control analysis random peak selection Supplementary Fig. 5b & 6b
parfor s=1:16;
[TFpeakall_erarnd{s}, TFpeakall_erhrnd{s}, TFpeakall_ehithrnd{s}] = findamypeaksEESERS_rnd(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nrarnd{s}, TFpeakall_nrhrnd{s}, TFpeakall_nhithrnd{s}] = findamypeaksEESERS_rnd(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efarnd{s}, TFpeakall_efhrnd{s}, TFpeakall_emisshrnd{s}] = findamypeaksEESERS_rnd(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfarnd{s}, TFpeakall_nfhrnd{s}, TFpeakall_nmisshrnd{s}] = findamypeaksEESERS_rnd(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
end
%% compute rsa for Encoding-Encoding similarity between amygdala and hippocampus
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 Hz;%

parfor subj=1:16;
[alleraht{subj}] = rsa_peaks(TFpeakall_era{subj},TFpeakall_erh{subj},tsamps,fsh);
[allnraht{subj}] = rsa_peaks(TFpeakall_nra{subj},TFpeakall_nrh{subj},tsamps,fsh);
[allefaht{subj}] = rsa_peaks(TFpeakall_efa{subj},TFpeakall_efh{subj},tsamps,fsh);
[allnfaht{subj}] = rsa_peaks(TFpeakall_nfa{subj},TFpeakall_nfh{subj},tsamps,fsh);
end

parfor subj=1:16;
[alleraht_c1{subj}] = rsa_peaks(TFpeakall_erac1{subj},TFpeakall_erhc1{subj},tsamps,fsh);
[allnraht_c1{subj}] = rsa_peaks(TFpeakall_nrac1{subj},TFpeakall_nrhc1{subj},tsamps,fsh);
[allefaht_c1{subj}] = rsa_peaks(TFpeakall_efac1{subj},TFpeakall_efhc1{subj},tsamps,fsh);
[allnfaht_c1{subj}] = rsa_peaks(TFpeakall_nfac1{subj},TFpeakall_nfhc1{subj},tsamps,fsh);
end

for subj=1:16;
[alleraht_rnd{subj}] = rsa_peaks(TFpeakall_erarnd{subj},TFpeakall_erhrnd{subj},tsamps,fsh);
[allnraht_rnd{subj}] = rsa_peaks(TFpeakall_nrarnd{subj},TFpeakall_nrhrnd{subj},tsamps,fsh);
[allefaht_rnd{subj}] = rsa_peaks(TFpeakall_efarnd{subj},TFpeakall_efhrnd{subj},tsamps,fsh);
[allnfaht_rnd{subj}] = rsa_peaks(TFpeakall_nfarnd{subj},TFpeakall_nfhrnd{subj},tsamps,fsh);
end

%% compute rsa for Encoding-Retrieval similarity between amygdala and hippocampus
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]

parfor subj=1:16;
[alleraht_ers{subj}] = rsa_peaks(TFpeakall_era{subj},TFpeakall_ehith{subj},tsamps,fsh);
[allnraht_ers{subj}] = rsa_peaks(TFpeakall_nra{subj},TFpeakall_nhith{subj},tsamps,fsh);
[allefaht_ers{subj}] = rsa_peaks(TFpeakall_efa{subj},TFpeakall_emissh{subj},tsamps,fsh);
[allnfaht_ers{subj}] = rsa_peaks(TFpeakall_nfa{subj},TFpeakall_nmissh{subj},tsamps,fsh);
end

parfor subj=1:16;
[alleraht_c1_ers{subj}] = rsa_peaks(TFpeakall_erac1{subj},TFpeakall_ehithc1{subj},tsamps,fsh);
[allnraht_c1_ers{subj}] = rsa_peaks(TFpeakall_nrac1{subj},TFpeakall_nhithc1{subj},tsamps,fsh);
[allefaht_c1_ers{subj}] = rsa_peaks(TFpeakall_efac1{subj},TFpeakall_emisshc1{subj},tsamps,fsh);
[allnfaht_c1_ers{subj}] = rsa_peaks(TFpeakall_nfac1{subj},TFpeakall_nmisshc1{subj},tsamps,fsh);
end

for subj=1:16;
[alleraht_rnd_ers{subj}] = rsa_peaks(TFpeakall_erarnd{subj},TFpeakall_ehithrnd{subj},tsamps,fsh);
[allnraht_rnd_ers{subj}] = rsa_peaks(TFpeakall_nrarnd{subj},TFpeakall_nhithrnd{subj},tsamps,fsh);
[allefaht_rnd_ers{subj}] = rsa_peaks(TFpeakall_efarnd{subj},TFpeakall_emisshrnd{subj},tsamps,fsh);
[allnfaht_rnd_ers{subj}] = rsa_peaks(TFpeakall_nfarnd{subj},TFpeakall_nmisshrnd{subj},tsamps,fsh);
end

%%
load Gatemp.mat

subjects_renamed = [60 13 15 16 160 25 32 11 33 6 600 10 2000 36 38 8];
method = 'oneel'

subjects_ers = [60 13 15 160 25 32 33 600 2000 36 38 8];
subjects_ersme = [60 13 15 160 25 32 11 33 600 10 2000 36 38 8]; %11 cannot be included in the int for nr
pos = find(ismember(subjects_renamed,subjects_ers));
posme = find(ismember(subjects_renamed,subjects_ersme));

GAs.time = []
GAs.time = TF_erh(1).values.time(151:221);

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:221);

GAs.individual = []
Ga_erah = GAs
Ga_nrah = GAs
Ga_efah = GAs
Ga_nfah = GAs

v=0
for s=pos;
v=v+1
Ga_erah.individual(v,1,:)= mean(alleraht{s},1);
Ga_nrah.individual(v,1,:)= mean(allnraht{s},1);
Ga_efah.individual(v,1,:)= mean(allefaht{s},1);
Ga_nfah.individual(v,1,:)= mean(allnfaht{s},1);
end
%%
Ga_erahc1 = GAs
Ga_nrahc1 = GAs
Ga_efahc1 = GAs
Ga_nfahc1 = GAs

v=0
for s=pos;
v=v+1
Ga_erahc1.individual(v,1,:)= mean(alleraht_c1{s},1);
Ga_nrahc1.individual(v,1,:)= mean(allnraht_c1{s},1);
Ga_efahc1.individual(v,1,:)= mean(allefaht_c1{s},1);
Ga_nfahc1.individual(v,1,:)= mean(allnfaht_c1{s},1);
end

%%
Ga_erah_rnd = GAs
Ga_nrah_rnd = GAs
Ga_efah_rnd = GAs
Ga_nfah_rnd = GAs
v=0
for s=pos;
v=v+1
Ga_erah_rnd.individual(v,1,:)= mean(alleraht_rnd{s},1);
Ga_nrah_rnd.individual(v,1,:)= mean(allnraht_rnd{s},1);
Ga_efah_rnd.individual(v,1,:)= mean(allefaht_rnd{s},1);
Ga_nfah_rnd.individual(v,1,:)= mean(allnfaht_rnd{s},1);
end

%%
URUFah = Ga_erah;
URUFah.individual = Ga_erah.individual - Ga_efah.individual
NRNFah = Ga_erah;
NRNFah.individual = Ga_nrah.individual - Ga_nfah.individual
%%
URUFahc1 = Ga_erahc1;
URUFahc1.individual = Ga_erahc1.individual - Ga_efahc1.individual
NRNFahc1 = Ga_erahc1;
NRNFahc1.individual = Ga_nrahc1.individual - Ga_nfahc1.individual
%%
URUFah_rnd = Ga_erah_rnd;
URUFah_rnd.individual = Ga_erah_rnd.individual - Ga_efah_rnd.individual
NRNFah_rnd = Ga_erah_rnd;
NRNFah_rnd.individual = Ga_nrah_rnd.individual - Ga_nfah_rnd.individual

%% cluter stat 

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05
cfg.numrandomization = 10000;
%cfg.latency     = [0 1.5]
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFah.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int = ft_timelockstatistics(cfg,URUFah, NRNFah)
intc1= ft_timelockstatistics(cfg,URUFahc1, NRNFahc1)
intrnd= ft_timelockstatistics(cfg,URUFah_rnd, NRNFah_rnd)

%% ERS results
GAs.time = []
GAs.time = TF_erh(1).values.time(151:221)

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:221)

GAs.individual = []
Ga_erahers = GAs
Ga_nrahers = GAs
Ga_efahers = GAs
Ga_nfahers = GAs

v=0
for s=pos;
v=v+1
Ga_erahers.individual(v,1,:)= mean(alleraht_ers{s},1);
Ga_nrahers.individual(v,1,:)= mean(allnraht_ers{s},1);
Ga_efahers.individual(v,1,:)= mean(allefaht_ers{s},1);
Ga_nfahers.individual(v,1,:)= mean(allnfaht_ers{s},1);
end


Ga_erahc1ers = GAs
Ga_nrahc1ers = GAs
Ga_efahc1ers = GAs
Ga_nfahc1ers = GAs

v=0
for s=pos;
v=v+1
Ga_erahc1ers.individual(v,1,:)= mean(alleraht_c1_ers{s},1);
Ga_nrahc1ers.individual(v,1,:)= mean(allnraht_c1_ers{s},1);
Ga_efahc1ers.individual(v,1,:)= mean(allefaht_c1_ers{s},1);
Ga_nfahc1ers.individual(v,1,:)= mean(allnfaht_c1_ers{s},1);
end


Ga_erah_rnders = GAs
Ga_nrah_rnders = GAs
Ga_efah_rnders = GAs
Ga_nfah_rnders = GAs
v=0
for s=pos;
v=v+1
Ga_erah_rnders.individual(v,1,:)= mean(alleraht_rnd_ers{s},1);
Ga_nrah_rnders.individual(v,1,:)= mean(allnraht_rnd_ers{s},1);
Ga_efah_rnders.individual(v,1,:)= mean(allefaht_rnd_ers{s},1);
Ga_nfah_rnders.individual(v,1,:)= mean(allnfaht_rnd_ers{s},1);
end

%%
URUFahers = Ga_erahers;
URUFahers.individual = Ga_erahers.individual - Ga_efahers.individual
NRNFahers = Ga_erahers;
NRNFahers.individual = Ga_nrahers.individual - Ga_nfahers.individual


URUFahc1ers = Ga_erahc1ers;
URUFahc1ers.individual = Ga_erahc1ers.individual - Ga_efahc1ers.individual
NRNFahc1ers = Ga_erahc1ers;
NRNFahc1ers.individual = Ga_nrahc1ers.individual - Ga_nfahc1ers.individual


URUFah_rnders = Ga_erah_rnders;
URUFah_rnders.individual = Ga_erah_rnders.individual - Ga_efah_rnders.individual
NRNFah_rnders = Ga_erah_rnders;
NRNFah_rnders.individual = Ga_nrah_rnders.individual - Ga_nfah_rnders.individual

ns=size(URUFahers.individual,1)

%% cluter stat 

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05
cfg.numrandomization = 10000;
%cfg.latency     = [0 1.5]
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFahers.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int_ers = ft_timelockstatistics(cfg,URUFahers, NRNFahers)
intc1_ers= ft_timelockstatistics(cfg,URUFahc1ers, NRNFahc1ers)
intrndnop_ers= ft_timelockstatistics(cfg,URUFah_rndnopers, NRNFah_rndnopers)
intrnd_ers= ft_timelockstatistics(cfg,URUFah_rnders, NRNFah_rnders)
%% to plot results in Fig.4 and Supplementary Fig. 5 & 6
clear all
close all
clc

restoredefaultpath
addpath '~/fieldtrip-20210212/'
ft_defaults

addpath('~/utils/reinstatement-analysis-master/additional_functions')
addpath('~/utils/functions')
addpath('~/utils/beeswarm-master')
addpath('~/utils')

%%
cd ~/data2github
load  dataEES.mat;
load  statEES.mat;

%% 
ns=size(Ga_erah.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int.posclusters) 
        clus = find([int.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int.posclusterslabelmat==i);
            stat_line(i,:) = [where(1) where(end)];
end

clusc1 =[]; i=[]; wherec1 =[]; 
    if ~isempty(intc1.posclusters) 
        clusc1 = find([intc1.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clusc1);
            wherec1 = find(intc1.posclusterslabelmat==i);
            stat_linec1(i,:) = [wherec1(1) wherec1(end)];
end

clus_rnd =[]; i=[]; where_rnd =[]; 
    if ~isempty(intrnd.posclusters) 
        clus_rnd = find([intrnd.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus_rnd);
            where_rnd = find(intrnd.posclusterslabelmat==i);
            stat_line_rnd(i,:) = [where_rnd(1) where_rnd(end)];
end

%%
semeR = std(Ga_erah.individual/sqrt(ns))
semeF = std(Ga_efah.individual/sqrt(ns))
semnR = std(Ga_nrah.individual/sqrt(ns))
semnF = std(Ga_nfah.individual/sqrt(ns))

semeRc1 = std(Ga_erahc1.individual/sqrt(ns))
semeFc1 = std(Ga_efahc1.individual/sqrt(ns))
semnRc1 = std(Ga_nrahc1.individual/sqrt(ns))
semnFc1 = std(Ga_nfahc1.individual/sqrt(ns))

semeRrnd = std(Ga_erah_rnd.individual/sqrt(ns))
semeFrnd = std(Ga_efah_rnd.individual/sqrt(ns))
semnRrnd = std(Ga_nrah_rnd.individual/sqrt(ns))
semnFrnd = std(Ga_nfah_rnd.individual/sqrt(ns))

%% Supplementary Fig. 6

time=[0:21:1500];
time=time(1,[1:71])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahc1.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,3)
plot(time,intc1.stat,'k','LineWidth',3)
line([time(stat_linec1(1,1)), time(stat_linec1(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time from amygdala peaks (s)')
ylabel ('t-values')
title('emotion by memory effect');
hold on
subplot(4,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah_rnd.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_rnd.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah_rnd.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_rnd.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time from amygdala peaks (s)')
hold on
subplot(4,3,6)
plot(time,intrnd.stat,'k','LineWidth',3)
line([time(stat_line_rnd(1,1)), time(stat_line_rnd(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time from amygdala peaks (s)')
ylabel ('t-values')
title('emotion by memory effect');
hold on
%%
cd ~/data2github
load dataERS.mat
load statERS.mat
%% 

clus =[]; i=[]; where =[]; 
    if ~isempty(int_ers.posclusters) 
        clus = find([int_ers.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int_ers.posclusterslabelmat==i);
            stat_lineers(i,:) = [where(1) where(end)];
end

clusc1 =[]; i=[]; wherec1 =[]; 
    if ~isempty(intc1_ers.posclusters) 
        clusc1 = find([intc1_ers.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clusc1);
            wherec1 = find(intc1_ers.posclusterslabelmat==i);
            stat_linec1(i,:) = [wherec1(1) wherec1(end)];
end

clus_rnd =[]; i=[]; where_rnd =[]; 
    if ~isempty(intrnd_ers.posclusters) 
        clus_rnd = find([intrnd_ers.posclusters.prob]'<0.08);% result is marked with a different color
    end 
    
for i = 1:length(clus_rnd);
            where_rnd = find(intrnd_ers.posclusterslabelmat==i);
            stat_line_rnd(i,:) = [where_rnd(1) where_rnd(end)];
end
%%
semeRers = std(Ga_erahers.individual/sqrt(ns))
semeFers = std(Ga_efahers.individual/sqrt(ns))
semnRers = std(Ga_nrahers.individual/sqrt(ns))
semnFers = std(Ga_nfahers.individual/sqrt(ns))

semeRc1 = std(Ga_erahc1ers.individual/sqrt(ns))
semeFc1 = std(Ga_efahc1ers.individual/sqrt(ns))
semnRc1 = std(Ga_nrahc1ers.individual/sqrt(ns))
semnFc1 = std(Ga_nfahc1ers.individual/sqrt(ns))

semeRrnd = std(Ga_erah_rnders.individual/sqrt(ns))
semeFrnd = std(Ga_efah_rnders.individual/sqrt(ns))
semnRrnd = std(Ga_nrah_rnders.individual/sqrt(ns))
semnFrnd = std(Ga_nfah_rnders.individual/sqrt(ns))

%% plot Supplementary Fig.5
time=[0:21:1500];
time=time(1,[1:71])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahc1ers.individual,1)), semeRc1, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1ers.individual,1)), semeFc1, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1ers.individual,1)), semnRc1, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1ers.individual,1)), semnFc1, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,3)
plot(time,intc1_ers.stat,'k','LineWidth',3)
line([time(stat_linec1(1,1)), time(stat_linec1(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
hold on
line([time(stat_linec1(2,1)), time(stat_linec1(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by memory int; inp=300 p=' num2str(intc1_ers.posclusters(1).prob,'%0.3f')])
hold on
subplot(4,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah_rnders.individual,1)), semeRrnd, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah_rnders.individual,1)), semeFrnd, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah_rnders.individual,1)), semnRrnd, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah_rnders.individual,1)), semnFrnd, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,6)
plot(time,intrnd_ers.stat,'k','LineWidth',3)
%line([time(stat_line_rnd(1,1)), time(stat_line_rnd(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
line([time(stat_line_rnd(1,1)), time(stat_line_rnd(1,2))],[-2, -2],'LineWidth',5, 'color', [0.9 0.5 0.6],'LineStyle','-'); 
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by memory int; rnd peak p=' num2str(intrnd_ers.posclusters(1).prob, '%0.3f')])


%% plot Fig. 4

time=[0:21:1500];
time=time(1,[1:71])/1000
r = figure;set(r,'position', [0 0 1200 400])
subplot(2,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahers.individual,1)), semeRers, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahers.individual,1)), semeFers, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahers.individual,1)), semnRers, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahers.individual,1)), semnFers, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,3)
plot(time,int_ers.stat,'k','LineWidth',3)
line([time(stat_lineers(1,1)), time(stat_lineers(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); %[0.9 0.5 0.6]
hold on
line([time(stat_lineers(2,1)), time(stat_lineers(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('hc time from amy peaks (s)')
ylabel ('t-values')
title('emotion by memory: ERS amy hc')
hold on
subplot(2,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erah.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efah.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrah.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfah.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('hc time from amy peaks (s)')
hold on
subplot(2,3,6)
plot(time,int.stat,'k','LineWidth',3)
line([time(stat_line(1,1)), time(stat_line(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('hc time from amy peaks (s)')
ylabel ('t-values')
title('emotion by memory: EES amy hc')
