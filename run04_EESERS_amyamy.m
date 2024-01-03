clear all
close all
clc
addpath('/Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/01_Paper/Githubcodes/Function')

%%
restoredefaultpath
addpath '/Volumes/ManuelaCTB/manuela/Desktop/CTB/000_Toolbox/fieldtrip-20210212/'
ft_defaults

oripath = '/Volumes/ManuelaCTB/manuela/Desktop/AversiveMemRetrieval'; 
addpath(fullfile(oripath,'utils'))
addpath(fullfile(oripath,'utils','beeswarm-master'))
addpath('/Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Function')

%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load TF_data_16sj.mat
load('data_16sj.mat')

subjects = [60 13 15 16 160 25 32 10 33 6 600 10 2 36 38 8];
%%
%[4:8,10,11,13,14] A1; [1:3,9,12,15,16] A2
% [1,3,6:13,16] h1;[2,4,5,14:15] h2
Amylist=[2 2 2 1 1 1 1 1 2 1 1 2 1 1 2 2];
Hclist=[1 2 1 2 2 1 1 1 1 1 1 1 1 2 2 1];
freq = [23:47]

%%
toi=[0 1.5]
parfor s=1:16;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_era(s),TF_ehita(s),Amylist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nra(s),TF_nhita(s),Amylist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efa(s),TF_emissa(s),Amylist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfa(s),TF_nmissa(s),Amylist(s),50,freq,toi);
end

parfor s=1:16;
[TFpeakall_erac1{s}, TFpeakall_erhc1{s}, TFpeakall_ehithc1{s}, nb_er{s}] = findamypeaksEESERS(eR_amy{s},TF_era(s),Amylist(s), TF_era(s),TF_ehita(s),Amylist(s),150,freq,toi);
[TFpeakall_nrac1{s}, TFpeakall_nrhc1{s}, TFpeakall_nhithc1{s}, nb_nr{s}] = findamypeaksEESERS(nR_amy{s},TF_nra(s),Amylist(s), TF_nra(s),TF_nhita(s),Amylist(s),150,freq,toi);
[TFpeakall_efac1{s}, TFpeakall_efhc1{s}, TFpeakall_emisshc1{s}, nb_ef{s}] = findamypeaksEESERS(eKF_amy{s},TF_efa(s),Amylist(s), TF_efa(s),TF_emissa(s),Amylist(s),150,freq,toi);
[TFpeakall_nfac1{s}, TFpeakall_nfhc1{s}, TFpeakall_nmisshc1{s}, nb_nf{s}] = findamypeaksEESERS(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfa(s),TF_nmissa(s),Amylist(s),150,freq,toi);
end

%%
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]
%%
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


%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load Gatemp.mat

        subjects = [60 13 15 16 160 25 32 10 33 6 600 10 2 36 38 8];
subjects_renamed = [60 13 15 16 160 25 32 11 33 6 600 10 2000 36 38 8];
method = 'oneel'
if strcmp(method,'allel')
pos=find(ismember(subjects_renamed,subjects_ers))
else
subjects_ers1 = [60 13 15 160 25 32 33 600 2000 36 38 8];
subjects_ers1me = [60 13 15 160 25 32 11 33 600 10 2000 36 38 8]; %11 cannot be included in the int for nr
pos = find(ismember(subjects_renamed,subjects_ers1))
posme = find(ismember(subjects_renamed,subjects_ers1me))
end

%%
GAs.time = []
GAs.time = TF_erh(1).values.time(151:221)

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:221)

GAs.individual = []

%%
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

%%
URUFahers = Ga_erahers;
URUFahers.individual = Ga_erahers.individual - Ga_efahers.individual
NRNFahers = Ga_erahers;
NRNFahers.individual = Ga_nrahers.individual - Ga_nfahers.individual


URUFahc1ers = Ga_erahc1ers;
URUFahc1ers.individual = Ga_erahc1ers.individual - Ga_efahc1ers.individual
NRNFahc1ers = Ga_erahc1ers;
NRNFahc1ers.individual = Ga_nrahc1ers.individual - Ga_nfahc1ers.individual

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
cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

ns=size(URUFahers.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int_ers = ft_timelockstatistics(cfg,URUFahers, NRNFahers)
intc1_ers= ft_timelockstatistics(cfg,URUFahc1ers, NRNFahc1ers)

%% plot Supplementary Fig.3
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load('statERS_Amygdala.mat')
load ('dataERSAmygdala.mat')
ns=size(Ga_erahers.individual,1)

%% 

clus =[]; i=[]; where =[]; 
    if ~isempty(int_ers.posclusters) 
        clus = find([int_ers.posclusters.prob]'<0.025);%before 0.05
    end 
    
for i = 1:length(clus);
            where = find(int_ers.posclusterslabelmat==i);
            stat_lineers(i,:) = [where(1) where(end)];
end
%%
clusc1 =[]; i=[]; wherec1 =[]; 
    if ~isempty(intc1_ers.posclusters) 
        clusc1 = find([intc1_ers.posclusters.prob]'<0.025);%before 0.05
    end 
    
for i = 1:length(clusc1);
            wherec1 = find(intc1_ers.posclusterslabelmat==i);
            stat_linec1(i,:) = [wherec1(1) wherec1(end)];
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

%%
time=[0:21:1500];
time=time(1,[1:71])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahers.individual,1)), semeRers, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahers.individual,1)), semeFers, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahers.individual,1)), semnRers, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahers.individual,1)), semnFers, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,3)
plot(time,int_ers.stat,'k','LineWidth',3)
%line([time(stat_lineers(1,1)), time(stat_lineers(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); %[0.9 0.5 0.6]
hold on
%line([time(stat_lineers(2,1)), time(stat_lineers(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); %[0.9 0.5 0.6],
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by Memory int; inp=100 p=' num2str(int_ers.negclusters(1).prob,'%0.3f')])
hold on
subplot(4,3,4)
shadedErrorBar(time,squeeze(mean(Ga_erahc1ers.individual,1)), semeRc1, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahc1ers.individual,1)), semeFc1, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,5)
shadedErrorBar(time,squeeze(mean(Ga_nrahc1ers.individual,1)), semnRc1, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahc1ers.individual,1)), semnFc1, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,6)
plot(time,intc1_ers.stat,'k','LineWidth',3)
%line([time(stat_linec1(1,1)), time(stat_linec1(1,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
hold on
%line([time(stat_linec1(2,1)), time(stat_linec1(2,2))],[-2, -2],'LineWidth',5, 'color', 'r','LineStyle','-'); 
xlabel('Time (s)')
ylabel ('t-values')
title(['Emotion by Memory int; inp=300 p=' num2str(intc1_ers.negclusters(1).prob,'%0.3f')])

