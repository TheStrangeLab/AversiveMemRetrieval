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
%%
subjects = [60 13 15 16 160 25 32 10 33 6 600 10 2 36 38 8];
%%
%[4:8,10,11,13,14] A1; [1:3,9,12,15,16] A2
% [1,3,6:13,16] h1;[2,4,5,14:15] h2
Amylist=[2 2 2 1 1 1 1 1 2 1 1 2 1 1 2 2];
Hclist=[1 2 1 2 2 1 1 1 1 1 1 1 1 2 2 1];
freq = [23:47]
toi=[0 1.5]
%%
parfor s=1:16;
[TFpeakall_era{s}, TFpeakall_erh{s}, TFpeakall_ehith{s}, nb_er{s}, pktime_er{s}, pktrltime_er{s}] = findamypeaksEESERS_peaktime(eR_amy{s},TF_era(s),Amylist(s), TF_erh(s),TF_ehith(s),Hclist(s),50,freq,toi);
[TFpeakall_nra{s}, TFpeakall_nrh{s}, TFpeakall_nhith{s}, nb_nr{s}, pktime_nr{s}, pktrltime_nr{s}] = findamypeaksEESERS_peaktime(nR_amy{s},TF_nra(s),Amylist(s), TF_nrh(s),TF_nhith(s),Hclist(s),50,freq,toi);
[TFpeakall_efa{s}, TFpeakall_efh{s}, TFpeakall_emissh{s}, nb_ef{s}, pktime_ef{s}, pktrltime_ef{s}] = findamypeaksEESERS_peaktime(eKF_amy{s},TF_efa(s),Amylist(s), TF_efh(s),TF_emissh(s),Hclist(s),50,freq,toi);
[TFpeakall_nfa{s}, TFpeakall_nfh{s}, TFpeakall_nmissh{s}, nb_nf{s}, pktime_nf{s}, pktrltime_nf{s}] = findamypeaksEESERS_peaktime(nKF_amy{s},TF_nfa(s),Amylist(s), TF_nfh(s),TF_nmissh(s),Hclist(s),50,freq,toi);
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
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.9;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_erh{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 hz;%
%%
for subj=1:16;
[alleraht{subj}] = rsa_peaks(TFpeakall_era{subj},TFpeakall_erh{subj},tsamps,fsh);
[allnraht{subj}] = rsa_peaks(TFpeakall_nra{subj},TFpeakall_nrh{subj},tsamps,fsh);
[allefaht{subj}] = rsa_peaks(TFpeakall_efa{subj},TFpeakall_efh{subj},tsamps,fsh);
[allnfaht{subj}] = rsa_peaks(TFpeakall_nfa{subj},TFpeakall_nfh{subj},tsamps,fsh);
end

%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load data_aheestime.mat
%% plot one subject example as in Fig. 5
figure 
x=[0 1.5]

for s=7
clim=[-0.5 0.5]
data= squeeze(alleraht{s})
y=1:size(data,1)

imagesc(x,y,data,clim)
colorbar
xlabel('time (s)')
ylabel('Amygdala peaks')
axis xy
end
%% plot line 13 to show the max value
figure
x=[0 1.5]

for s=7
clim=[-0.5 0.5]
data= squeeze(alleraht{s})
imagesc(x,1,data(13,:),clim)
colorbar
title(sprintf('%s - Subject: %d', num2str(s)));
xlabel('time (s)')
ylabel('peaks number')
axis xy
end
%% calculate the delay 
testval=max(data(13,:))
[i,j]=find(data(13,:)==testval)
pktimeamy_p7pk13=pktrltime_er{7}{3}(1)
delay_p7pk13=35*1.5/150
hctime_p7pk13= pktimeamy_p7pk13+delay_p7pk13

%%
figure
subplot(1,2,1)
clims=[-2 2]
x = [-0.05 0.05];
y = [90 150];
datapattern=TF_erh(11).values.powspctrm(3,:,23:47,228:238)%233
datapattern=squeeze(datapattern)
imagesc(x,y,datapattern,clims)
set(gca,'YDir','normal')
xlabel ('time (s)')
ylabel ('frequency (Hz)')

subplot(1,2,2)
clims=[-1.5 1.5]
x = [-0.05 0.05];
y = [90 150];
datapattern_rec=TF_ehith(11).values.powspctrm(3,:,23:47,226:236)%233
datapattern_rec=squeeze(datapattern_rec)
imagesc(x,y,datapattern_rec,clims)
set(gca,'YDir','normal')
xlabel ('time (s)')
ylabel ('frequency (Hz)')
% colorbar

%%
win = 5;
freq=[23:47];
t=[0 1.5];
pt1 = nearest(TF_erh(1).values.time,t(1));
pt2 = nearest(TF_erh(1).values.time,t(2));
        
for subj=1:16
clear dataeRsj maxval
dataeRsj= squeeze(alleraht{subj});% this is ees amy hc
ch=Hclist(subj);% channel to use most lateral

for trl=1:size(nb_er{subj},2);
for nbpk=1:(nb_er{subj}{trl});
for pk=1:size(dataeRsj,1);
    clear delay hctime hcdatapoint
    val=max(dataeRsj(pk,:));
    [i,j]=find(dataeRsj(pk,:)==val);
    maxval(pk)=j;
    allval(pk)= val
    amypktime=pktrltime_er{subj}{trl}(nbpk);
    delay=j*1.5/151;
    hctime= amypktime+delay;
    hcdatapoint= nearest(TF_erh(subj).values.time,hctime);
    eRencpattern{trl}(nbpk,:,:,:) = TF_erh(subj).values.powspctrm(1,ch,freq,hcdatapoint-win:hcdatapoint+win);
    eRrecpattern{trl}(nbpk,:,:,:) = TF_ehith(subj).values.powspctrm(1,ch,freq,pt1:pt2);

end
end
end
pkeR{subj}=dataeRsj;
maxvalpkeR{subj}=maxval;
maxalleR{subj}=allval;
TFpeak_erhh{subj}= vertcat(eRencpattern{:})
TFpeak_ehithh{subj}= vertcat(eRrecpattern{:})
end
%%

for subj=1:16
clear dataeFsj maxval allval
dataeFsj= squeeze(allefaht{subj});
ch=Hclist(subj);

for trl=1:size(nb_ef{subj},2);
for nbpk=1:(nb_ef{subj}{trl});
for pk=1:size(dataeFsj,1);
    clear delay hctime hcdatapoint

    val=max(dataeFsj(pk,:));
    [i,j]=find(dataeFsj(pk,:)==val);
    maxval(pk)=j;
    allval(pk)= val

    amypktime=pktrltime_ef{subj}{trl}(nbpk);
    delay=j*1.5/151;
    hctime= amypktime+delay;
    hcdatapoint= nearest(TF_erh(subj).values.time,hctime);
    eFencpattern{trl}(nbpk,:,:,:) = TF_efh(subj).values.powspctrm(1,ch,freq,hcdatapoint-win:hcdatapoint+win);
    eFrecpattern{trl}(nbpk,:,:,:) = TF_emissh(subj).values.powspctrm(1,ch,freq,pt1:pt2);

end
end
end
pkeF{subj}=dataeFsj;
maxvalpkeF{subj}=maxval;
maxalleF{subj}=allval;

TFpeak_efhh{subj}= vertcat(eFencpattern{:})
TFpeak_emisshh{subj}= vertcat(eFrecpattern{:})
end

%%
win = 5;
freq=[23:47];
t=[0 1.5];
pt1 = nearest(TF_erh(1).values.time,t(1));
pt2 = nearest(TF_erh(1).values.time,t(2));
        
for subj=1:16
clear datanRsj maxval allval
datanRsj= squeeze(allnraht{subj});% this is ees amy hc
ch=Hclist(subj);% channel to use most lateral


for trl=1:size(nb_nr{subj},2);
for nbpk=1:(nb_nr{subj}{trl});
for pk=1:size(datanRsj,1);
    clear delay hctime hcdatapoint
    val=max(datanRsj(pk,:));
    [i,j]=find(datanRsj(pk,:)==val);
    maxval(pk)=j;
    allval(pk)= val

    amypktime=pktrltime_nr{subj}{trl}(nbpk);
    delay=j*1.5/151;
    hctime= amypktime+delay;
    hcdatapoint= nearest(TF_nrh(subj).values.time,hctime);
    nRencpattern{trl}(nbpk,:,:,:) = TF_nrh(subj).values.powspctrm(1,ch,freq,hcdatapoint-win:hcdatapoint+win);
    nRrecpattern{trl}(nbpk,:,:,:) = TF_nhith(subj).values.powspctrm(1,ch,freq,pt1:pt2);

end
end
end
pknR{subj}=datanRsj;
maxvalpknR{subj}=maxval;
maxallnR{subj}=allval;

TFpeak_nrhh{subj}= vertcat(nRencpattern{:})
TFpeak_nhithh{subj}= vertcat(nRrecpattern{:})
end
%%

for subj=1:16
clear datanFsj maxval allval
datanFsj= squeeze(allnfaht{subj});
ch=Hclist(subj);

for trl=1:size(nb_nf{subj},2);
for nbpk=1:(nb_nf{subj}{trl});
for pk=1:size(datanFsj,1);
    clear delay hctime hcdatapoint

    val=max(datanFsj(pk,:));
    [i,j]=find(datanFsj(pk,:)==val);
    maxval(pk)=j;
    allval(pk)= val

    amypktime=pktrltime_nf{subj}{trl}(nbpk);
    delay=j*1.5/151;
    hctime= amypktime+delay;
    hcdatapoint= nearest(TF_nrh(subj).values.time,hctime);
    nFencpattern{trl}(nbpk,:,:,:) = TF_nfh(subj).values.powspctrm(1,ch,freq,hcdatapoint-win:hcdatapoint+win);
    nFrecpattern{trl}(nbpk,:,:,:) = TF_nmissh(subj).values.powspctrm(1,ch,freq,pt1:pt2);

end
end
end
pknF{subj}=datanFsj;
maxvalpknF{subj}=maxval;
maxallnF{subj}=allval;

TFpeak_nfhh{subj}= vertcat(nFencpattern{:})
TFpeak_nmisshh{subj}= vertcat(nFrecpattern{:})
end

%% calculate rsa
nsmp    = size(TFpeakall_era{1},3)-1;
overlap = 0.8;
nshift  = round((1-overlap)*nsmp);
endsample = size(TFpeakall_ehith{1},3);

begs = (1:nshift:(endsample-nsmp))';
ends = begs+nsmp;
tsamps= [begs ends]

fsh=[1:25]%this is 90-150 hz;%
%%
parfor subj=1:16;
[allerahh{subj}] = rsa_peaks(TFpeak_erhh{subj},TFpeak_ehithh{subj},tsamps,fsh);
[allnrahh{subj}] = rsa_peaks(TFpeak_nrhh{subj},TFpeak_nhithh{subj},tsamps,fsh);
[allefahh{subj}] = rsa_peaks(TFpeak_efhh{subj},TFpeak_emisshh{subj},tsamps,fsh);
[allnfahh{subj}] = rsa_peaks(TFpeak_nfhh{subj},TFpeak_nmisshh{subj},tsamps,fsh);
end

%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load data_ahh_ers.mat
%%
for subj=1:16
[i,j,v]=find(maxalleR{subj} < 0)
allerahh_pos{subj}=allerahh{subj};
allerahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:16
[i,j,v]=find(maxalleF{subj} < 0)
allefahh_pos{subj}=allefahh{subj};
allefahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:16
[i,j,v]=find(maxallnR{subj} < 0)
allnrahh_pos{subj}=allnrahh{subj};
allnrahh_pos{subj}(v,:)=[];
end

clear v
for subj=1:16
[i,j,v]=find(maxallnF{subj} < 0)
allnfahh_pos{subj}=allnfahh{subj};
allnfahh_pos{subj}(v,:)=[];
end
%%

cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load Gatemp.mat

%%
GAs.time = []
GAs.time = TF_erh(1).values.time(151:221)% 80% overalp; 100 ms

GAs.stderr = []
GAs.stderr = TF_erh(1).values.time(151:221)

GAs.individual = []
%%
Ga_erahh = GAs
Ga_nrahh = GAs
Ga_efahh = GAs
Ga_nfahh = GAs

v=0
for s=pos;
v=v+1
Ga_erahh.individual(v,1,:)= mean(allerahh{s},1);
Ga_nrahh.individual(v,1,:)= mean(allnrahh{s},1);
Ga_efahh.individual(v,1,:)= mean(allefahh{s},1);
Ga_nfahh.individual(v,1,:)= mean(allnfahh{s},1);
end


%%
URUFahh = Ga_erahh;
URUFahh.individual = Ga_erahh.individual - Ga_efahh.individual
NRNFahh = Ga_erahh;
NRNFahh.individual = Ga_nrahh.individual - Ga_nfahh.individual

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

ns=size(Ga_erahh.individual,1)
cfg.design = [ones(1,ns) ones(1,ns).*2;[1:ns] [1:ns]];

int = ft_timelockstatistics(cfg,URUFahh, NRNFahh)

%%
cd /Volumes/ManuelaCTB/manuela/Desktop/CTB/00000_IAPS_ERS/001Submission/22December2023/Githubcodes/Data
load stat_ahhERS.mat

ns=size(Ga_erahh.individual,1);

clus =[]; i=[]; where =[]; 
    if ~isempty(int.posclusters) 
        clus = find([int.posclusters.prob]'<0.025);
    end 
    
for i = 1:length(clus);
            where = find(int.posclusterslabelmat==i);
            stat_line(i,:) = [where(1) where(end)];
end

for i = 1:length(clus);
            where = find(int.negclusterslabelmat==i);
            stat_lineneg(i,:) = [where(1) where(end)];
end
%%
semeR = std(Ga_erahh.individual/sqrt(ns))
semeF = std(Ga_efahh.individual/sqrt(ns))
semnR = std(Ga_nrahh.individual/sqrt(ns))
semnF = std(Ga_nfahh.individual/sqrt(ns))

%%
time=[0:21:1500];
time=time(1,[1:71])/1000

r=figure;set(r,'position', [0 0 1200 1500])
subplot(4,3,1)
shadedErrorBar(time,squeeze(mean(Ga_erahh.individual,1)), semeR, 'r', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_efahh.individual,1)), semeF, 'm', 1); hold on
%legend('er','ef')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,2)
shadedErrorBar(time,squeeze(mean(Ga_nrahh.individual,1)), semnR, 'b', 1); hold on
shadedErrorBar(time,squeeze(mean(Ga_nfahh.individual,1)), semnF, 'k', 1); hold on
%legend('nr','nf')
ylabel('Corr Coef')
xlabel('Time (s)')
hold on
subplot(4,3,3)
plot(time,int.stat,'k','LineWidth',3)
line([time(stat_line(1,1)), time(stat_line(1,2))],[-5, -5],'LineWidth',5, 'color', 'r','LineStyle','-'); 
ylabel ('t-values')
title(['Emotion by memory int;  p=' num2str(int.posclusters(1).prob,'%0.3f')])
xlabel('Time (s)')
hold on
