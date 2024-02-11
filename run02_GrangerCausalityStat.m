% This code reproduce the granger causality analysis in the direction
% amygdala to hippocampus and hippocampus to amygdala as in Fig.2
% Fig 2
clear all
close all
clc

%% To reproduce Figure 2
load Grangerdata.mat
load Grangerstat.mat

h = figure;set(h,'position', [0 0 1000 400]) 
cfg = [];
cfg.figure = 'gcf';
cfg.xlim =[0 1.5]
cfg.ylim = [2 34]
cfg.zlim = [-4 4];
subplot(2,2,1)
logRelative1 = emoneuold_H2A
logRelative1.mask = emoneuold_H2A.mask;
logRelative1.powspctrm = emoneuold_H2A.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('emo vs neu H2A'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title('emotional vs neutral h2a');

hold on
subplot(2,2,2)
E = emoneuold_A2H
E.powspctrm = emoneuold_A2H.stat;
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('emo vs neu A2H'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title('emotional vs neutral a2h');

subplot(2,2,3)
logRelative1 = memhitkmiss_H2A
logRelative1.mask = memhitkmiss_H2A.mask;
logRelative1.powspctrm = memhitkmiss_H2A.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('hit vs hitkmiss H2A'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title('hit vs hitkmiss h2a')

hold on
subplot(2,2,4)
logRelative1 = memhitkmiss_A2H
logRelative1.mask = memhitkmiss_A2H.mask;
logRelative1.powspctrm = memhitkmiss_A2H.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);colorbar;  tit=title('hit vs hitkmiss A2H'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title('hit vs hitkmiss ah2')
