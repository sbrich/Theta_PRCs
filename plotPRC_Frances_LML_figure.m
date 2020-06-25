close all
clear
clc

Iapp=20;
repnum=74;
period=zeros(repnum,1);
c=jet(repnum);
PRCs=zeros(repnum,101);

for i=0:repnum-1
    
    str4=sprintf('InhibitoryNetwork_PRC_LML_%d_%d.csv', i, Iapp);
    str5=sprintf('InhibitoryNetwork_Parameters_LML_%d_%d.csv', i, Iapp);
    
    % trackvariables=csvread(str1);
    % spikes=csvread(str2);
    % currents=csvread(str3);
    PRCdata=csvread(str4);
    variables=csvread(str5);
    
    period(i+1)=PRCdata(1,3);
    PRCs(i+1,:)=PRCdata(1:end,2);
    
%     plot(PRCdata(1:end-1, 1)./variables(1,3), PRCdata(1:end-1,2), 'color', c(i+1,:));
%     axis([0 1 -.5 .1])
%     hold on
end

PRC_std=std(PRCs);
PRC_mean=mean(PRCs);
times=(PRCdata(1:end, 1)./variables(1,3))';

curve1=PRC_mean+PRC_std;
curve2=PRC_mean-PRC_std;
times2 = [times, fliplr(times)];
inBetween = [curve1, fliplr(curve2)];

freq=1./(period./1000);
freq_mean=mean(freq);
freq_STD=std(freq);

figure('units','normalized','position',[0 0 1 1])
h=fill(times2, inBetween, 'b', 'LineStyle','none');
set(h,'facealpha',.25)
hold on
plot(times, PRC_mean, 'b-', 'LineWidth', 3);
hold on
plot(PRCdata(1:end, 1)./variables(1,3), zeros(1,length(PRCdata(1:end,1))), 'k-', 'LineWidth', 2);
axis([0 1 -.5 .15])

str1 = sprintf('Mean Frequency=%1.3f Hz', freq_mean);
str2= sprintf('STD=%1.3f Hz', freq_STD);
str={str1, str2};
dim = [.55 .7 .2 .2];
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize', 34);



set(gca, 'FontSize', 34);
xlabel('Normalized Phase', 'FontSize', 30, 'FontWeight','bold')
ylabel('Normalized Phase Difference', 'FontSize', 30, 'FontWeight','bold')
str=sprintf('LML, Iapp=%d pA', Iapp);
title(str, 'FontSize', 26);


% filename=sprintf('PRC_Iapp%1.0f_gsyn%1.0f_state%1d_alteredmodel.fig', variables(1,2), variables(1,1), variables(1,4));
filename=sprintf('PRCs_figure_LML_Iapp%d_gsyn%1.0f.fig', Iapp, variables(1,1));
savefig(filename);
filename2=sprintf('PRCs_figure_LML_Iapp%d_gsyn%1.0f.png', Iapp, variables(1,1));
set(gcf,'PaperPositionMode','auto')
print(filename2, '-dpng', '-r0');



% freq=1./(period./1000);
% for i=1:repnum
%     plot(i,freq(i), '.', 'MarkerSize', 10, 'color', c(i,:))
%     hold on
% end
% axis([0 repnum-1 0 max(freq)+1])
% set(gca, 'FontSize', 30);
% xlabel('Model #', 'FontSize', 24, 'FontWeight','bold')
% ylabel('Firing Frequency (10th Spike), Hz', 'FontSize', 24, 'FontWeight','bold')
% str=sprintf('LML, Iapp=%d pA', Iapp);
% title(str, 'FontSize', 26);
% filename=sprintf('Freqs_LML_Iapp%d_gsyn%1.0f.fig', Iapp, variables(1,1));
% savefig(filename);

