%% Calculate selection coefficients and z-scores for cytometry competition experiments
addpath('../../../Scriptoma/')
temp = dir('../001_DataOrganization/PL*');

[num txt] = xlsread('../001_RawData/MoBYComplexes_plate_guide.xlsx');
ind = 1:80;
names = txt(2:81,2:6);
for ii = 1:6%1:length(temp)
    fname = strcat('../002_DataOrganized/',temp(ii).name);

    [num txt] = xlsread(fname);
    data(ii).RoF = num(ind,:);
    pos = names(ind,1);
    data(ii).orf = strcat(pos,'-',names(ind,2)); 
    data(ii).genename = names(ind,3);
    data(ii).complex = names(ind,4);
    
    
    
end
% Cytometry data analysis: selection coefficient calculation and curves
close all
generations = (1:6)*log2(160/10);%calcular el numero de generaciones con la dil 5 en 155ul
inn = [1:20;21:40;41:60;61:80];
%figures plate 2
x = generations; sets = [1:3;4:6];
clrs = 'bbbrrr';

counter = 0;
for wii = 1:2 %number of sets (2 sets,3replicas)
    plts = sets(wii,:);
    
    for wi = 1:4
    counter = counter+1;
    figure(counter);clf
    subplot1(4,5,'Gap',[0.01 0.01],'XTickL','Margin','YTickL','Margin','FontS',10);
    
        for i = 1:20 % numero de graficas por figura
               ti = inn(wi,i);
               subplot1(i); 

               for ii = 1:6
               y = data(ii).RoF(ti,:); 
                   if sum(isnan(y))>1
                   plot(x(1:length(y)),y,'.','Color',clrs(ii));hold on
                   data(ii).scoeff(ti,1:2)= [NaN NaN];
                   else
                   [brob,stats] = robustfit(x(1:length(y)),y);
                   plot(x(1:length(y)),y,'.','Color',clrs(ii));hold on
                   plot(x,brob(1)+brob(2)*x,'color',clrs(ii),'LineWidth',2)
                   text(3,-5+ii, num2str(stats.se(2),2),'color',clrs(ii))
                   data(ii).scoeff(ti,1:2)= brob';
                   data(ii).serr(ti,1)= stats.se(2)';
                   end
               end
               
               
                
                   text(3,-6,data(ii).orf(ti),'BackgroundColor','y')
                   ylim([-7 3])
                   xlim([0 max(generations)+3])

               if i == 16
               xlabel('Generations')
               ylabel('#Cells log(RFP/YFP)')
               else
               end
         end
   end

end

%% Plot with standard error of the mean

alls = [];allerr = [];
complex = names(:,4);
prot = find(strcmp('Proteasome',complex))
RNA = find(strcmp('RNApol',complex))
retr = find(strcmp('retromer',complex))
ref = find(strcmp('ref',complex))
neutral= find(strcmp('neutral',complex))

for i = 1:6
   temp = data(i).scoeff(:,2);
   alls = [alls temp]
   temp = data(i).serr;
   allerr = [allerr temp]
    
    
end 

figure(100);clf
fts = 16
temp = alls(:,1:3);
x = nanmean(temp,2);
SEM_x =  std(temp,0,2)/sqrt(3);
temp = alls(:,4:6);
y = nanmean(temp,2);
SEM_y =  nanmean(allerr(:,4:6),2);
plot(x(prot),y(prot),'ro','markerfacecolor','r');hold on
plot(x(RNA),y(RNA),'go','markerfacecolor','g');
plot(x(retr),y(retr),'bo','markerfacecolor','b');
plot(x(ref),y(ref),'ko','markerfacecolor','k');


plot([-.1 .2],[-.1 .2],'k:')
axis square
xlim([-.1 .2])
ylim([-.1 .2])
xlabel('SC,scoeff','fontsize',fts)
ylabel('Peroxide,scoeff','fontsize',fts)
%%
cols = {'SC1','SC2','SC3','Peroxide1','Peroxide2','Peroxide3'}
labels = {'Proteasome', 'RNApol','retromer','refs','neutral'};

figure(101);clf
notBoxPlotAA(alls)
ylabel('scoeff','fontsize',fts)
set(gca,'xtick',1:6,'xticklabel',cols,'fontsize',fts)

temp = nan(80,5);
temp(prot,1) = x(prot);
temp(RNA,2) = x(RNA);
temp(retr,3) = x(retr);
tt = alls(ref(1:4),1:3); tt = tt(:);
temp(1:length(tt),4) = tt;
tt = alls(neutral,1:3); tt = tt(:);
temp(1:length(tt),5) = tt;

temp2 = nan(80,5);
temp2(prot,1) = y(prot);
temp2(RNA,2) = y(RNA);
temp2(retr,3) = y(retr);
tt = alls(retr(1:4),4:6); tt = tt(:);
temp2(1:length(tt),4) = tt';
tt = alls(neutral,4:6); tt = tt(:);
temp2(1:length(tt),5) = tt;


pos1 = (1:5)-0.2;
pos2 = (1:5)+0.2;

fts = 18;
cmap = linspecer(4,'qualitative');
gray = [.6 .6 .6];
figure(102);clf
h1 = notBoxPlotAA(temp,pos1,.15,'patch',.5,cmap(1,:),gray,'k',5);hold on
h2 = notBoxPlotAA(temp2,pos2,.15,'patch',.5,cmap(2,:),gray,'k',5);

ylabel('scoeff','fontsize',fts)
set(gca,'xtick',1:5,'xticklabel',labels,'fontsize',fts)
text(4.5,.17,'SC','fontsize',fts,'color',cmap(1,:))
text(4.5,.16,'Peroxide','fontsize',fts,'color',cmap(2,:))


xi = x;
%% all data in boxplot
temp = nan(80,5);
temp(prot,1) = xi(prot);
temp(RNA,2) = xi(RNA);
temp(retr,3) = xi(retr);
tt = alls(ref,1:3); tt = tt(:);
temp(1:length(tt),4) = tt';
tt = alls(neutral,1:3); tt = tt(:);
temp(1:length(tt),5) = tt';

fts = 14;
cmap = linspecer(2,'qualitative');

figure(103);clf
notBoxPlotAA(temp,1:5,.3,'sdline',.5,cmap(1,:),[1 1 1],cmap(2,:),5);hold on

xl = [1:5]-0.2; yl = repmat(-.042,1,5);
ylabel('S coeff','fontsize',fts)
set(gca,'xtick',1:5,'xticklabel','','fontsize',fts)
text(xl,yl,labels,'rotation',-25,'fontsize',fts)
text(4.5,.17,'SC','fontsize',fts,'color',cmap(1,:))
text(4.5,.16,'Peroxide','fontsize',fts,'color',cmap(2,:))
ylim([-.04 .08])
