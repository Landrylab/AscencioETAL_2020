%% Normalize selection coeffients and z-scores for fluorometry based competition experiments
%Normalize data with the each plate media in normal condition
controls=[12 91 63 27 55 80];
x96=1:96;
aa=setdiff(x96,controls);
load dataAll_nominal
s_norm=[];s=[];sc_norm=[];sc=[];
for i=1:10
    sTemp=bgdataAll(i).s4(1,aa);%mutant strains aprox 90 per plate
    mTemp=nanmedian(sTemp);
    snormTemp=sTemp-mTemp;
    s=[s sTemp];
    s_norm=[s_norm;snormTemp];
    
    scTemp=bgdataAll(i).s4(1,controls);%control strain aprox 6 per plate
    mcTemp=nanmedian(scTemp);
    scnormTemp=scTemp-mcTemp;
    sc=[sc scTemp];
    sc_norm=[sc_norm;scnormTemp];
    
end 

%for the reposition PL11
    sTemp=bgdataAll(11).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp11=sTemp-mTemp;
    s11_norm=snormTemp11;
 %normalization for the two plates of wild types
    sTemp=bgdataAll(12).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp12=sTemp%-mTemp;
    s12_norm=snormTemp12;
    sTemp=bgdataAll(13).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp13=sTemp-mTemp;
    s13_norm=snormTemp13;
    
%Join all the data in a single matrix
sn=nan(13,96);
sn(1:10,aa)=s_norm;
sn(1:10,controls)=sc_norm;
sn(11,:)=s11_norm; sn(12,:)=s12_norm;sn(13,:)=s13_norm;
   

%
% Normalize data with each plate media for NaCl data 

load dataAll_nacl
s_norm=[];s=[];sc_norm=[];sc=[];
for i=1:10
    sTemp=bgdataAll(i).s4(1,aa);%mutant strains aprox 90 per plate
    mTemp=nanmedian(sTemp);
    snormTemp=sTemp-mTemp;
    s=[s sTemp];
    s_norm=[s_norm;snormTemp];
    
    scTemp=bgdataAll(i).s4(1,controls);%control strain 6 per plate
    mcTemp=nanmedian(scTemp);
    scnormTemp=scTemp-mcTemp;
    sc=[sc scTemp];
    sc_norm=[sc_norm;scnormTemp];
    
end 
%for the reposition PL11
    sTemp=bgdataAll(11).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp11=sTemp-mTemp;
    s11_norm=snormTemp11;
 %normalization for the two plates of wild types
    sTemp=bgdataAll(12).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp12=sTemp-mTemp;
    s12_norm=snormTemp12;
    sTemp=bgdataAll(13).s4(1,:);
    mTemp=nanmedian(sTemp);
    snormTemp13=sTemp-mTemp;
    s13_norm=snormTemp13;
    
%Join all the data in a single matrix
snacl=nan(13,96);
snacl(1:10,aa)=s_norm;
snacl(1:10,controls)=sc_norm;
snacl(11,:)=s11_norm; snacl(12,:)=s12_norm;snacl(13,:)=s13_norm;
snacl_exp2=snacl;sn_exp2=sn;



%% Calculate z-scores.
Zt=3; % El Z-score cutoff.

%Calculate Zscores for nominal

figure(4);clf
subplot(2,2,1)
sc=sn(1:10,controls);
sx=sn(1:10,aa);
sdtc=nanstd(sc(:));
ZscoreN=sx(:)/sdtc;

bin=-20:1:20;
hx=histc(ZscoreN(:),bin);
b=bar(bin,hx,'Facecolor',[0.5 0.5 0.5],'barwidth',1);
hold on
xlim([-25 20])
plot([-4.7619 -4.7619],[0 150],'r--');plot([4.7619 4.7619],[0 150],'r--')
title('Nominal, Zscore')
axis square
grid on

temp=ZscoreN>Zt;
beneficialN=ZscoreN(temp);
temp=ZscoreN<-Zt;
deleteriousN=ZscoreN(temp);
freq1=[length(ZscoreN(:)) length(beneficialN) length(deleteriousN)];

subplot(2,2,3)
        explode=[1 1 1];
        h=pie(freq1,explode);
        cmap=[0 0 0;0 0.7 0.9;0.9 0 0.4];
        colormap(cmap)

        textObjs = findobj(h,'Type','text');
        oldStr = get(textObjs,{'String'});
        val = get(textObjs,{'Extent'});
        oldExt = cat(1,val{:});
        Names = {'Neutral: ';'Beneficial: ';'Deleterious: '};
        newStr = strcat(Names,oldStr);
        set(textObjs,{'String'},newStr);

%for NaCl
sc=snacl(1:10,controls);
sx=snacl(1:10,aa);
sdtcN=nanstd(sc(:));
ZscoreNacl=sx(:)/sdtcN;

temp=ZscoreNacl>Zt;
beneficialNacl=ZscoreNacl(temp);
temp=ZscoreNacl<-Zt;
deleteriousNacl=ZscoreNacl(temp);
freq1=[length(ZscoreNacl(:)) length(beneficialNacl) length(deleteriousNacl)];
%% Load genenames in a similar matrix

orfs=[];
for i=1:13
    fin=i*96;
    in=fin-95;
    temp=temptext(in:fin);
    orfs=[orfs;temp'];
      
end 

genenames=orfs;
    for ii=1:13
        orfnames=orfs(ii,:);
        tempp=[];
    for i=1:length(orfnames)
%         if strcmp('WT',orfnames(i));
%             tmp='WT';
%         elseif strcmp('empty',orfnames(i))|strcmp('EMPTY',orfnames(i));
%             tmp='empty';
%         else
            f=find(strcmp(orfnames(i),namesdata(:,1)));
            if ~isempty(namesdata(f,2))
                tmp=namesdata(f,2);
                        else
                tmp=orfnames(i);
            end
    %    end
        tempp=[tempp tmp];
    end
    genenames(ii,:)=tempp;
    end 
    orfs_exp2=orfs;
    
    save data_exp2 orfs_exp2 sn_exp2 snacl_exp2 orfs ZscoreN ZscoreNacl


%% find genes extremely neutral
temp_s = sn_exp2(1,:); temp_names = orfs_exp2(1,:);
ff = find(~isnan(temp_s));
s = temp_s(ff); snames = temp_names(ff);
[y,I] = sort(s); nn = snames(I);
x= 1/length(y):1/(length(y)):1;
temp_wt = find(strcmp('WT',nn));
ff1 = find(y<0.0005);ff2 = find(y>-0.0005); ff = intersect(ff1,ff2);

figure(100);clf
    scatter(x,y,12,y,'filled','marker','s');hold on
    plot(x(temp_wt),y(temp_wt),'ko')
    plot(x(ff),y(ff),'bo')
    legend('all strains','wt','neutral')
    ylim([-0.25,0.15])



ff1 = find(temp_s<0.0005);ff2 = find(temp_s>-0.0005); ff = intersect(ff1,ff2);
neutral_genes_plt1 = pos(ff,:)
temp = temp_names(ff)'
neutral_s = temp_s(ff)'
neutral_sna = snacl_exp2(ff)'

ylim([-0.1,0.1])



