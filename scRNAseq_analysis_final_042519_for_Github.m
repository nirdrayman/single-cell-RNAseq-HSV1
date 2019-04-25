tic
%% add folder to path
addpath(genpath('M:\Nir Drayman\Papers\19-Transcriptmics of HSV1 infection\scripts\eLife_2019_final'));

%% Data upload and pre-processing

clearvars; close all; clc;
cd('M:\Nir Drayman\Papers\19-Transcriptmics of HSV1 infection\scripts\eLife_2019_final');
load Files_for_analysis_intitation.mat

mock=mock_orig; wt=wt_orig; d0=d0_orig;
names_mock=gene_names_mock_orig; names_wt=gene_names_wt_orig; names_d0=gene_names_d0_orig;

% Filter out cells will less than 2000 or more than 10,000 UMI
min_num_transcripts=2000; max_num_transcripts=1e4;
remove_cells=sum(mock)<min_num_transcripts | sum(mock) >max_num_transcripts;
mock(:,remove_cells)=[];
remove_cells=sum(wt)<min_num_transcripts | sum(wt) >max_num_transcripts;
wt(:,remove_cells)=[];
remove_cells=sum(d0)<min_num_transcripts | sum(d0) >max_num_transcripts;
d0(:,remove_cells)=[];

% Filter genes not expressed in any cell
remove_genes=sum(mock,2)==0;
mock(remove_genes,:)=[];
names_mock(remove_genes)=[];
remove_genes=sum(wt,2)==0;
wt(remove_genes,:)=[];
names_wt(remove_genes)=[];
remove_genes=sum(d0,2)==0;
d0(remove_genes,:)=[];
names_d0(remove_genes)=[];

% Filter genes expressed in less than 10 cells
remove_genes=sum(mock>0,2)<=10;
mock(remove_genes,:)=[];
names_mock(remove_genes)=[];
remove_genes=sum(wt>0,2)<=10;
wt(remove_genes,:)=[];
names_wt(remove_genes)=[];
remove_genes=sum(d0>0,2)<=10;
d0(remove_genes,:)=[];
names_d0(remove_genes)=[];

%annotate gene sets (host, mitochondrial and viral)
%mock:
mvgenes=find(contains(names_mock,'HSV1-'));
mmtgenes=find(contains(names_mock,'MT-'));
mhgenes=find(contains(names_mock,'HSV1-')==0 & contains(names_mock,'MT-')==0);
mnonvgenes=[mmtgenes;mhgenes];
mnonv_counts=sum(mock(mnonvgenes,:));
mv_counts=sum(mock(mvgenes,:));
%wt:
wvgenes=find(contains(names_wt,'HSV1-'));
wmtgenes=find(contains(names_wt,'MT-'));
whgenes=find(contains(names_wt,'HSV1-')==0 & contains(names_wt,'MT-')==0);
wnonvgenes=[wmtgenes;whgenes];
wnonv_counts=sum(wt(wnonvgenes,:));
wv_counts=sum(wt(wvgenes,:));
%dICP0:
dvgenes=find(contains(names_d0,'HSV1-'));
dmtgenes=find(contains(names_d0,'MT-'));
dhgenes=find(contains(names_d0,'HSV1-')==0 & contains(names_d0,'MT-')==0);
dnonvgenes=[dmtgenes;dhgenes];
dnonv_counts=sum(d0(dnonvgenes,:));
dv_counts=sum(d0(dvgenes,:));

%mito fraction
mock_mito = sum(mock(mmtgenes,:)) ./ mnonv_counts;
wt_mito = sum(wt(wmtgenes,:)) ./ wnonv_counts;
d0_mito = sum(d0(dmtgenes,:)) ./ dnonv_counts;

%viral fraction
mock_viral = sum(mock(mvgenes,:)) ./ sum(mock);
wt_viral = sum(wt(wvgenes,:)) ./ sum(wt);
d0_viral = sum(d0(dvgenes,:)) ./ sum(d0);

% Filter cells with high mito fraction
remove_cells=mock_mito>0.2;
mock(:,remove_cells)=[];
mock_mito(remove_cells)=[];
mock_viral(remove_cells)=[];

remove_cells=wt_mito>0.2;
wt(:,remove_cells)=[];
wt_mito(remove_cells)=[];
wt_viral(remove_cells)=[];

remove_cells=d0_mito>0.4;
d0(:,remove_cells)=[];
d0_mito(remove_cells)=[];
d0_viral(remove_cells)=[];

% Filter genes not in dataset
remove_genes=sum(mock,2)==0;
mock(remove_genes,:)=[];
names_mock(remove_genes)=[];
remove_genes=sum(wt,2)==0;
wt(remove_genes,:)=[];
names_wt(remove_genes)=[];
remove_genes=sum(d0,2)==0;
d0(remove_genes,:)=[];
names_d0(remove_genes)=[];

% Re-annotate gene sets
%mock:
mvgenes=find(contains(names_mock,'HSV1-'));
mmtgenes=find(contains(names_mock,'MT-'));
mhgenes=find(contains(names_mock,'HSV1-')==0 & contains(names_mock,'MT-')==0);
mnonvgenes=[mmtgenes;mhgenes];
mnonv_counts=sum(mock(mnonvgenes,:));
mv_counts=sum(mock(mvgenes,:));
%wt:
wvgenes=find(contains(names_wt,'HSV1-'));
wmtgenes=find(contains(names_wt,'MT-'));
whgenes=find(contains(names_wt,'HSV1-')==0 & contains(names_wt,'MT-')==0);
wnonvgenes=[wmtgenes;whgenes];
wnonv_counts=sum(wt(wnonvgenes,:));
wv_counts=sum(wt(wvgenes,:));
%d0:
dvgenes=find(contains(names_d0,'HSV1-'));
dmtgenes=find(contains(names_d0,'MT-'));
dhgenes=find(contains(names_d0,'HSV1-')==0 & contains(names_d0,'MT-')==0);
dnonvgenes=[dmtgenes;dhgenes];
dnonv_counts=sum(d0(dnonvgenes,:));
dv_counts=sum(d0(dvgenes,:));

%remove cells that now have less than 2000 transcripts  
remove_cells=sum(mock)<min_num_transcripts;
mock(:,remove_cells)=[];
mock_mito(remove_cells)=[];
mock_viral(remove_cells)=[];
remove_cells=sum(wt)<min_num_transcripts;
wt(:,remove_cells)=[];
wt_mito(remove_cells)=[];
wt_viral(remove_cells)=[];
remove_cells=sum(d0)<min_num_transcripts;
d0(:,remove_cells)=[];
d0_mito(remove_cells)=[];
d0_viral(remove_cells)=[];

%Re-annotate gene sets
%mock:
mvgenes=find(contains(names_mock,'HSV1-'));
mmtgenes=find(contains(names_mock,'MT-'));
mhgenes=find(contains(names_mock,'HSV1-')==0 & contains(names_mock,'MT-')==0);
mnonvgenes=[mmtgenes;mhgenes];
mnonv_counts=sum(mock(mnonvgenes,:));
mv_counts=sum(mock(mvgenes,:));
%wt:
wvgenes=find(contains(names_wt,'HSV1-'));
wmtgenes=find(contains(names_wt,'MT-'));
whgenes=find(contains(names_wt,'HSV1-')==0 & contains(names_wt,'MT-')==0);
wnonvgenes=[wmtgenes;whgenes];
wnonv_counts=sum(wt(wnonvgenes,:));
wv_counts=sum(wt(wvgenes,:));
%d0:
dvgenes=find(contains(names_d0,'HSV1-'));
dmtgenes=find(contains(names_d0,'MT-'));
dhgenes=find(contains(names_d0,'HSV1-')==0 & contains(names_d0,'MT-')==0);
dnonvgenes=[dmtgenes;dhgenes];
dnonv_counts=sum(d0(dnonvgenes,:));
dv_counts=sum(d0(dvgenes,:));

%%%% Done Pre-processing data %%%
%% Analysis of HSV-1 viral genes

%plots for Figure 1
%violin plot - % viral transcripts in cells
%Note - distributionPlot is a function that makes Viloin plots and can be downlaoded from: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m

%Fig 1E
data={log10(mock_viral*1e4+1),log10(wt_viral*1e4+1),log10(d0_viral*1e4+1)};
figure; distributionPlot(data,'showMM',0,'histOpt',1.1)
h=get(gca);
h.Children(1).FaceColor=colors(2,:);
h.Children(2).FaceColor=colors(1,:);
ylabel 'HSV-1 transcripts (%)'
xticklabels({'mock','wt','\DeltaICP0'});
ylim([0,4]);
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
yticks(0:1:4);
yticklabels({'0','0.1','1','10','100'})
title 'HSV-1 gene expression in single cells'

%Fig 1F - tSNE of only viral genes, wt infection
rng(1) % for reproducibility
vwt=wt(wvgenes,:);
Y=tsne(vwt');
figure; scatter(-Y(:,1),-Y(:,2),40,log10(wt_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k';
h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
title 'tSNE based on viral gene expresion (wt)'

dwt=d0(dvgenes,:);
rng(1) % for reproducibility
Y2=tsne(dwt');
figure; scatter(Y2(:,1),Y2(:,2),40,log10(d0_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k';
h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
title 'tSNE based on viral gene expresion (\DeltaICP0)'

clear h

% classes of viral genes
IE={'HSV1-UL54','HSV1-RS1','HSV1-US1', 'HSV1-RL2_1', 'HSV1-RL2_2', 'HSV1-RL2_3', 'HSV1-US12'};
E={'HSV1-UL2','HSV1-UL5','HSV1-UL8','HSV1-UL9','HSV1-UL12','HSV1-UL23',...
    'HSV1-UL29','HSV1-UL30','HSV1-UL39','HSV1-UL40','HSV1-UL42','HSV1-UL50','HSV1-UL52'};
L1={'HSV1-RL1','HSV1-UL47','HSV1-UL46','HSV1-UL37','HSV1-UL51','HSV1-UL16',...
    'HSV1-UL20','HSV1-UL27','HSV1-UL34','HSV1-US8','HSV1-US10','HSV1-UL7',...
    'HSV1-UL18','HSV1-UL24','HSV1-US7','HSV1-US3','HSV1-US8A','HSV1-US4','HSV1-UL41'};
L2={'HSV1-UL36','HSV1-UL38','HSV1-UL25','HSV1-UL49A','HSV1-UL32','HSV1-UL44',...
    'HSV1-UL33','HSV1-UL45','HSV1-UL14','HSV1-US11','HSV1-UL3','HSV1-UL10',...
    'HSV1-UL56','HSV1-UL4','HSV1-US2','HSV1-UL35','HSV1-UL55'};

wt_IE=zeros(1,size(vwt,2));
for ii=1:length(IE)
    gene=find(strcmp(names_wt(wvgenes),IE(ii)));
    if ~isempty(gene)
    wt_IE=wt_IE+vwt(gene,:);
    end
end
wt_IE=wt_IE./sum(wt);

wt_E=zeros(1,size(vwt,2));
for ii=1:length(E)
    gene=find(strcmp(names_wt(wvgenes),E(ii)));
    if ~isempty(gene)
    wt_E=wt_E+vwt(gene,:);
    end
end
wt_E=wt_E./sum(wt);

wt_L1=zeros(1,size(vwt,2));
for ii=1:length(L1)
    gene=find(strcmp(names_wt(wvgenes),L1(ii)));
    if ~isempty(gene)
    wt_L1=wt_L1+vwt(gene,:);
    end
end
wt_L1=wt_L1./sum(wt);

wt_L2=zeros(1,size(vwt,2));
for ii=1:length(L2)
    gene=find(strcmp(names_wt(wvgenes),L2(ii)));
    if ~isempty(gene)
    wt_L2=wt_L2+vwt(gene,:);
    end
end
wt_L2=wt_L2./sum(wt);

%d0
d0_IE=zeros(1,size(dwt,2));
for ii=1:length(IE)
    gene=find(strcmp(names_d0(dvgenes),IE(ii)));
    if ~isempty(gene)
    d0_IE=d0_IE+dwt(gene,:);
    end
end
d0_IE=d0_IE./sum(d0);

d0_E=zeros(1,size(dwt,2));
for ii=1:length(E)
    gene=find(strcmp(names_d0(dvgenes),E(ii)));
    if ~isempty(gene)
    d0_E=d0_E+dwt(gene,:);
    end
end
d0_E=d0_E./sum(d0);

d0_L1=zeros(1,size(dwt,2));
for ii=1:length(L1)
    gene=find(strcmp(names_d0(dvgenes),L1(ii)));
    if ~isempty(gene)
    d0_L1=d0_L1+dwt(gene,:);
    end
end
d0_L1=d0_L1./sum(d0);

d0_L2=zeros(1,size(dwt,2));
for ii=1:length(L2)
    gene=find(strcmp(names_d0(dvgenes),L2(ii)));
    if ~isempty(gene)
    d0_L2=d0_L2+dwt(gene,:);
    end
end
d0_L2=d0_L2./sum(d0);

%Fraction of classes
wt_classes=[wt_IE;wt_E;wt_L1;wt_L2];
wt_classes=wt_classes';
wt_sum=sum(wt_classes,2);
wt_viral2=wt_sum;
[~,ind_wt]=sort(wt_viral2);

%find cells with at least 15 viral counts
start_wt=find(wt_viral2(ind_wt).*sum(wt)'>=15,1);

%Fig 1G
figure;
bar(wt_classes(ind_wt(start_wt:end),:)./wt_sum(ind_wt(start_wt:end)),'stacked')
ylim([0 1]);
xticks('');
xlabel 'Infected cells (sorted by % viral transcripts)'
ylabel '% of viral transcripts'
yticks(0:0.25:1)
yticklabels(0:25:100)
hold;
plot(1:length(wt_viral2)-start_wt+1,wt_viral2(ind_wt(start_wt:end)),'k','LineWidth',2)
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
title 'viral gene classes (wt)'

%Fraction of classes d0
d0_classes=[d0_IE;d0_E;d0_L1;d0_L2];
d0_classes=d0_classes';
d0_sum=sum(d0_classes,2);
d0_viral2=d0_sum;
[~,ind_d0]=sort(d0_viral2);
%find cells with at least 15 viral counts
start_d0=find(d0_viral2(ind_d0).*sum(d0)'>=15,1);

%Fig S2B
figure;
bar(d0_classes(ind_d0(start_d0:end),:)./d0_sum(ind_d0(start_d0:end)),'stacked')
ylim([0 1]);
xticks('');
xlabel 'Infected cells (sorted by % viral transcripts)'
ylabel '% of viral transcripts'
yticks(0:0.25:1)
yticklabels(0:25:100)
hold;
plot(1:length(d0_viral2)-start_d0+1,d0_viral2(ind_d0(start_d0:end)),'k','LineWidth',2)
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
title 'viral gene classes (\DeltaICP0)'


% Fig 1H, Fig S2C - correlations between viral gene clasees and total viral gene expression

figure;
subplot(2,4,1)
scatter(wt_viral2(ind_wt(start_wt:end)),wt_classes(ind_wt(start_wt:end),1)./wt_sum(ind_wt(start_wt:end)),40,'filled','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'IE genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20); 
subplot(2,4,2)
scatter(wt_viral2(ind_wt(start_wt:end)),wt_classes(ind_wt(start_wt:end),2)./wt_sum(ind_wt(start_wt:end)),40,'filled','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'E genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);
subplot(2,4,3)
scatter(wt_viral2(ind_wt(start_wt:end)),wt_classes(ind_wt(start_wt:end),3)./wt_sum(ind_wt(start_wt:end)),40,'filled','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'L1 genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);
subplot(2,4,4)
scatter(wt_viral2(ind_wt(start_wt:end)),wt_classes(ind_wt(start_wt:end),4)./wt_sum(ind_wt(start_wt:end)),40,'filled','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'L2 genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);

subplot(2,4,5)
scatter(d0_viral2(ind_d0(start_d0:end)),d0_classes(ind_d0(start_d0:end),1)./d0_sum(ind_d0(start_d0:end)),40,'filled','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'IE genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20); 
subplot(2,4,6)
scatter(d0_viral2(ind_d0(start_d0:end)),d0_classes(ind_d0(start_d0:end),2)./d0_sum(ind_d0(start_d0:end)),40,'filled','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'E genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);
subplot(2,4,7)
scatter(d0_viral2(ind_d0(start_d0:end)),d0_classes(ind_d0(start_d0:end),3)./d0_sum(ind_d0(start_d0:end)),40,'filled','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'L1 genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);
subplot(2,4,8)
scatter(d0_viral2(ind_d0(start_d0:end)),d0_classes(ind_d0(start_d0:end),4)./d0_sum(ind_d0(start_d0:end)),40,'filled','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
set(gca,'XScale','log'); xlim([0.003,0.5]); xticks([0.01,0.03,0.1,0.3]); xticklabels({'1','3','10','30'});
xlabel 'HSV-1 trancripts (%)'; ylabel 'L2 genes (Fraction)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20);

title 'Correlation between viral gene classes and total viral gene expression'

%Pearson correlations
a=[]; p=[];
[a(1),p(1)]=corr(log10(wt_viral2(ind_wt(start_wt:end))),wt_classes(ind_wt(start_wt:end),1)./wt_sum(ind_wt(start_wt:end)));
[a(2),p(2)]=corr(log10(wt_viral2(ind_wt(start_wt:end))),wt_classes(ind_wt(start_wt:end),2)./wt_sum(ind_wt(start_wt:end)));
[a(3),p(3)]=corr(log10(wt_viral2(ind_wt(start_wt:end))),wt_classes(ind_wt(start_wt:end),3)./wt_sum(ind_wt(start_wt:end)));
[a(4),p(4)]=corr(log10(wt_viral2(ind_wt(start_wt:end))),wt_classes(ind_wt(start_wt:end),4)./wt_sum(ind_wt(start_wt:end)));
[a(5),p(5)]=corr(log10(d0_viral2(ind_d0(start_d0:end))),d0_classes(ind_d0(start_d0:end),1)./d0_sum(ind_d0(start_d0:end)));
[a(6),p(6)]=corr(log10(d0_viral2(ind_d0(start_d0:end))),d0_classes(ind_d0(start_d0:end),2)./d0_sum(ind_d0(start_d0:end)));
[a(7),p(7)]=corr(log10(d0_viral2(ind_d0(start_d0:end))),d0_classes(ind_d0(start_d0:end),3)./d0_sum(ind_d0(start_d0:end)));
[a(8),p(8)]=corr(log10(d0_viral2(ind_d0(start_d0:end))),d0_classes(ind_d0(start_d0:end),4)./d0_sum(ind_d0(start_d0:end)));


%% Analysis of Host and viral genes

% count UMI per cell
wt_umi=sum(wt);
d0_umi=sum(d0);

%log normalize data (divide each gene by the cell's total UMI. Multiply by 10,000, add 1 and take the log
wt=log(wt./sum(wt)*1e4+1);
d0=log(d0./sum(d0)*1e4+1);

%regress out number of UMI
wt_umi_corrected=wt;
d0_umi_corrected=d0;

for ii=1:size(wt,1)
try eq=fitlm(wt_umi',wt(ii,:)'); catch; continue; end
new_y=eq.Coefficients{1,1}+eq.Coefficients{2,1}*wt_umi;    
wt_umi_corrected(ii,:)=wt(ii,:)-new_y;
end

for ii=1:size(d0,1)
try eq=fitlm(d0_umi',d0(ii,:)'); catch; continue; end
new_y=eq.Coefficients{1,1}+eq.Coefficients{2,1}*d0_umi;    
d0_umi_corrected(ii,:)=d0(ii,:)-new_y;
end

%z-score each gene
wt_umi_corrected_z=(wt_umi_corrected-mean(wt_umi_corrected,2))./std(wt_umi_corrected,[],2);
d0_umi_corrected_z=(d0_umi_corrected-mean(d0_umi_corrected,2))./std(d0_umi_corrected,[],2);

% G2/M score
g2m={ 'HMGB2'   'CDK1'    'NUSAP1'  'UBE2C'   'BIRC5'   'TPX2'    'TOP2A'   'NDC80'   'CKS2'    'NUF2'    'CKS1B' ...  
 'MKI67'   'TMPO'    'CENPF'};
    
wt_g2m=zeros(1,size(wt_umi_corrected_z,2));
wt_binary=wt_umi_corrected_z>0;

d0_g2m=zeros(1,size(d0_umi_corrected_z,2));
d0_binary=d0_umi_corrected_z>0;

for ii=1:length(g2m)
    gene=find(strcmp(names_wt,g2m(ii)));
    if ~isempty(gene)
    wt_g2m=wt_g2m+wt_binary(gene,:);
    end
end

for ii=1:length(g2m)
    gene=find(strcmp(names_d0,g2m(ii)));
    if ~isempty(gene)
    d0_g2m=d0_g2m+d0_binary(gene,:);
    end
end

%tSNE of host+viral genes wt
%First PCA
[~,score] = pca(wt_umi_corrected_z');
%now tSNE
rng(1);
Y3=tsne(score(:,1:20));
idx = kmeans(score(:,1:20),3);

figure; title 'tSNE of host+viral genes before cell-cycle regression (wt)'
subplot(131); title 'clusters'
scatter(-Y3(:,1),Y3(:,2),40,idx,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(132); title 'viral gene expression'
scatter(-Y3(:,1),Y3(:,2),40,log10(wt_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(133); title 'G2/M score'
scatter(-Y3(:,1),Y3(:,2),40,wt_g2m,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)

%tSNE of host+viral genes dICP0
%find most correlated/anti-correlated genes with viral load

gene_corr=zeros(size(d0_umi_corrected_z,1),1);

for ii=1:size(d0_umi_corrected_z,1)
    gene_corr(ii)=corr(d0_umi_corrected_z(ii,:)',log10(d0_viral*1e4+1)');
end

[~,ind_gene_corr]=sort(gene_corr);
genes_for_tsne=[ind_gene_corr(1:100),ind_gene_corr(end-99:end)];

mat_for_tsne=d0_umi_corrected_z(genes_for_tsne,:);

%tSNE of host+viral genes dICP0
%First PCA
[~,score2] = pca(mat_for_tsne');
%now tSNE
rng(1)
Y4=tsne(score2(:,1:20));
idx2 = kmeans(score2(:,1:10),4);

% The effect of the cell-cycle on viral gene expression

%viral load vs cell-cycle
bins2=linspace(min(wt_g2m),max(wt_g2m),14);
bins3=linspace(min(d0_g2m),max(d0_g2m),14);

group2=[]; load2=[]; percent_infected2=[];
group3=[]; load3=[]; percent_infected3=[];

for ii=1:length(bins2)-1
    group2{ii}=find(wt_g2m>=bins2(ii) & wt_g2m<bins2(ii+1));
    load2(ii,1)=mean(wt_viral(group2{ii})*100);
    load2(ii,2)=std(wt_viral(group2{ii})*100)/sqrt(length(group2{ii})); 
    percent_infected2(ii,1)=sum(wt_viral(group2{ii})>1e-3)./length(group2{ii});
end

for ii=1:length(bins3)-1
    group3{ii}=find(d0_g2m>=bins3(ii) & d0_g2m<bins3(ii+1));
    load3(ii,1)=mean(d0_viral(group3{ii})*100);
    load3(ii,2)=std(d0_viral(group3{ii})*100)/sqrt(length(group3{ii})); 
    percent_infected3(ii,1)=sum(d0_viral(group3{ii})>1e-3)./length(group3{ii});
end

% Fig S3C,D - correlation between viral gene expression and G2_M score.
% Note that to make plots you will need the boundedline function from: https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m  

figure; title 'cell cycle effect on HSV-1 gene expression (wt)'
boundedline(bins2(1:end-1),smooth(load2(:,1)),smooth(load2(:,2)));
h=get(gca,'Children');
h(1).Color='k'; h(1).LineWidth=2; h(2).FaceColor='k';
h(2).FaceAlpha=0.2;
xlabel 'G2/M score'; ylabel 'HSV-1 viral transcripts (%)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([-0.5 12.5])
position=[680   455   592   523];
h=get(gca,'Parent'); h.Position=position;
hold;
scatter(wt_g2m,wt_viral*100,3,'filled','MarkerFaceColor','none', 'MarkerEdgeColor','k')
h.Children.YScale='log';
h.Children.YTickLabel={0.01,0.1,1,10,100};

figure; title 'cell cycle effect on HSV-1 gene expression (\DeltaICP0)'
boundedline(bins3(1:end-1),smooth(load3(:,1)),smooth(load3(:,2)));
h=get(gca,'Children');
h(1).Color='k'; h(1).LineWidth=2; h(2).FaceColor='k';
h(2).FaceAlpha=0.2;
xlabel 'G2/M score'; ylabel 'HSV-1 viral transcripts (%)'
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([-0.5 12.5])
position=[680   455   592   523];
h=get(gca,'Parent'); h.Position=position;
hold;
scatter(d0_g2m,d0_viral*100,3,'filled','MarkerFaceColor','none', 'MarkerEdgeColor','k')
h.Children.YScale='log';
ylim([0 30])
yticks([0.01 0.1 1 10]);
yticks([0.01 0.1 1 10]);

%regress out G2M score
wt_g2m_corrected=wt_umi_corrected_z;
for ii=1:size(wt_umi_corrected_z,1)
try eq=fitlm(wt_g2m',wt_umi_corrected_z(ii,:)'); catch; continue; end
new_y=eq.Coefficients{1,1}+eq.Coefficients{2,1}*wt_g2m;    
wt_g2m_corrected(ii,:)=wt_umi_corrected_z(ii,:)-new_y;
end

d0_g2m_corrected=d0_umi_corrected_z;
for ii=1:size(d0_umi_corrected_z,1)
try eq=fitlm(d0_g2m',d0_umi_corrected_z(ii,:)'); catch; continue; end
new_y=eq.Coefficients{1,1}+eq.Coefficients{2,1}*d0_g2m;    
d0_g2m_corrected(ii,:)=d0_umi_corrected_z(ii,:)-new_y;
end
 
%z-score each gene
wt_g2m_corrected_z=(wt_g2m_corrected-mean(wt_g2m_corrected,2))./std(wt_g2m_corrected,[],2);
d0_g2m_corrected_z=(d0_g2m_corrected-mean(d0_g2m_corrected,2))./std(d0_g2m_corrected,[],2);

%tSNE of host+viral genes
%First PCA
[~,score3] = pca(wt_g2m_corrected_z');
%now tSNE
rng(1)
Y5=tsne(score3(:,1:20));
idx5 = kmeans(score3(:,1:20),2);

%Fig S3A,B - host+viral tSNE of wt, before and after cell-cycle regression
figure; title 'tSNE host+viral, wt, before and after cell-cycle regression'
subplot(2,3,1); title 'clusters'
scatter(-Y3(:,1),Y3(:,2),40,idx,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(2,3,2); title 'viral gene expression'
scatter(-Y3(:,1),Y3(:,2),40,log10(wt_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(2,3,3); title 'G2/M score'
scatter(-Y3(:,1),Y3(:,2),40,wt_g2m,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(2,3,4)
scatter(-Y5(:,1),-Y5(:,2),40,idx,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(2,3,5)
scatter(-Y5(:,1),-Y5(:,2),40,log10(wt_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(2,3,6)
scatter(-Y5(:,1),-Y5(:,2),40,wt_g2m,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)


%tSNE of host+viral genes dICP0
%find most correlated/anti-correlated genes with viral load

gene_corr=zeros(size(d0_g2m_corrected_z,1),1);

for ii=1:size(d0_g2m_corrected_z,1)
    gene_corr(ii)=corr(d0_g2m_corrected_z(ii,:)',log10(d0_viral*1e4+1)');
end

[~,ind_gene_corr]=sort(gene_corr);
genes_for_tsne=[ind_gene_corr(1:100),ind_gene_corr(end-99:end)];

mat_for_tsne=d0_g2m_corrected_z(genes_for_tsne,:);

%First PCA
[~,score4] = pca(mat_for_tsne');
%now tSNE
rng(1)
Y6=tsne(score4(:,1:20));
idx6 = kmeans(score4(:,1:10),4);

figure; title 'tSNE host+viral, \DeltaICP0, after cell-cycle regression'
subplot(131); title 'viral gene expression'
scatter(Y6(:,1),Y6(:,2),40,log10(d0_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(132); title 'G2/M score'
scatter(Y6(:,1),Y6(:,2),40,d0_g2m,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(133); title 'clusters'
scatter(Y6(:,1),Y6(:,2),40,idx6,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Specific genes expression

%anti-viral, wt
wt_cluster_id=zeros(length(idx5),3);
for ii=1:size(wt_cluster_id,1)
if idx5(ii)==1
    wt_cluster_id(ii,:)=colors(4,:);
else
wt_cluster_id(ii,:)=colors(5,:);
end
end

ifit2=wt(find(strcmp(names_wt,'IFIT2')),:); ifit2_expressd=ifit2>0;
ifit3=wt(find(strcmp(names_wt,'IFIT3')),:); ifit3_expressd=ifit3>0;

%Fig 2A,B,C
figure; title 'tSNE viral+host, wt'
subplot(141); title 'clusters'
scatter(-Y5(:,1),-Y5(:,2),40,wt_cluster_id,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(142); title 'viral gene expression'
scatter(-Y5(:,1),-Y5(:,2),40,log10(wt_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(143); title 'IFIT2 expression'
scatter(-Y5(:,1),-Y5(:,2),40,[0.8 0.8 0.8],'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
hold;
scatter(-Y5(ifit2_expressd,1),-Y5(ifit2_expressd,2),40,'r','filled'); colormap('jet')
h=get(gca); h.Children(1).MarkerEdgeColor='k'; h.Children(1).MarkerEdgeAlpha=0.2; h.Children(1).MarkerFaceAlpha=0.8;
xticks(''); yticks('');  set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(144); title 'IFIT3 expression'
scatter(-Y5(:,1),-Y5(:,2),40,[0.8 0.8 0.8],'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
hold;
scatter(-Y5(ifit3_expressd,1),-Y5(ifit3_expressd,2),40,'r','filled'); colormap('jet')
h=get(gca); h.Children(1).MarkerEdgeColor='k'; h.Children(1).MarkerEdgeAlpha=0.2; h.Children(1).MarkerFaceAlpha=0.8;
xticks(''); yticks(''); set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)

%re-order the cluster in d0
idx7=zeros(length(idx6),1);
idx7(idx6==1)=2;
idx7(idx6==2)=3;
idx7(idx6==3)=4;
idx7(idx6==4)=1;
idx6=idx7;

%anti-viral, d0
d0_cluster_id=zeros(length(idx6),3);
for ii=1:size(d0_cluster_id,1)
if idx6(ii)==1
    d0_cluster_id(ii,:)=colors(1,:);
elseif idx6(ii)==2
    d0_cluster_id(ii,:)=colors(2,:);
elseif idx6(ii)==3
    d0_cluster_id(ii,:)=colors(3,:);    
else    
d0_cluster_id(ii,:)=colors(4,:);
end
end

ifit1_d0=d0(find(strcmp(names_d0,'IFIT1')),:); ifit1_d0_expressd=ifit1_d0>0;
mx2_d0=d0(find(strcmp(names_d0,'MX2')),:); mx2_d0_expressd=mx2_d0>0;
oasl_d0=d0(find(strcmp(names_d0,'OASL')),:); oasl_d0_expressd=oasl_d0>0;


%tSNE of d0 infection (Fig. 3A,B,C)
figure; title 'tSNE host+viral, \DeltaICP0'
subplot(151); title 'clusters'
scatter(Y6(:,1),Y6(:,2),40,d0_cluster_id,'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); xlabel 'tSNE1'; ylabel 'tSNE2'; set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(152); title 'viral gene expression'
scatter(Y6(:,1),Y6(:,2),40,log10(d0_viral*1e4+1),'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
xticks(''); yticks(''); set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(153); title 'IFIT1 expression'
scatter(Y6(:,1),Y6(:,2),40,[0.8 0.8 0.8],'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
hold;
scatter(Y6(ifit1_d0_expressd,1),Y6(ifit1_d0_expressd,2),40,'r','filled'); colormap('jet')
h=get(gca); h.Children(1).MarkerEdgeColor='k'; h.Children(1).MarkerEdgeAlpha=0.2; h.Children(1).MarkerFaceAlpha=0.8;
xticks(''); yticks('');  set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(154); title 'MX2 expression'
scatter(Y6(:,1),Y6(:,2),40,[0.8 0.8 0.8],'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
hold;
scatter(Y6(mx2_d0_expressd,1),Y6(mx2_d0_expressd,2),40,'r','filled'); colormap('jet')
h=get(gca); h.Children(1).MarkerEdgeColor='k'; h.Children(1).MarkerEdgeAlpha=0.2; h.Children(1).MarkerFaceAlpha=0.8;
xticks(''); yticks(''); set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)
subplot(155); title 'OASL expression'
scatter(Y6(:,1),Y6(:,2),40,[0.8 0.8 0.8],'filled'); colormap('jet')
h=get(gca); h.Children.MarkerEdgeColor='k'; h.Children.MarkerEdgeAlpha=0.2; h.Children.MarkerFaceAlpha=0.8;
hold;
scatter(Y6(oasl_d0_expressd,1),Y6(oasl_d0_expressd,2),40,'r','filled'); colormap('jet')
h=get(gca); h.Children(1).MarkerEdgeColor='k'; h.Children(1).MarkerEdgeAlpha=0.2; h.Children(1).MarkerFaceAlpha=0.8;
xticks(''); yticks('');  set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',20)

%Expression levels by cluster (Fig. 3D)
figure; hold; title 'gene expression by clusters, \DeltaICP0'

%viral genes
data={mock_viral,d0_viral(idx6==1),d0_viral(idx6==2),d0_viral(idx6==3),d0_viral(idx6==4)};

subplot(221); hold;
scatter(linspace(0.75,1.25,4500), data{1}*1e4+1,100,'filled','MarkerFaceColor','k', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(1.75,2.25,sum(idx6==1)), data{2}*1e4+1,100,'filled','MarkerFaceColor',colors(1,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(2.75,3.25,sum(idx6==2)), data{3}*1e4+1,100,'filled','MarkerFaceColor',colors(2,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(3.75,4.25,sum(idx6==3)), data{4}*1e4+1,100,'filled','MarkerFaceColor',colors(3,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(4.75,5.25,sum(idx6==4)), data{5}*1e4+1,100,'filled','MarkerFaceColor',colors(4,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
set(gca,'YScale','log')
xlabel 'Cluster'
ylabel '% HSV-1 transcripts'
xticks(1:5);
xticklabels({'mock','1','2','3','4'})
box('off')
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',25); xlim([0.5,5.5]);
yticks([1,10,100,1000,10000]); yticklabels({'0','0.1','1','10','100'})


%IFIT1
data={mock(find(strcmp(names_mock,'IFIT1')),:),d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),idx6==1),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),idx6==2),d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),idx6==3),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),idx6==4)};

miny=min(d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),:));
maxy=max(d0_g2m_corrected_z(find(strcmp(names_d0,'IFIT1')),:));

subplot(222); hold;
scatter(linspace(0.75,1.25,4500), zeros(4500,1),100,'filled','MarkerFaceColor','k', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(1.75,2.25,sum(idx6==1)), data{2},100,'filled','MarkerFaceColor',colors(1,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(2.75,3.25,sum(idx6==2)), data{3},100,'filled','MarkerFaceColor',colors(2,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(3.75,4.25,sum(idx6==3)), data{4},100,'filled','MarkerFaceColor',colors(3,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(4.75,5.25,sum(idx6==4)), data{5},100,'filled','MarkerFaceColor',colors(4,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
xlabel 'Cluster'
ylabel 'IFIT1 expression'
xticks(1:5);
xticklabels({'mock','1','2','3','4'})
box('off')
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',25); xlim([0.5,5.5]); ylim([miny maxy])

    %chi^2 test
    %make binray:
    low_infected=[data{2},data{3},data{4}]; 
    low_infected(low_infected>0)=1; low_infected(low_infected<=0)=0;
    high_infected=data{5};
    high_infected(high_infected>0)=1; high_infected(high_infected<=0)=0;
   [~,pval]=vartest2(low_infected,high_infected)

%MX2
data={mock(find(strcmp(names_mock,'MX2')),:),d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),idx6==1),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),idx6==2),d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),idx6==3),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),idx6==4)};

miny=min(d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),:));
maxy=max(d0_g2m_corrected_z(find(strcmp(names_d0,'MX2')),:));

subplot(223); hold;
scatter(linspace(0.75,1.25,4500), zeros(4500,1),100,'filled','MarkerFaceColor','k', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(1.75,2.25,sum(idx6==1)), data{2},100,'filled','MarkerFaceColor',colors(1,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(2.75,3.25,sum(idx6==2)), data{3},100,'filled','MarkerFaceColor',colors(2,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(3.75,4.25,sum(idx6==3)), data{4},100,'filled','MarkerFaceColor',colors(3,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(4.75,5.25,sum(idx6==4)), data{5},100,'filled','MarkerFaceColor',colors(4,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
xlabel 'Cluster'
ylabel 'MX2 expression'
xticks(1:5);
xticklabels({'mock','1','2','3','4'})
box('off')
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',25); xlim([0.5,5.5]); ylim([miny maxy])

%chi^2 test
    %make binray:
    low_infected=[data{2},data{3},data{4}]; 
    low_infected(low_infected>0)=1; low_infected(low_infected<=0)=0;
    high_infected=data{5};
    high_infected(high_infected>0)=1; high_infected(high_infected<=0)=0;
   [~,pval]=vartest2(low_infected,high_infected)
   
%OASL
data={mock(find(strcmp(names_mock,'OASL')),:),d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),idx6==1),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),idx6==2),d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),idx6==3),...
    d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),idx6==4)};

miny=min(d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),:));
maxy=max(d0_g2m_corrected_z(find(strcmp(names_d0,'OASL')),:));

subplot(224); hold;
scatter(linspace(0.75,1.25,4500), zeros(4500,1),100,'filled','MarkerFaceColor','k', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(1.75,2.25,sum(idx6==1)), data{2},100,'filled','MarkerFaceColor',colors(1,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(2.75,3.25,sum(idx6==2)), data{3},100,'filled','MarkerFaceColor',colors(2,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(3.75,4.25,sum(idx6==3)), data{4},100,'filled','MarkerFaceColor',colors(3,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
scatter(linspace(4.75,5.25,sum(idx6==4)), data{5},100,'filled','MarkerFaceColor',colors(4,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
xlabel 'Cluster'
ylabel 'OASL expression'
xticks(1:5);
xticklabels({'mock','1','2','3','4'})
box('off')
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',25); xlim([0.5,5.5]); ylim([miny maxy])
%chi^2 test
    %make binray:
    low_infected=[data{2},data{3},data{4}]; 
    low_infected(low_infected>0)=1; low_infected(low_infected<=0)=0;
    high_infected=data{5};
    high_infected(high_infected>0)=1; high_infected(high_infected<=0)=0;
   [~,pval]=vartest2(low_infected,high_infected)


%% Identify differently expressed genes between clusters

%wt infection
%differently expressed by Wilcoxson ranksum test (idx5==1 : low infeciton - cluster 2. idx5==2, high infection, cluster 1)
p=[];
for ii=1:size(wt_g2m_corrected_z,1)
    p(ii,1)=ranksum(wt_g2m_corrected_z(ii,idx5==1),wt_g2m_corrected_z(ii,idx5==2));
end
p(:,2)=mafdr(p(:,1),'BHFDR',true);

hits=find(p(:,2)<0.05);
hits_matrix=wt_g2m_corrected_z(hits,:);
names_wt_hits=names_wt(hits);
hits_p_value=p(hits,2);
hits_medians=[];
hits_medians(:,1)=median(hits_matrix(:,idx5==2),2); %cluster 1
hits_medians(:,2)=median(hits_matrix(:,idx5==1),2); %cluster 2
up_hits=find(hits_medians(:,1)>hits_medians(:,2));
down_hits=find(hits_medians(:,1)<hits_medians(:,2));

%keep only host genes
remove=contains(names_wt_hits(up_hits),'HSV1');
up_hits2=up_hits(~remove);

%keep only expressed in >10 cells in cluster 1
up_matrix=hits_matrix(up_hits2,:);
remove=sum(up_matrix(:,idx5==2)>0,2)<10;
up_hits2(remove)=[];

%choose 40 cells to present
sub_matrix=zeros(length(up_hits2),80);
for ii=1:length(up_hits2)
    vec1=wt_umi_corrected_z(hits(up_hits2(ii)),idx5==2);
    vec2=wt_umi_corrected_z(hits(up_hits2(ii)),idx5==1);
    [~,ind1]=sort(vec1,'descend');
    [~,ind2]=sort(vec2,'ascend');
    sub_matrix(ii,41:80)=vec1(ind1(1:40));
    sub_matrix(ii,1:40)=vec2(ind2(1:40));
end

%Fig 4B - Heatmap of genes up-regualted in highly infected cells

cc=clustergram(sub_matrix,'RowLabels',names_wt_hits(up_hits2),'Standardize','row','Cluster','column','Colormap',redbluecmap);
h1=cc.plot; hits_data=h1.Children(1).CData; gene_up=cc.RowLabels; close;

figure; imagesc(hits_data); colormap(redbluecmap); title 'gene up-regulated in highly infected cells, wt'
h=get(gca);
h.XAxis.TickLength=[0 0];
h.YAxis.TickLength=[0 0];
set(gca,'YTick',1:length(gene_up));
set(gca,'YTickLabel',gene_up);
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',10)

%Supplementary Table 1
T=table(names_wt_hits(up_hits2),hits_medians(up_hits2,1),hits_medians(up_hits2,2),hits_p_value(up_hits2),...
         'VariableNames',{'Gene','Median_expression_cluster_1','Median_expression_cluster_2','pvalue_FDR_corrected'});
writetable(T,'Sup_Table_1_wt_scRNAseq_up_in_high.xlsx');

% dICP0
%differently expressed by Wilcoxson test 
%cluster 4 vs 1
p=[];
for ii=1:size(d0_g2m_corrected_z,1)
    p(ii,1)=ranksum(d0_g2m_corrected_z(ii,idx6==4),d0_g2m_corrected_z(ii,idx6==1));
end
p(:,2)=mafdr(p(:,1),'BHFDR',true);

hits_d=find(p(:,2)<0.05);
hits_d_matrix=d0_g2m_corrected_z(hits_d,:);
names_d0_hits_d=names_d0(hits_d);
hits_d_p_value=p(hits_d,2);
hits_d_medians=[];
hits_d_medians(:,1)=median(hits_d_matrix(:,idx6==4),2);
hits_d_medians(:,2)=median(hits_d_matrix(:,idx6==1),2);
up_hits_d=find(hits_d_medians(:,1)>hits_d_medians(:,2));
down_hits_d=find(hits_d_medians(:,1)<hits_d_medians(:,2));

%keep only host genes
remove=contains(names_d0_hits_d(up_hits_d),'HSV1');
up_hits_d2=up_hits_d(~remove);

%keep only expressed in >10 cells in cluster 1
up_matrix_d=hits_d_matrix(up_hits_d2,:);
remove=sum(up_matrix_d(:,idx6==4)>0,2)<10;
up_hits_d2(remove)=[];


%choose 40 cells to present
sub_matrix_d=zeros(length(up_hits_d2),80);
for ii=1:length(up_hits_d2)
    vec1=d0_umi_corrected_z(hits_d(up_hits_d2(ii)),idx6==4);
    vec2=d0_umi_corrected_z(hits_d(up_hits_d2(ii)),idx6==1);
    [~,ind1]=sort(vec1,'descend');
    [~,ind2]=sort(vec2,'ascend');
    sub_matrix_d(ii,41:80)=vec1(ind1(1:40));
    sub_matrix_d(ii,1:40)=vec2(ind2(1:40));
end

cc=clustergram(sub_matrix_d,'RowLabels',names_d0_hits_d(up_hits_d2),'Standardize','row','Cluster','column','Colormap',redbluecmap);
h1=cc.plot; hits_d_data=h1.Children(1).CData; gene_up=cc.RowLabels; close;
figure; imagesc(hits_d_data); colormap(redbluecmap);
h=get(gca);
h.XAxis.TickLength=[0 0];
h.YAxis.TickLength=[0 0];
set(findall(gcf,'-property','FontName'),'FontName','Arial'); set(findall(gcf,'-property','FontSize'),'FontSize',10)
title 'genes up-regualted in highly infected, \DeltaICP0'

T2=table(names_d0_hits_d(up_hits_d2),hits_d_medians(up_hits_d2,1),hits_d_medians(up_hits_d2,2),hits_d_p_value(up_hits_d2),...
         'VariableNames',{'Gene','Median_expression_cluster_1','Median_expression_cluster_2','pvalue_FDR_corrected'});
writetable(T2,'Sup_Table_11_scRNAseq_d0_up_in_high.xlsx');

%% Analysis of ISGs in high vs low infected cells

%wt infection
%List of ISGs expresssed in wt-infected cells
genes_wt={'IFIT3','ISG15','OAS3','IFI6','IFI44','IFI44L',...
       'MYD88','TRIM5','TRIM14','TRIM21','STAT1'};
   
wt_g2m_corrected_z_binary=wt_g2m_corrected_z>0; %make binary to check if ISG is expressed or not

%bin to equal size, by viral expression
t1w=prctile(wt_viral,25);
t3w=prctile(wt_viral,75);

bin1w=wt_viral<t1w;
bin4w=wt_viral>=t3w;

figure; title 'ISGs expression in wt-infected cells, low vs high'
subplot(4,6,1); hold;
scatter(linspace(0.6,1.4,sum(bin1w)),log10(wt_viral(bin1w)*1e4+1),12,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
scatter(linspace(1.6,2.4,sum(bin4w)),log10(wt_viral(bin4w)*1e4+1),12,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
title ('%HSV-1 genes')
xticks(1:2); xticklabels({'low','high'});
xlim([0.4 2.6])
ylim([0,4]); yticks(0:4);  yticklabels([0,0.1,1,10,100]);
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',12)    

for ii=1:length(genes_wt)
subplot(4,6,ii+1)
hold;
bar(1,sum(wt_g2m_corrected_z_binary(find(strcmp(names_wt,genes_wt{ii})),bin1w))/sum(bin1w)*100);
bar(2,sum(wt_g2m_corrected_z_binary(find(strcmp(names_wt,genes_wt{ii})),bin4w))/sum(bin4w)*100);
title (genes_wt{ii})
xticks(1:2); xticklabels({'low','high'});
xlim([0.4 2.6])
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
end

pvals_wt=[];
for ii=1:length(genes_wt)    
    [~,pvals_wt(ii,1)]=vartest2(wt_g2m_corrected_z_binary(find(strcmp(names_wt,genes_wt{ii})),bin1w),...
                                wt_g2m_corrected_z_binary(find(strcmp(names_wt,genes_wt{ii})),bin4w));
end
pvals_wt(:,2)=mafdr(pvals_wt(:,1),'BHFDR',true);

%dICP0 infection
%List of ISGs expresssed in dICP0-infected cells
genes={'IFIT1','IFIT3','IFITM1','ISG15','ISG20','MX1','MX2','OAS1','OAS2','OAS3','OASL','IFI6','IFI44','IFI44L',...
       'IRF7','MYD88','RSAD2','TRIM5','TRIM14','TRIM21','TRIM22','TRIM69','STAT1'};

d0_g2m_corrected_z_binary=d0_g2m_corrected_z>0; %make binary to check if ISG is expressed or not

%bin to equal size, by viral expression
t1=prctile(d0_viral,25);
t3=prctile(d0_viral,75);

bin1=d0_viral<t1;
bin4=d0_viral>=t3;

figure; title 'ISGs expression in \DeltaICP0-infected cells, low vs high'
subplot(4,6,1); hold;
scatter(linspace(0.6,1.4,sum(bin1)),log10(d0_viral(bin1)*1e4+1),12,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
scatter(linspace(1.6,2.4,sum(bin4)),log10(d0_viral(bin4)*1e4+1),12,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
title ('%HSV-1 genes')
xticks(1:2); xticklabels({'low','high'});
xlim([0.4 2.6])
ylim([0,4]); yticks(0:4);  yticklabels([0,0.1,1,10,100]);
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',12)    

for ii=1:length(genes)
subplot(4,6,ii+1)
hold;
bar(1,sum(d0_g2m_corrected_z_binary(find(strcmp(names_d0,genes{ii})),bin1))/sum(bin1)*100);
bar(2,sum(d0_g2m_corrected_z_binary(find(strcmp(names_d0,genes{ii})),bin4))/sum(bin4)*100);
title (genes{ii})
xticks(1:2); xticklabels({'low','high'});
xlim([0.4 2.6])
set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
end

pvals_d0=[];
for ii=1:length(genes)    
    [~,pvals_d0(ii,1)]=vartest2(d0_g2m_corrected_z_binary(find(strcmp(names_d0,genes{ii})),bin1),...
                                d0_g2m_corrected_z_binary(find(strcmp(names_d0,genes{ii})),bin4));
end
pvals_d0(:,2)=mafdr(pvals_d0(:,1),'BHFDR',true); 
toc





