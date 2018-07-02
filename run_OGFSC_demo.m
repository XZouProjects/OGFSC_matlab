clear;
close all;
warning off;

%% load data and pre-processing
[~, ~, raw] = xlsread('GSE59739_DataTable.xlsx', 1); 
data = cell2mat(raw(3:end, 2:end));
pickingSession = raw(1, 2:end);
cell_IDs = raw(2, 2:end);
geneList = raw(3:end, 1);
cell_type_unique = {'NF1','NF2','NF3', 'NF4','NF5','NP1', 'NP2', 'NP3', 'PEP1', 'PEP2', 'TH'};
cell_types = zeros(size(cell_IDs));
for i = 1:length(cell_type_unique)
    idx1 = strcmp(cell_IDs, cell_type_unique{i}); 
    idx2 = find(idx1==1); % idx2==1, matched cell type
    cell_types(idx2) = i;
end
idx = find(cell_types==0);
cell_types(idx) = [];
data(:,idx) = [];

% % removing batch effect
pickingSession(idx) = [];
pickingSession_uni = unique(pickingSession);
pickingSe = zeros(size(pickingSession));
for i = 1:length(pickingSession_uni)
    idx = find(strcmp(pickingSession_uni{i}, pickingSession));
    pickingSe(idx) = i;
end
P = ones(size(data,1),1);
for i = 1:length(P)
    P(i) = anova1(data(i,:), pickingSe,'off');
end
[~, ~, adj_p]=fdr_bh(P);
idx = find(adj_p<0.05);
data(idx,:) = [];
geneList(idx) = [];
sigma2 = var(data,0,2);
%% log transform
log2Data = log2(data+1);
%% gene filtering by OGFSC and save result
[OGFSC_idx, idx_output, cv2_threshold] = OGFSC(log2Data, 'plot_option', 1);
dataBuffer_OGFSC = log2Data(OGFSC_idx,:);
temp = num2cell(data(OGFSC_idx,:));
temp = [geneList(OGFSC_idx), temp];
delete('OGFSCFilteredData.xlsx');
xlswrite('OGFSCFilteredData.xlsx', temp);
%% scatter plots
X = dataBuffer_OGFSC';
[~,SCORE] = pca(X);
ydata = tsne(X,[],SCORE(:,1:5), 30);
figure;
colormap(jet);
Color=get(gcf,'Colormap');
stepsize = floor(64/max(cell_types));
buffer = cell(1, max(cell_types));
for i = 1:max(cell_types)
    idx = find(cell_types==i);
    plot(ydata(idx,1), ydata(idx,2), 'o', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', Color(64-(i-1)*stepsize,:), 'MarkerSize', 3);
    if i == 1
        hold on;
    end
end
xlabel('tSNE_1');
ylabel('tSNE_2');
title('OGFSC');
legend(cell_type_unique);
set(gca, 'fontsize', 14);

%%
N1 = round(0.06*size(data,2));
N2 = round(0.94*size(data,2));
idx1 = find(sum(data>2,2)<N1);
idx2 = find(sum(data>0,2)>N2);
Idx_ref = unique([idx1;idx2]);

% % scatter plot
X = log2Data(Idx_ref,:)';
[coeff_Ref,SCORE] = pca(X);
ydata = tsne(X,[],SCORE(:,1:5), 30);
figure;
% colormap(jet);
% Color=get(gcf,'Colormap');
stepsize = floor(64/max(cell_types));
buffer = cell(1, max(cell_types));
for i = 1:max(cell_types)
    idx = find(cell_types==i);
    plot(ydata(idx,1), ydata(idx,2), 'o', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', Color(64-(i-1)*stepsize,:), 'MarkerSize', 4);
    if i == 1
        hold on;
    end
end
xlabel('tSNE_1');
ylabel('tSNE_2');
legend(cell_type_unique);
set(gca, 'fontsize', 14);
