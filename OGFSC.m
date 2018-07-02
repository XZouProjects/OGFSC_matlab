function [OGFSC_idx, idx_output, cv2_threshold] = OGFSC(data, varargin)

% OGFSC to perform optimized gene filtering for single-cell RNA-seq data
% Xin Zou & Jie Hao
% 2018, SJTU, China

warning off;
%% parameters setting up
varList = {'nBins', 'minBinSize', 'LR_p', 'LR_thresholds', 'TW_threshold', 'plot_option'};
convgThrh = {60, 100, 0.01, [0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999], 0.0001, 1}; % default values
if ~isempty(varargin)
    L = length(varargin);
    for i = 1:2:L-1
        idx = find(strcmp(varargin{i}, varList));
        convgThrh{idx} = varargin{i+1};
    end
end

nBins = convgThrh{1}; % number of bins
minBinSize = convgThrh{2}; % minimeansDatam bin size
LR_p = convgThrh{3}; % p-value threshold to identify valid linear regression model
alpha = convgThrh{4}; % list of different threshold, i.e., the upper boundary of confidence interval
TW_threshold = convgThrh{5}; % the threshold to determine the number of eigenvalues on TW distribution
plot_option = convgThrh{6}; % plot the outputs of ODFSC, by default 0. 

%% remove non-sense or dropout values
data_ori = data;
meansData = mean(data, 2);
varsData = var(data, [], 2);
cv2Data = varsData./meansData.^2;
idx = [find(isnan(cv2Data)); find(cv2Data==0); find(meansData==0)];
II = 1:size(data,1);
cv2Data(idx) = [];
meansData(idx) = [];
varsData(idx) = [];
data(idx,:) = [];
II(idx) = [];
%% Binning genes and construct regression model using MLM method
log10mean = log10(meansData);
minV = min(log10mean);
maxV = max(log10mean);
stepsize = (maxV-minV)/nBins;
boundaries = minV:stepsize:maxV;
numBins = length(boundaries)-2;
Bins = zeros(numBins,2);
BinSize = zeros(numBins,1);
BinIdx = cell(numBins,1);
for i = 1:numBins
    Bins(i,1) = boundaries(i);
    Bins(i,2) = boundaries(i+2);
    BinIdx{i} = find(log10mean>=Bins(i,1) & log10mean<Bins(i,2));
    BinSize(i) = length(BinIdx{i});
end
IdxValidBin = find(BinSize>=minBinSize);

B = zeros(2,length(IdxValidBin));
P = ones(1,length(IdxValidBin));
for i = 1:length(IdxValidBin)
    idx_temp = BinIdx{IdxValidBin(i)};
    mean_temp = meansData(idx_temp);
    cv2_temp = cv2Data(idx_temp);
    X = [ones(size(mean_temp)), 1./mean_temp];
    [B(:,i),~, ~, ~, stats] = regress(cv2_temp,X);
    P(i) = stats(3);
end

idx3 = find(P>LR_p);
M = length(P)/2;
III = max(find(idx3<M));
idx4 = idx3(III);
III = min(find(idx3>M));
idx5 = idx3(III);
if ~isempty(idx4) && ~isempty(idx5)
    IdxValidBin([1:idx4,idx5:end]) = [];
    B(:,[1:idx4,idx5:end]) = [];
elseif ~isempty(idx4) && isempty(idx5)
    IdxValidBin(1:idx4) = [];
    B(:,1:idx4) = [];
elseif isempty(idx4) && ~isempty(idx5)
    IdxValidBin(idx5:end) = [];
    B(:,idx5:end) = [];
end

idx = find([IdxValidBin;max(IdxValidBin)]-[0;IdxValidBin]>1);
if ~isempty(idx)
    idx_temp = 1:min(idx)-1;
    IdxValidBin(idx_temp) = [];
    B(:,idx_temp) = [];
end
lowBoundary = Bins(min(IdxValidBin), 1);
Bins = Bins(IdxValidBin,:);
BinIdx = BinIdx(IdxValidBin);
BinSize = BinSize(IdxValidBin);
idx4 = [];
for i = 1:size(Bins, 1)
    idx4 = [idx4; BinIdx{i}];
end
idx4 = unique(idx4);% the list of genes have valid expression levels
BinAssignment = zeros(size(Bins, 1), length(idx4)); % assign each gene to one or two bins, as bins have overlap
for i = 1:size(Bins, 1)
    [~,Locb] = ismember(BinIdx{i}, idx4);
    BinAssignment(i,Locb) = 1; % 1 means assignment to the bin, 0 means no assignment
end
cv2_hat = zeros(1,length(idx4));
for i = 1:length(idx4) % estimate average noise cv2Data for each gene
    idx5 = find(BinAssignment(:,i)==1);
    temp = 0;
    for j = 1:length(idx5) % for each LR model
        temp = temp+B(1,idx5(j))+B(2,idx5(j))/meansData(idx4(i));
    end
    cv2_hat(i) = temp/j;
end

idx_output = II(idx4);
cv2_output = cv2Data(idx4)';

%% plot regression curve if plot_option = 1
if plot_option == 1
    figure;
    plot(log10(meansData), log10(cv2Data), '.', 'MarkerEdgeColor',[0.8, 0.8, 0.8], 'markerfaceColor', [0.8, 0.8, 0.8], 'markersize', 3);
    hold on
    plot(log10(meansData(idx4)), log10(cv2_hat), 'r.','markersize', 3);
    xtick = floor(min(log10(meansData))):1:ceil(max(log10(meansData)));
    xticklabel = num2cell(10.^xtick);
    xticklabel = cellfun(@num2str, xticklabel, 'UniformOutput', false);
    set(gca,'XTick',xtick);
    set(gca,'XTickLabel',xticklabel);
    xlabel('\mu');
    ytick = floor(min(log10(cv2Data))):1:ceil(max(log10(cv2Data)));
    yticklabel = num2cell(10.^ytick);
    yticklabel = cellfun(@num2str, yticklabel, 'UniformOutput', false);
    set(gca,'YTick',ytick);
    set(gca,'YTickLabel',yticklabel);
    ylabel('CV^2');
    set(gca, 'fontsize', 14);
end

%% Constructing candiation thresholding curves and evaluting on T-W distribution 
data = data_ori(idx_output,:);
df = size(data,2)-1;
N = zeros(size(alpha));
selectedGenesIdx = cell(size(alpha));
for i = 1:length(alpha)
    idx = find(cv2_output>cv2_hat*chi2inv(alpha(i),df)/df);
    selectedData = data(idx,:)';
    [~,~,latent] = pca(selectedData);
    [p,n] = size(selectedData);
    [~, s, mu_np, sigma_np] = TW_trace_ratio_threshold(p,n,1,TW_threshold);
    N(i) = length(find((latent/(sum(latent)/p)-mu_np)/sigma_np>s));
    selectedGenesIdx{i} = idx_output(idx);
end
idx = find(N==max(N));
OGFSC_idx = selectedGenesIdx{idx(end)};
cv2_threshold = cv2_hat*chi2inv(alpha(idx(end)),df)/df;
%% plot OGFSC metrics if plot_option = 1
if plot_option == 1
    LN = N/max(N);
    figure
    plot(1:length(alpha), LN,'s-b','markersize', 10, 'markerfacecolor', 'b');
    xtick = 1:length(alpha);
    xticklabel = num2cell(alpha);
    xticklabel = cellfun(@num2str, xticklabel, 'UniformOutput', false);
    set(gca,'XTick',xtick);
    set(gca,'XTickLabel',xticklabel);
    xlabel('\alpha');
    ylabel('Normalized No. \lambda');
    set(gca, 'fontsize', 10);
end
