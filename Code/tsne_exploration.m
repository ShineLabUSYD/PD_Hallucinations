%% Load in data

% 400 schaefer x 1050 TRs x 59 subjects
load('D:\PhD\corts_ts2a\cort_ts1c.mat');

% Reshape into 3 dimensions
%cort_ts1c = reshape(cort_ts1c,[400 1050 60]);

% As subject 60 data is crappy, exclude
%cort_ts1c = cort_ts1c(:,:,1:59);

%% Data exploration (time series)

% Plot mean time series
mean_ts = mean(cort_ts1c,3);
figure; plot(mean_ts');
% Check out the outlier time series
ts_std = std(mean_ts,1,2);
min_std = min(ts_std);
max_std = max(ts_std);
for ii = 1:length(ts_std)
    if ts_std(ii) == min_std
        min_index = ii;
    elseif ts_std(ii) == max_std
        max_index = ii;
    end
end
% Plot regions with min and max movement
figure; plot(mean_ts(min_index,:)');
figure; plot(mean_ts(max_index,:)');

%% Zscore time series
ts_zscore = zeros(400,1050,59);
for ii = 1:size(cort_ts1c,3)
    sub_ts = cort_ts1c(:,:,ii); %raw time series
    sub_tsz = zeros(400,size(cort_ts1c,2)); %predefine
    for j = 1:400 %for each region in turn
        sub_tsz(j,:) = zscore(sub_ts(j,:)); %normalize the data so that each region has mean = 0 and standard deviation = 1
    end
    ts_zscore(:,:,ii) = sub_tsz;
end
% Plot
mean_tszscore = mean(ts_zscore,3);
plot(mean_tszscore');

%% Calculate correlation matrices

% Storage of correlation matrices
storage(size(ts_zscore,3)) = struct('ts_corr',1);
% Calculate correlation
for ii = 1:size(ts_zscore,3)
    temp = ts_zscore(:,:,ii)';
    ts_corr = corr(temp);
    storage(ii).ts_corr = ts_corr;
end

% Calculate average correlation matrix
for ii = 1:length(storage)
    if ii == 1
        sum_corr = storage(ii).ts_corr;
    else
        sum_corr = sum_corr + storage(ii).ts_corr;
    end
end
mean_corr = sum_corr/length(storage);
% Plot avg correlation matrix
figure; imagesc(mean_corr); colorbar;

%% Check networks
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\grad_schaef.mat');
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
imagesc_ids(mean_corr,grad_schaef);

% Visualise on cortical surface
mean_corr_vector = mean(mean_corr,2);
surf_schaef(schaef_id,'schaefer 400'); %template
surf_schaef(mean_corr_vector,'test data'); %data

% Visualise on Margulies gradient
avgGM = GradientMaps('kernel','na','approach','dm');
avgGM = avgGM.fit(mean_corr,'sparsity',90,'alpha',0.5);
surf_schaef(avgGM.gradients{1}(:,1),'test Gradient');
flipgradient = -(avgGM.gradients{1});
surf_schaef(flipgradient(:,1),'test Gradient');

% Using Enigma toolbox
test = parcel_to_surface(schaef_id,'schaefer_400_conte69'); %template
test = parcel_to_surface(flipgradient(:,1),'schaefer_400_conte69'); %data
f = figure; plot_cortical(test, 'surface_name', 'conte69', 'cmap', 'bb_g');

%% Calculate False nearest neighbours (FNN) across time for each region time series

fnn = zeros(10,4,400);
for ii = 1:size(mean_ts,1)
    fnn(:,:,ii) = false_nearest(mean_ts(ii,:),1,10,1,10);
end

% Plot
% All regions
figure; plot(squeeze(fnn(:,2,:)))
xlabel('number of PCs')
ylabel('% false nearest neighbours')
yline(0.1)

% Number of PCs = 3 for ideal FNN

%% t-SNE default parameters
% Perplexity = 30, Exaggeration = 4

rng(1234);
[twoDim_tsne,loss] = tsne(mean_corr,'Algorithm','barneshut','Distance','euclidean','NumDimensions',2,'NumPCAComponents',3);
rng(1234);
[threeDim_tsne,loss2] = tsne(mean_corr,'Algorithm','barneshut','Distance','euclidean','NumDimensions',3,'NumPCAComponents',3);

% Single tsne test for 1 subject
tsne_storage(10) = struct('data',0,'loss',0);
for ii = 1:10
    [tsne_out,loss] = tsne(meanCon_matrix,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);
            tsne_storage(ii).data = tsne_out;
            tsne_storage(ii).loss = loss;
end
for ii = 1:8
    subplot(2,4,ii)
    scatter3(tsne_storage(ii).data(:,1),tsne_storage(ii).data(:,2),tsne_storage(ii).data(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    colorbar;
end

            
            
% Plot
figure; gscatter(twoDim_tsne(:,1),twoDim_tsne(:,2),schaef_id);
figure; scatter3(threeDim_tsne(:,1),threeDim_tsne(:,2),threeDim_tsne(:,3),10,schaef_id,'filled');
% Loss comparison
fprintf('2-D embedding has loss %g, and 3-D embedding has loss %g.\n',loss,loss2)

% k-nearest neighbours setting k = sqrt(N) (how many points to look for)
Idx1 = knnsearch(mean_corr,mean_corr,'K',21,'IncludeTies',true,'Distance','euclidean');
Idx2 = knnsearch(twoDim_tsne,twoDim_tsne,'K',21,'IncludeTies',true,'Distance','euclidean');
% Remove 1st element (self) from each vector
for ii = 1:length(Idx1)
    Idx1{ii} = Idx1{ii}(2:20);
    Idx2{ii} = Idx2{ii}(2:20);
end

% Calculate fraction of k-nearest neighbours kept
total_knn = 20*400;
total_idx = 0;
for ii = 1:length(Idx1)
    idx = bsxfun(@eq,sort(Idx1{ii}),sort(Idx2{ii}(:)));
    same_idx = nnz(cumsum(idx,2) == cumsum(idx,1) & idx);
    total_idx = total_idx + same_idx;
end
fraction_knn = total_idx/total_knn;

% Spearman correlations between pairwise distance in high-dimensional space
% and in embedding
test = pdist2(perplexity_testing(1).data(1,:),perplexity_testing(1).data(2,:),'euclidean');


% Use loss output to compare how true to original embedding is 

%% t-SNE iteration of all combination for perplexity, learning rate, and exaggeration parameters

perplexity = [1 10 20 30 40 50 60 70 80 90 100];
learningRate = [1 100 200 300 400 500 600 700 800 900 1000];
exaggeration = [1 10 20 30 40 50 60 70 80 90 100];
for a = 1:length(perplexity)
    for b = 1:length(learningRate)
        for c = 1:length(exaggeration)
            rng(1234);
            [tsne_out,loss] = tsne(mean_corr,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',perplexity(a),'LearnRate',learningRate(b),...
                'Exaggeration',exaggeration(c));
            tsne_storage{a,b,c} = tsne_out;
            tsne_loss{a,b,c} = loss;
        end
    end
end

%% Euclidean Distance for every t-SNE, calculating within/between connections

load("C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\PD_Hallucinations\tsne_iteration.mat");
test = cell2mat(tsne_storage);
test2 = reshape(test,400,3993);
[coeff,score,latent,tsquared,explained,mu] = pca(test2);

% Calculate euclidean distance for each tsne output
for a = 1:length(tsne_storage)
    for b = 1:length(tsne_storage)
        for c = 1:length(tsne_storage)
            temp = cell2mat(tsne_storage(a,b,c)); % 400x3 tsne plot
            for d = 1:length(temp)
                for e = 1:length(temp)
                     dist(d,e) = pdist2(temp(d,:),temp(e,:),'euclidean');
                end
            end
            tsne_eucd{a,b,c} = dist; % 400x400 eucd
        end
    end
end
% Create template within and between matrices
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
for jj = 1:400
    for kk = 1:400
        if schaef_id(jj)==schaef_id(kk)
            within(jj,kk) = 1;
        else
            between(jj,kk) = 1;
        end
    end
end

% Match every eucd of t-SNE with both the within and between templates
% Calculate ratio of within vs between eucd
for a = 1:size(tsne_eucd,1)
    for b = 1:size(tsne_eucd,2)
        for c = 1:size(tsne_eucd,3)
            temp = cell2mat(tsne_eucd(a,b,c)); % 400x400 eucd
            idx_within = 1; % counter for within
            idx_between = 1; % counter for between
            for d  = 1:size(within,1)
                for e = 1:size(within,2)
                    if within(d,e) == 1
                        temp_within(idx_within) = temp(d,e); % vector of ED for within
                        idx_within = idx_within + 1;
                    elseif within(d,e) == 0
                        temp_between(idx_between) = temp(d,e); % vector of ED for between
                        idx_between = idx_between + 1;
                    end
                end
            end
            mean_within = mean(temp_within);
            mean_between = mean(temp_between);
            ratio = mean_within/mean_between;
            tsne_ratio{a,b,c} = ratio; % within vs between ED
            fprintf('a=%i b=%i c=%i\n',a,b,c)
        end
    end
end
% High ratio (>1) = ED within network > ED between network -> less distinct
% networks
% Low ratio (<1) = ED within network < ED between network -> distinct
% networks
tsne_ratio = cell2mat(tsne_ratio);

% Min and Max t-SNE ratios
min_ratio = min(tsne_ratio,[],'all');
max_ratio = max(tsne_ratio,[],'all');
% Visualisation of t-SNE ratio
for ii = 1:size(tsne_ratio,3)
    subplot(4,3,ii)
    imagesc(tsne_ratio(:,:,ii))
    xlabel('Learning Rate')
    ylabel('Perplexity')
    title(['Exaggeration = ',num2str(ii)])
    caxis([min_ratio max_ratio])
end

%% Intial Visualisations based off t-SNE ratios

% Look at distribution of ratios
figure; histogram(tsne_ratio,'Normalization','probability','BinWidth',0.01);
title('Histogram of T-SNE ratios')
xlabel('T-SNE Ratio')
ylabel('Probability')

% Look at distribution of correlations
figure; histogram(tsne_meanCorr,'Normalization','probability','BinWidth',0.01);
title('Histogram of T-SNE correlations')
xlabel('T-SNE correlation')
ylabel('Probability')

% Find index of t-SNE plots that had highest probability for lower ratios
min_ratio = min(tsne_ratio,[],'all');
max_ratio = max(tsne_ratio,[],'all');
ii = 1;
for a = 1:size(tsne_ratio,1)
    for b = 1:size(tsne_ratio,2)
        for c = 1:size(tsne_ratio,3)
            if round(tsne_ratio(a,b,c),3) == round(min_ratio,3)
                idx1{ii} = [a b c];
                ii = ii + 1;
            end
        end
    end
end

load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end
% Plot all idx1 onto one panel
for ii = 1:4
    temp = cell2mat(tsne_storage(idx1{1,ii}(1,1),idx1{1,ii}(1,2),idx1{1,ii}(1,3)));
    subplot(2,2,ii)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    colorbar
end

% Find index of t-SNE plots that had highest probability for higher ratios
ii = 1;
for a = 1:size(tsne_ratio,1)
    for b = 1:size(tsne_ratio,2)
        for c = 1:size(tsne_ratio,3)
            if round(tsne_ratio(a,b,c),2) == round(max_ratio,2)
                idx2{ii} = [a b c];
                ii = ii + 1;
            end
        end
    end
end

% Plot all idx2 onto one panel
for ii = 1:1
    temp = cell2mat(tsne_storage(idx2{1,ii}(1,1),idx2{1,ii}(1,2),idx2{1,ii}(1,3)));
    subplot(1,1,ii)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
end

%% Pearsons correlation of euclidean distances low vs high dimension

% Calculate pearsons correlation between pairwise distances of low and high
% dimensional embedding
% High dim embedding
high_distance = zeros(400,400);
for ii = 1:length(mean_corr)
    roi = mean_corr(ii,:);
    for jj = 1:length(roi)
        distance = pdist2(roi,mean_corr(jj,:),'euclidean');
        high_distance(ii,jj) = distance;
    end
end

% Pearsons correlation between low and high dimensional euclidean distance
% High correlation means more faithful embedding
for a = 1:size(tsne_eucd,1)
    for b = 1:size(tsne_eucd,2)
        for c = 1:size(tsne_eucd,3)
            low_dist = cell2mat(tsne_eucd(a,b,c));
            correlations = zeros(400,1);
            for d = 1:length(low_dist)
                correlation = corr(low_dist(:,d),high_distance(:,d));
                correlations(d) = correlation;
            end
            tsne_corr{a,b,c} = correlations;
        end
    end
end

% Calculate average correlation for each tsne iteration
for a = 1:size(tsne_corr,1)
    for b = 1:size(tsne_corr,2)
        for c = 1:size(tsne_corr,3)
            tsne_meanCorr{a,b,c} = mean(cell2mat(tsne_corr(a,b,c)));
        end
    end
end

% Calculate min and max correlations
tsne_meanCorr = cell2mat(tsne_meanCorr);
min_meanCorr = min(tsne_meanCorr,[],'all');
max_meanCorr = max(tsne_meanCorr,[],'all');
% Find tsne plots that are close to min and max
ii = 1;
jj = 1;
for a = 1:size(tsne_meanCorr,1)
    for b = 1:size(tsne_meanCorr,2)
        for c = 1:size(tsne_meanCorr,3)
            if round(tsne_meanCorr(a,b,c),2) > 0.80%round(min_meanCorr,2)
                idx_min{ii} = [a b c];
                ii = ii + 1;
            elseif round(tsne_meanCorr(a,b,c),2) == round(max_meanCorr,2)
                idx_max{jj} = [a b c];
                jj = jj + 1;
            end
        end
    end
end

% Visualise tsne plots for distance correlations
% Plot all idx1 onto one panel
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end
% Plot tsne plots from idx_max
for ii = 1:10
    temp = cell2mat(tsne_storage(idx_max{1,ii}(1,1),idx_max{1,ii}(1,2),idx_max{1,ii}(1,3)));
    subplot(2,5,ii)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
end
% Plot tsne plots from idx_min
a = 1;
for ii = 542:549
    temp = cell2mat(tsne_storage(idx_min{1,ii}(1,1),idx_min{1,ii}(1,2),idx_min{1,ii}(1,3)));
    subplot(2,4,a)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    a = a + 1;
end

%% Spearman correlation of low-D ED with static FC in high-D

% Spearman correlation between low and high dimensional euclidean distance
% High correlation means more faithful embedding
for a = 1:size(tsne_eucd,1)
    for b = 1:size(tsne_eucd,2)
        for c = 1:size(tsne_eucd,3)
            low_dist = cell2mat(tsne_eucd(a,b,c));
            correlations = zeros(400,1);
            for d = 1:length(low_dist)
                correlation = corr(low_dist(:,d),mean_corr(:,d));
                correlations(d) = correlation;
            end
            tsne_corr2{a,b,c} = correlations;
        end
    end
end

% Calculate average correlation for each tsne iteration
for a = 1:size(tsne_corr2,1)
    for b = 1:size(tsne_corr2,2)
        for c = 1:size(tsne_corr2,3)
            tsne_distfcCorr{a,b,c} = mean(cell2mat(tsne_corr2(a,b,c)));
        end
    end
end
tsne_distfcCorr = cell2mat(tsne_distfcCorr);

% Check results
figure; histogram(tsne_distfcCorr);
for ii = 1:size(tsne_distfcCorr,3)
    subplot(3,4,ii)
    imagesc(tsne_distfcCorr(:,:,ii))
    colorbar
    title(['Exaggeration = ',num2str(ii)])
end
% Compare with meanCorr results
for ii = 1:size(tsne_meanCorr,3)
    subplot(3,4,ii)
    imagesc(tsne_meanCorr(:,:,ii))
    colorbar
    title(['Exaggeration = ',num2str(ii)])
end

% Plot tsne for min and max
min_corr = min(tsne_distfcCorr,[],'all');
max_corr = max(tsne_distfcCorr,[],'all');
% Find tsne plots that are close to min and max
ii = 1;
jj = 1;
for a = 1:size(tsne_distfcCorr,1)
    for b = 1:size(tsne_distfcCorr,2)
        for c = 1:size(tsne_distfcCorr,3)
            if round(tsne_distfcCorr(a,b,c),2) < -0.67%round(min_meanCorr,2)
                idx_min{ii} = [a b c];
                ii = ii + 1;
            elseif round(tsne_distfcCorr(a,b,c),2) > 0.01
                idx_max{jj} = [a b c];
                jj = jj + 1;
            end
        end
    end
end

% Visualise tsne plots for distance correlations
% Plot all idx1 onto one panel
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end
% Plot tsne plots from idx_max
a = 1;
for ii = 1:3
    temp = cell2mat(tsne_storage(idx_max{1,ii}(1,1),idx_max{1,ii}(1,2),idx_max{1,ii}(1,3)));
    subplot(1,3,a)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    a = a + 1;
end
% Plot tsne plots from idx_min
a = 1;
for ii = 1:8
    temp = cell2mat(tsne_storage(idx_min{1,ii}(1,1),idx_min{1,ii}(1,2),idx_min{1,ii}(1,3)));
    subplot(2,4,a)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    a = a + 1;
end


%% k-means clustering of t-SNE results
% m x n matrix, m = observations n = variables
% 11 x 121, every 11 = 1 index of exaggeration

% Formatting of tsne ratio 
tsne_ratio2 = tsne_ratio(2:11,:,:);
tsne_ratio_flat = reshape(tsne_ratio2,11^2*10,1);

% Formatting of tsne distance correlations
tsne_meanCorr2 = tsne_meanCorr(2:11,:,:);
tsne_meanCorr_flat = reshape(tsne_meanCorr2,11^2*10,1);

% Run for range of k values, with many iterations
tsne_kmeans = zeros(1210,100,19);
a = 2;
for ii = 1:19 
    for jj = 1:100
        %[idx,C] = kmeans(tsne_ratio_flat,a,'MaxIter',100);
        [idx,C] = kmeans(tsne_meanCorr_flat,a,'MaxIter',100);
        tsne_kmeans(:,jj,ii) = idx;
    end
    a = a + 1;
end

% corr values of each cluster 
test5 = reshape(tsne_kmeans(:,1,4),10,11,11); % Adjust for specific k value
cluster1 = tsne_meanCorr2(test5==1);
cluster2 = tsne_meanCorr2(test5==2);
cluster3 = tsne_meanCorr2(test5==3);
cluster4 = tsne_meanCorr2(test5==4);
cluster5 = tsne_meanCorr2(test5==5);
mean_cluster1 = mean(cluster1);
mean_cluster2 = mean(cluster2);
mean_cluster3 = mean(cluster3);
mean_cluster4 = mean(cluster4);
mean_cluster5 = mean(cluster5);

tsne_storage = tsne_storage(2:11,:,:);
cluster1_tsne = tsne_storage(test5==2);

% Plot results
% Visualise tsne plots for distance correlations
% Plot all idx1 onto one panel
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end
% Plot tsne plots from tsne cluster
a = 1;
for ii = 1:4
    temp = cluster1_tsne{ii,1};
    subplot(2,2,a)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    colorbar;
    grid off
    a = a + 1;
end

% Check within/between ratios of cluster 3
cluster1_ratio = tsne_ratio(test5==2);
mean_ratio = mean(cluster1_ratio);


%% AMI - Adjusted mutual information
% Run AMI for each iteration against every other iteration within a k value
tsne_ami = zeros(100,100,19);
for ii = 1:size(tsne_kmeans,3)
    sim_all = zeros(100,100);
    for jj = 1:size(tsne_kmeans,2)
       for kk = 1:size(tsne_kmeans,2)
           temp = ami(tsne_kmeans(:,jj,ii),tsne_kmeans(:,kk,ii));
           sim_all(kk,jj) = temp;
       end  
    end
    tsne_ami(:,:,ii) = sim_all;
end

% Plot mean and SD
tsne_ami2 = reshape(tsne_ami,100^2,19);
figure; plot(mean(tsne_ami2));
figure; shadedErrorBar(2:1:20,mean(tsne_ami2),std(tsne_ami2,[],1));
figure; shadedErrorBar(2:1:20,mean(tsne_ami2),std(tsne_ami2,[],1)/100);

% Looking at clustering k results
test5 = reshape(tsne_kmeans(:,1,4),10,11,11); % Adjust for specific k value
test8 = reshape(tsne_kmeans(:,1,17),10,11,11); % Adjust for specific k value
for ii = 1:11
    subplot(4,3,ii)
    imagesc(test5(:,:,ii))
    caxis([1 5]) % Adjust for no. of clusters
    colormap('jet')
    colorbar
end

% Prepping data to run statistic tests between clusters
cluster1_ratio = tsne_ratio2(test3==1);
cluster2_ratio = tsne_ratio2(test3==2);
cluster3_ratio = tsne_ratio2(test3==3);
cluster4_ratio = tsne_ratio2(test3==4);

cluster_all3(3) = struct('Cluster',0);
cluster_all3(1).cluster = cluster1_ratio;
cluster_all3(2).cluster = cluster2_ratio;
cluster_all3(3).cluster = cluster3_ratio;

cluster1_ratio = tsne_ratio2(test8==1);
cluster2_ratio = tsne_ratio2(test8==2);
cluster3_ratio = tsne_ratio2(test8==3);
cluster4_ratio = tsne_ratio2(test8==4);
cluster5_ratio = tsne_ratio2(test8==5);
cluster6_ratio = tsne_ratio2(test8==6);
cluster7_ratio = tsne_ratio2(test8==7);
cluster8_ratio = tsne_ratio2(test8==8);

mean_cluster(1) = mean(cluster1_ratio);
mean_cluster(2) = mean(cluster2_ratio);
mean_cluster(3) = mean(cluster3_ratio);
mean_cluster(4) = mean(cluster4_ratio);
mean_cluster(5) = mean(cluster5_ratio);
mean_cluster(6) = mean(cluster6_ratio);
mean_cluster(7) = mean(cluster7_ratio);
mean_cluster(8) = mean(cluster8_ratio);
figure; plot(mean_cluster);
[h,p] = ttest2(cluster1_ratio,cluster4_ratio);

cluster_all8(8) = struct('Cluster',0);
cluster_all8(1).Cluster = cluster1_ratio;
cluster_all8(2).Cluster = cluster2_ratio;
cluster_all8(3).Cluster = cluster3_ratio;
cluster_all8(4).Cluster = cluster4_ratio;
cluster_all8(5).Cluster = cluster5_ratio;
cluster_all8(6).Cluster = cluster6_ratio;
cluster_all8(7).Cluster = cluster7_ratio;
cluster_all8(8).Cluster = cluster8_ratio;
[sig,pval] = perm_1d_delta(cluster1_ratio,cluster4_ratio,5000);

for ii = 1:length(cluster_all)
    for jj = 1:length(cluster_all)
        [sig,pval] = perm_1d_delta(cluster_all(ii).Cluster,cluster_all(jj).Cluster,5000);
        pvalues{jj} = [sig,pval];
    end
    cluster_all(ii).ptest = pvalues;
end

map = brewermap(8,'Set1');
figure; histogram(cluster1_ratio);
hold on
histogram(cluster2_ratio);
histogram(cluster3_ratio);
hold off


map = brewermap(8,'Set1');
figure; histogram(cluster1_ratio);
hold on
histogram(cluster2_ratio);
histogram(cluster3_ratio);
histogram(cluster4_ratio);
histogram(cluster5_ratio);
histogram(cluster6_ratio);
histogram(cluster7_ratio);
histogram(cluster8_ratio);
hold off

% AMI between top 2 k values
idx3 = tsne_kmeans(:,1,2); % Change to whichever k value
idx8 = tsne_kmeans(:,1,7); % Change to whichever k value
for kk = 1:3
    for jj = 1:8
        temp1 = double(idx3==kk);
        temp2 = double(idx8==jj);
        temp1(temp1==0) = 2;
        temp2(temp2==0) = 2;
        ami_3_8(kk,jj) = ami(temp1,temp2);
    end
end

% Visualise ami between clusters
figure;imagesc(ami_3_8);colorbar;
xlabel('k=8')
ylabel('k=3')
figure;plot(ami_5_20');



%% Bohm et al. Method

% Check distance correlations for lowest ratios
low_ratio = [tsne_meanCorr(10,1,5) tsne_meanCorr(10,1,6) tsne_meanCorr(10,1,7) tsne_meanCorr(11,1,8)];

% Find tsne plots that are within certain threshold
ii = 1;
for a = 1:size(tsne_meanCorr,1)
    for b = 1:size(tsne_meanCorr,2)
        for c = 1:size(tsne_meanCorr,3)
            if round(tsne_meanCorr(a,b,c),2) >= 0.83 && round(tsne_ratio(a,b,c),2) <= 0.54
                sweetspot{ii} = [a,b,c];
                ii = ii + 1;
            end
        end
    end
end

% Plot results
% Visualise tsne plots for distance correlations
% Plot all idx1 onto one panel
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
% 3-Dimensional plots
hexMap = [0 0.4470 0.7410 0.8500 0.3250 0.0980 0.9290 0.6940 0.1250 0.4940 0.1840 0.5560 0.4660 0.6740 0.1880 0.3010 0.7450 0.9330 0.6350 0.0780 0.1840];
myColorMap = zeros(length(hexMap)/3, 3); % Preallocate
a = 1;
for k = 1 : length(myColorMap)
	r = hexMap(a);
	g = hexMap(a+1);
	b = hexMap(a+2);
	myColorMap(k, :) = [r, g, b];
    a = a + 3;
end
% Plot tsne plots from idx_max
a = 1;
for ii = 1:6
    temp = cell2mat(tsne_storage(sweetspot{1,ii}(1,1),sweetspot{1,ii}(1,2),sweetspot{1,ii}(1,3)));
    subplot(2,3,a)
    scatter3(temp(:,1),temp(:,2),temp(:,3),10,schaef_id,'filled');
    colormap(myColorMap);
    colorbar;
    a = a + 1;
end

% Compare distance correlation for 'stingray' plots vs other plots - Bohm
% et al
% Compare distance correlation between ideal cluster and other clusters -
% Bohm et al
% Get euclidean distances for tsne plots
tsne_eucd = tsne_eucd(2:11,:,:);
cluster1_eucd = tsne_eucd(test5==1);
clusterother_eucd = tsne_eucd(test5~=1);
% Calculate correlation between tsne plots from sweetspot
% Randomly sample 5 plots and calculate correlation
within_corr(20) = struct('data',0);
test_within = datasample(cluster1_eucd,20,'Replace',false);
for ii = 1:length(test_within)
    temp_group = zeros(400,20);
    for a = 1:length(test_within{ii,1})
        for b = 1:length(test_within)
            temp_corr = corr(test_within{ii,1}(:,a),test_within{b,1}(:,a));
            temp_group(a,b) = temp_corr;
        end
    end
    within_corr(ii).data = temp_group;
end
% Average correlation for each within plot
for ii = 1:length(within_corr)
    within_corr(ii).mean = mean(within_corr(ii).data,[1 2]);
end

% Calculate correlation of each plot with every other plot not in sweetspot
% Randomly sample 5 plots and calculate correlation of within to these
test_between = datasample(clusterother_eucd,20,'Replace',false);
% Calculate correlation of each within tsne with between tsne
between_corr(20) = struct('data',0,'mean',0);
for ii = 1:length(test_within)
    temp_group = zeros(400,20);
    for a = 1:length(test_within{ii,1})
        for b = 1:length(test_between)
            temp_corr = corr(test_within{ii,1}(:,a),test_between{b,1}(:,a));
            temp_group(a,b) = temp_corr;
        end
    end
    between_corr(ii).data = temp_group;
end

% Average correlation for each between plot
for ii = 1:length(between_corr)
    between_corr(ii).mean = mean(between_corr(ii).data,[1 2]);
end

% Visualise results
for ii = 1:length(between_corr)
    subplot(4,5,ii)
    imagesc(between_corr(ii).data)
    colorbar
end

% Average
corr_all = [vertcat(within_corr.mean), vertcat(between_corr.mean)];
g = {'within','between'};
figure; boxplot(corr_all,g);
ylabel('Mean correlation')

[sig,pval] = perm_1d_delta(corr_all(:,1),corr_all(:,2),5000);
[sig,pval] = permutationTest(corr_all(:,1),corr_all(:,2),5000);