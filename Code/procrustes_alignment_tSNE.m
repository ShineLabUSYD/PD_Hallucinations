%% Plotting
% Plotting
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

%% Setting up data folder pathway and checking contents

% Specify the folder where the files lives.
% Patient + Controls
myFolder = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity';
% Only patient data
%myFolder = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity\Patient';
% Check to make sure that folder actually exists. Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
        % User clicked Cancel
        return;
    end
end

%Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*ts_corr_schaef*.mat');
theFiles = dir(filePattern); % Storing all filenames

%Create object to hold all the outputs
storage(length(theFiles)) = struct('name',1,'data',1);
%Loop to load in the data, store name and values of zscore
for ii = 1:length(theFiles)
    baseFileName = theFiles(ii).name; 
    fullFileName = fullfile(theFiles(ii).folder, baseFileName); %making absolute path to file
    ts_corr_schaef = load(fullFileName); %loads in data
    subjectname = extractBefore(baseFileName,'-ts_corr_schaef');
    subjectname = strrep(subjectname,'-','_');
    storage(ii).name = subjectname; %stores sub_id under .name
    storage(ii).data = ts_corr_schaef.ts_corr; %stores data under .data
end

%% Setting up VH IDs

% Folder pathway
myFolderVH = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity\Patient\VH';
if ~isfolder(myFolderVH)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderVH);
    uiwait(warndlg(errorMessage));
    myFolderVH = uigetdir(); % Ask for a new one.
    if myFolderVH == 0
        % User clicked Cancel
        return;
    end
end

% List of VH 
% Get a list of all files in the folder with the desired file name pattern.
filePatternVH = fullfile(myFolderVH, '*ts_corr_schaef*.mat'); % Looking for the 400x400 matrices
theFilesVH = dir(filePatternVH); % Storing all filenames

% Extracting subject IDs and converting to appropriate format
for a = 1:length(theFilesVH)
    fullsubjectVH = theFilesVH(a).name;
    subjectnameVH = extractBefore(fullsubjectVH,'-ts_corr');
    subjectnameVH = strrep(subjectnameVH,'-','_');
    subjectsVH(a) = {subjectnameVH};
end

% Matching VH subject IDs and storing index 
% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsVH)
    for c = 1:length(storage)
        if strcmp(subjectsVH(b),storage(c).name)
            storage(c).Index = 'VH';
        end
    end
end

%% NonVH IDs

%Setting up Non-VH folder pathway
myFolderNVH = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity\Patient\nonVH';
if ~isfolder(myFolderNVH)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderNVH);
    uiwait(warndlg(errorMessage));
    myFolderNVH = uigetdir(); % Ask for a new one.
    if myFolderNVH == 0
        % User clicked Cancel
        return;
    end
end

%Making list of subject IDs with Non-VH Status
% Get a list of all files in the folder with the desired file name pattern.
filePatternNVH = fullfile(myFolderNVH, '*ts_corr_schaef*.mat'); % Looking for the 400x400 matrices
theFilesNVH = dir(filePatternNVH); % Storing all filenames

% Extracting subject IDs and converting to appropriate format
for a = 1:length(theFilesNVH)
    fullsubjectNVH = theFilesNVH(a).name;
    subjectnameNVH = extractBefore(fullsubjectNVH,'-ts_corr');
    subjectnameNVH = strrep(subjectnameNVH,'-','_');
    subjectsNVH(a) = {subjectnameNVH};
end

% Matching Non-VH subject IDs and storing index 
% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsNVH)
    for c = 1:length(storage)
        if strcmp(subjectsNVH(b),storage(c).name)
            storage(c).Index = 'NonVH';
        end
    end
end

%% Controls ID

% Setting up Controls folder pathway
myFolderCon = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity\Controls';
if ~isfolder(myFolderCon)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderNVH);
    uiwait(warndlg(errorMessage));
    myFolderCon = uigetdir(); % Ask for a new one.
    if myFolderCon == 0
        % User clicked Cancel
        return;
    end
end

% Making list of subject IDs with Control Status
% Get a list of all files in the folder with the desired file name pattern.
filePatternCon = fullfile(myFolderCon, '*ts_corr_schaef*.mat'); % Looking for the 400x400 matrices
theFilesCon = dir(filePatternCon); % Storing all filenames

% Extracting subject IDs and converting to appropriate format
for a = 1:length(theFilesCon)
    fullsubjectCon = theFilesCon(a).name;
    subjectnameCon = extractBefore(fullsubjectCon,'-ts_corr');
    subjectnameCon = strrep(subjectnameCon,'-','_');
    subjectsCon(a) = {subjectnameCon};
end

% Matching Control subject IDs and storing index 
% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsCon)
    for c = 1:length(storage)
        if strcmp(subjectsCon(b),storage(c).name)
           storage(c).Index = 'Control';
        end
    end
end

%% Run t-SNE for all subjects

for ii = 1:length(storage)
    rng(1234);
    [tsne_out,loss] = tsne(storage(ii).data,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);
    storage(ii).tsne = tsne_out;
    storage(ii).tsne_loss = loss;
end

%% Make population average tSNE

% All subjects
all_corr = vertcat(storage.data);
all_corr = reshape(all_corr,400,400,96);
mean_corr = mean(all_corr,3);
rng(1234);
[tsne_mean,loss] = tsne(mean_corr,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);

% Visualise tsne_mean
figure; scatter3(tsne_mean(:,1),tsne_mean(:,2),tsne_mean(:,3),10,schaef_id,'filled');colormap(myColorMap);colorbar;


% Only Patients
patient_corr = zeros(400,400);
a = 0;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,"VH")
        patient_corr = patient_corr + storage(ii).data;
        a = a + 1;
    elseif isequal(storage(ii).Index,"NonVH")
        patient_corr = patient_corr + storage(ii).data;
        a = a + 1;
    end
end
mean_patient_corr = patient_corr/a;
rng(1234);
[tsne_patient_mean,loss] = tsne(mean_patient_corr,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);

% Visualise tsne_mean
figure; scatter3(tsne_patient_mean(:,1),tsne_patient_mean(:,2),tsne_patient_mean(:,3),10,schaef_id,'filled');colormap(myColorMap);colorbar;
           


%% Align each tsne to the population average
% Using procrustes() from MATLAB
% All subjects
for ii = 1:length(storage)
    [d,align,transform] = procrustes(tsne_mean,storage(ii).tsne);
    storage(ii).tsneAlign = align;
    storage(ii).dist = d;
    storage(ii).transform = transform;
end

% Patients only
for ii = 1:length(storage)
    if isequal(storage(ii).Index,"VH")
        [d,align,transform] = procrustes(tsne_patient_mean,storage(ii).tsne);
        storage(ii).tsneAlign = align;
        storage(ii).dist = d;
        storage(ii).transform = transform;
    elseif isequal(storage(ii).Index,"NonVH")
        [d,align,transform] = procrustes(tsne_patient_mean,storage(ii).tsne);
        storage(ii).tsneAlign = align;
        storage(ii).dist = d;
        storage(ii).transform = transform;
    end
end

%% Prepare tSNE for PCA + LDA analysis

for ii = 1:length(storage)
    storage(ii).tsneAlign = reshape(storage(ii).tsneAlign,1200,1);
end

% All subjects
a = 1;
b = 1;
c = 1;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,'VH')
        group_VH(:,a) = storage(ii).tsneAlign;
        a = a + 1;
    elseif isequal(storage(ii).Index,'NonVH')
        group_NonVH(:,b) = storage(ii).tsneAlign;
        b = b + 1;
    elseif isequal(storage(ii).Index,'Control')
        group_Control(:,c) = storage(ii).tsneAlign;
        c = c + 1;
    end
end
group_tsne = [group_Control,group_NonVH,group_VH];

% Patients Only
a = 1;
b = 1;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,'VH')
        group_VH(:,a) = storage(ii).tsneAlign;
        a = a + 1;
    elseif isequal(storage(ii).Index,'NonVH')
        group_NonVH(:,b) = storage(ii).tsneAlign;
        b = b + 1;
    end
end
group_tsne_patients = [group_NonVH,group_VH];


%% Reshape tSNE back to 400x3

for ii = 1:length(storage)
    storage(ii).tsneAlign = reshape(storage(ii).tsneAlign,400,3);
end

%% Visualise original tSNE plots

% VH 
figure; scatter3(storage(1).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,'filled'); %colormap(myColorMap);
hold on
scatter3(storage(2).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,'filled'); %colormap(myColorMap);
scatter3(storage(5).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,'filled'); %colormap(myColorMap);

% Control
figure; scatter3(storage(1).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,schaef_id,'filled'); 
hold on
scatter3(storage(5).tsne(:,1),storage(5).tsne(:,2),storage(5).tsne(:,3),10,schaef_id,'filled'); 
scatter3(storage(10).tsne(:,1),storage(10).tsne(:,2),storage(10).tsne(:,3),10,schaef_id,'filled'); 

% Across groups
figure; scatter3(storage(1).tsneAlign(:,1),storage(1).tsneAlign(:,2),storage(1).tsneAlign(:,3),10,'filled'); 
hold on
scatter3(storage(20).tsneAlign(:,1),storage(20).tsneAlign(:,2),storage(20).tsneAlign(:,3),10,'filled'); 
scatter3(storage(50).tsneAlign(:,1),storage(50).tsneAlign(:,2),storage(50).tsneAlign(:,3),10,'filled');
scatter3(storage(80).tsneAlign(:,1),storage(80).tsneAlign(:,2),storage(80).tsneAlign(:,3),10,'filled');

figure; scatter3(storage(1).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,'filled'); 
hold on
scatter3(storage(20).tsne(:,1),storage(20).tsne(:,2),storage(20).tsne(:,3),10,'filled'); 
scatter3(storage(50).tsne(:,1),storage(50).tsne(:,2),storage(50).tsne(:,3),10,'filled');
scatter3(storage(80).tsne(:,1),storage(80).tsne(:,2),storage(80).tsne(:,3),10,'filled');

figure; scatter3(storage(14).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,schaef_id,'filled'); colormap(myColorMap);
hold on
scatter3(storage(1).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,schaef_id,'filled'); colormap(myColorMap);
scatter3(storage(3).tsne(:,1),storage(1).tsne(:,2),storage(1).tsne(:,3),10,schaef_id,'filled'); colormap(myColorMap);

%% Visualise aligned tSNE plots
load('D:\PhD\PD_Hallucinations\group_tsne_all.mat');

% Recreate 3D matrix of tsne plots
tsne_3D = zeros(400,3,96);
for ii = 1:size(group_tsne,2)
    temp = group_tsne(:,ii);
    temp_3D = reshape(temp,400,3);
    tsne_3D(:,:,ii) = temp_3D;
end

figure; scatter3(tsne_3D(:,1,1),tsne_3D(:,2,1),tsne_3D(:,3,1),10,'filled');
hold on
scatter3(tsne_3D(:,1,20),tsne_3D(:,2,20),tsne_3D(:,3,20),10,'filled'); 
scatter3(tsne_3D(:,1,50),tsne_3D(:,2,50),tsne_3D(:,3,50),10,'filled'); 
scatter3(tsne_3D(:,1,80),tsne_3D(:,2,80),tsne_3D(:,3,80),10,'filled');


figure; scatter3(tsne_3D(:,1,1),tsne_3D(:,2,1),tsne_3D(:,3,1),10,schaef_id,'filled'); colormap(myColorMap);
hold on
scatter3(tsne_3D(:,1,20),tsne_3D(:,2,20),tsne_3D(:,3,20),10,schaef_id,'filled'); colormap(myColorMap);
scatter3(tsne_3D(:,1,50),tsne_3D(:,2,50),tsne_3D(:,3,50),10,schaef_id,'filled'); colormap(myColorMap);
scatter3(tsne_3D(:,1,80),tsne_3D(:,2,80),tsne_3D(:,3,80),10,schaef_id,'filled'); colormap(myColorMap);

%% Group Visualisation

% Load data
load('D:\PhD\PD_Hallucinations\staticFC_meanControl');
load('D:\PhD\PD_Hallucinations\staticFC_meanNonVH');
load('D:\PhD\PD_Hallucinations\staticFC_meanVH');

% NonVH t-SNE
rng(1234);
[tsne_nvh,loss] = tsne(meanNVH_matrix,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);

% VH t-SNE
rng(1234);
[tsne_vh,loss] = tsne(meanVH_matrix,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);

% Control t-SNE
rng(1234);
[tsne_con,loss] = tsne(meanCon_matrix,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);

% Visualise unaligned t-SNE plots
figure; 
scatter3(tsne_vh(:,1),tsne_vh(:,2),tsne_vh(:,3),10,'filled');
hold on
scatter3(tsne_nvh(:,1),tsne_nvh(:,2),tsne_nvh(:,3),10,'filled');
scatter3(tsne_con(:,1),tsne_con(:,2),tsne_con(:,3),10,'filled');

% Procrustes alignment of groups to population
[d,align_con,transform] = procrustes(tsne_mean,tsne_con);
[d,align_nvh,transform] = procrustes(tsne_mean,tsne_nvh);
[d,align_vh,transform] = procrustes(tsne_mean,tsne_vh);

% Visualise aligned t-SNE plots
figure; 
scatter3(align_vh(:,1),align_vh(:,2),align_vh(:,3),10,schaef_id,'filled'); 
colormap(myColorMap); colorbar;

figure; 
scatter3(align_nvh(:,1),align_nvh(:,2),align_nvh(:,3),10,schaef_id,'filled'); 
colormap(myColorMap); colorbar;

figure; 
scatter3(align_con(:,1),align_con(:,2),align_con(:,3),10,schaef_id,'filled'); 
colormap(myColorMap); colorbar;

% Visualise aligned t-SNE plots together
figure; 
scatter3(align_vh(:,1),align_vh(:,2),align_vh(:,3),10,'filled');
hold on
scatter3(align_nvh(:,1),align_nvh(:,2),align_nvh(:,3),10,'filled');
scatter3(align_con(:,1),align_con(:,2),align_con(:,3),10,'filled');


%% Within vs between network Euclidean distances

load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\between.mat');
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\within.mat');

% Calculate distance maps for each group average
% Control
for ii = 1:size(tsne_con,1)
    for jj = 1:size(tsne_con,1)
        eucd_con(ii,jj) = pdist2(tsne_con(ii,:),tsne_con(jj,:),'euclidean');
    end
end

% NonVH
for ii = 1:size(tsne_nvh,1)
    for jj = 1:size(tsne_nvh,1)
        eucd_nvh(ii,jj) = pdist2(tsne_nvh(ii,:),tsne_nvh(jj,:),'euclidean');
    end
end

% VH
for ii = 1:size(tsne_vh,1)
    for jj = 1:size(tsne_vh,1)
        eucd_vh(ii,jj) = pdist2(tsne_vh(ii,:),tsne_vh(jj,:),'euclidean');
    end
end

% Take upper triangles
idx = logical(triu(ones(size(between)),1));
within_upper = within(idx)';
between_upper = between(idx)';
con_upper = eucd_con(idx)';
nvh_upper = eucd_nvh(idx)';
vh_upper = eucd_vh(idx)';

% Isolate distances within a network
tsne_within(:,1) = con_upper(within_upper==1);
tsne_within(:,2) = nvh_upper(within_upper==1);
tsne_within(:,3) = vh_upper(within_upper==1);

% Isolate distances between newtorks
tsne_between(:,1) = con_upper(between_upper==1);
tsne_between(:,2) = nvh_upper(between_upper==1);
tsne_between(:,3) = vh_upper(between_upper==1);

% Plot histograms
% Find min max values
min(tsne_within,[],'all') % 0.1176
min(tsne_between,[],'all') % 0.1860
max(tsne_within,[],'all') % 17.0623
max(tsne_between,[],'all') % 18.1999

min_val = 0.1;
max_val = 19;

% Control histograms
figure; histogram(tsne_within(:,1),'Normalization','probability',...
    'BinWidth',0.5);
xlim([min_val max_val]);

figure; histogram(tsne_between(:,1),'Normalization','probability',...
    'BinWidth',0.5,'FaceColor','#D95319');
xlim([min_val max_val])
set(gca, 'YDir','reverse')

% NonVH histograms
figure; histogram(tsne_within(:,2),'Normalization','probability',...
    'BinWidth',0.5);
xlim([min_val max_val]);

figure; histogram(tsne_between(:,2),'Normalization','probability',...
    'BinWidth',0.5,'FaceColor','#D95319');
xlim([min_val max_val])
set(gca, 'YDir','reverse')

% VH histograms
figure; histogram(tsne_within(:,3),'Normalization','probability',...
    'BinWidth',0.5);
xlim([min_val max_val]);

figure; histogram(tsne_between(:,3),'Normalization','probability',...
    'BinWidth',0.5,'FaceColor','#D95319');
xlim([min_val max_val])
set(gca, 'YDir','reverse')