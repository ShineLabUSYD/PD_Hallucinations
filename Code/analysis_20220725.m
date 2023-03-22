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

% Setting up VH IDs
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

% NonVH IDs
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

% Controls ID
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

% Run t-SNE for all subjects
for ii = 1:length(storage)
    rng(1234);
    [tsne_out,loss] = tsne(storage(ii).data,'Algorithm','barneshut','Distance',...
                'euclidean','NumDimensions',3,'NumPCAComponents',3,...
                'Perplexity',90,'LearnRate',500,...
                'Exaggeration',50);
    storage(ii).tsne = tsne_out;
    storage(ii).tsne_loss = loss;
end

% Calculate euclidean distance between ROIs
for ii = 1:length(storage)
    temp = storage(ii).tsne;
    eucd_map = zeros(400,400);
    for a = 1:length(temp)
        for b = 1:length(temp)
            eucd_map(a,b) = pdist2(temp(a,:),temp(b,:),'euclidean');
        end
    end
    storage(ii).distmap = eucd_map;
end

% Creating group matrices
% Correlation Matrices
a = 1;
b = 1;
c = 1;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,"Control")
        control_corr(:,:,a) = storage(ii).data;
        a = a + 1;
    elseif isequal(storage(ii).Index,"NonVH")
        nonvh_corr(:,:,b) = storage(ii).data;
        b = b + 1;
    elseif isequal(storage(ii).Index,"VH")
        vh_corr(:,:,c) = storage(ii).data;
        c = c + 1;
    end
end


% Euclidean Distance
a = 1;
b = 1;
c = 1;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,"Control")
        control_eucd(:,:,a) = storage(ii).distmap;
        a = a + 1;
    elseif isequal(storage(ii).Index,"NonVH")
        nonvh_eucd(:,:,b) = storage(ii).distmap;
        b = b + 1;
    elseif isequal(storage(ii).Index,"VH")
        vh_eucd(:,:,c) = storage(ii).distmap;
        c = c + 1;
    end
end

% Patient groups
patient_corr = cat(3,nonvh_corr,vh_corr);
patient_eucd = cat(3,nonvh_eucd,vh_eucd);

% Upper triangle
for ii = 1:length(storage)
    idx = logical(triu(ones(size(storage(ii).data)),1));
    v = storage(ii).data(idx)';
    storage(ii).upper_corr = v;
end

for ii = 1:length(storage)
    idx = logical(triu(ones(size(storage(ii).distmap)),1));
    v = storage(ii).distmap(idx)';
    storage(ii).upper_eucd = v;
end

% Make group upper matrices
a = 1;
b = 1;
c = 1;
for ii = 1:length(storage)
    if isequal(storage(ii).Index,"Control")
        control_corr_upper(:,a) = storage(ii).upper_corr';
        %control_eucd_upper(:,a) = storage(ii).upper_eucd';
        a = a + 1;
    elseif isequal(storage(ii).Index,"NonVH")
        nonvh_corr_upper(:,b) = storage(ii).upper_corr';
        %nonvh_eucd_upper(:,b) = storage(ii).upper_eucd';
        b = b + 1;
    elseif isequal(storage(ii).Index,"VH")
        vh_corr_upper(:,c) = storage(ii).upper_corr';
        %vh_eucd_upper(:,c) = storage(ii).upper_eucd';
        c = c + 1;
    end
end
patient_corr_upper = [nonvh_corr_upper vh_corr_upper];
patient_eucd_upper = [nonvh_eucd_upper vh_eucd_upper];


%% Formatting Matrices

% Upper Triangle groups

% Concatenated matrices
cat_control_corr = permute(control_corr,[1 3 2]);
cat_control_corr = reshape(cat_control_corr,400*19,400);
cat_patient_corr = permute(patient_corr,[1 3 2]);
cat_patient_corr = reshape(cat_patient_corr,400*77,400);

cat_control_eucd = permute(control_eucd,[1 3 2]);
cat_control_eucd = reshape(cat_control_eucd,400*19,400);
cat_patient_eucd = permute(patient_eucd,[1 3 2]);
cat_patient_eucd = reshape(cat_patient_eucd,400*77,400);

cat_nonvh_corr = permute(nonvh_corr,[1 3 2]);
cat_nonvh_corr = reshape(cat_nonvh_corr,400*46,400);
cat_vh_corr = permute(vh_corr,[1 3 2]);
cat_vh_corr = reshape(cat_vh_corr,400*31,400);

cat_nonvh_eucd = permute(nonvh_eucd,[1 3 2]);
cat_nonvh_eucd = reshape(cat_nonvh_eucd,400*46,400);
cat_vh_eucd = permute(vh_eucd,[1 3 2]);
cat_vh_eucd = reshape(cat_vh_eucd,400*31,400);

% Group average matrices
mean_control_corr = mean(control_corr,3);
mean_control_eucd = mean(control_eucd,3);
mean_patient_corr = mean(patient_corr,3);
mean_patient_eucd = mean(patient_eucd,3);
mean_diff_conpd_corr = mean_patient_corr - mean_control_corr;
mean_diff_conpd_eucd = mean_patient_eucd - mean_control_eucd;

mean_nonvh_corr = mean(nonvh_corr,3);
mean_nonvh_eucd = mean(nonvh_eucd,3);
mean_vh_corr = mean(vh_corr,3);
mean_vh_eucd = mean(vh_eucd,3);
mean_diff_nvhvh_corr = mean_vh_corr - mean_nonvh_corr;
mean_diff_nvhvh_eucd = mean_vh_eucd - mean_nonvh_eucd;

% Group average vectors
mean_control_corr2 = mean(control_corr,[2 3]);
mean_control_eucd2 = mean(control_eucd,[2 3]);
mean_patient_corr2 = mean(patient_corr,[2 3]);
mean_patient_eucd2 = mean(patient_eucd,[2 3]);

mean_nonvh_corr2 = mean(nonvh_corr,[2 3]);
mean_nonvh_eucd2 = mean(nonvh_eucd,[2 3]);
mean_vh_corr2 = mean(vh_corr,[2 3]);
mean_vh_eucd2 = mean(vh_eucd,[2 3]);

mean_diff_conpd_corr2 = mean(mean_diff_conpd_corr,2);
mean_diff_conpd_eucd2 = mean(mean_diff_conpd_eucd,2);
mean_diff_nvhvh_corr2 = mean(mean_diff_nvhvh_corr,2);
mean_diff_nvhvh_eucd2 = mean(mean_diff_nvhvh_eucd,2);


%% Permutation testing

% Control vs Patient
% Correlation Matrix
corr_perm = zeros(size(control_corr_upper,1),2);
for ii = 1:length(corr_perm)
    disp(ii);
    [sig,pval] = perm_1d_delta(control_corr_upper(ii,:)',patient_corr_upper(ii,:)',5000);
    corr_perm(ii,1) = sig;
    corr_perm(ii,2) = pval;
end

% Euclidean Matrix
eucd_perm = zeros(size(control_eucd_upper,1),2);
for ii = 1:length(eucd_perm)
    disp(ii);
    [sig,pval] = perm_1d_delta(control_eucd_upper(ii,:)',patient_eucd_upper(ii,:)',5000);
    eucd_perm(ii,1) = sig;
    eucd_perm(ii,2) = pval;
end

% NonVH vs VH
% Correlation Matrix
corr_perm_nvhvh = zeros(size(nonvh_corr_upper,1),2);
for ii = 1:length(corr_perm_nvhvh)
    disp(ii);
    [sig,pval] = perm_1d_delta(nonvh_corr_upper(ii,:)',vh_corr_upper(ii,:)',5000);
    corr_perm_nvhvh(ii,1) = sig;
    corr_perm_nvhvh(ii,2) = pval;
end

% Euclidean Matrix
eucd_perm_nvhvh = zeros(size(nonvh_eucd_upper,1),2);
for ii = 1:length(eucd_perm_nvhvh)
    disp(ii);
    [sig,pval] = perm_1d_delta(nonvh_eucd_upper(ii,:)',vh_eucd_upper(ii,:)',5000);
    eucd_perm_nvhvh(ii,1) = sig;
    eucd_perm_nvhvh(ii,2) = pval;
end


%% Reformatting into square matrix

% Load data
load('D:\PhD\PD_Hallucinations\corr_perm_conpd.mat');
load('D:\PhD\PD_Hallucinations\eucd_perm_conpd.mat');
load('D:\PhD\PD_Hallucinations\corr_perm_nvhvh.mat');
load('D:\PhD\PD_Hallucinations\eucd_perm_nvhvh.mat');


% Control vs Patient
% Correlation Matrix
temp = zeros(400,400);
roi_load = corr_perm(:,1);
idx_upper = logical(triu(ones(size(temp)),1));
temp(idx_upper) = roi_load;
corr_perm2 = temp + temp';
figure; imagesc(corr_perm2); %colormap(bluewhitered());
title("Control vs Patient Correlation")
% Corr matrix ordered by Yeo 7 Network
[~,order] = sort(schaef_id);
corr_perm2_order = corr_perm2(order,order);

for xx = 1:7
    for yy = 1:7
        random(xx,yy) = mean(mean(corr_perm2(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(random); %colormap(bluewhitered()); colorbar;
%title("Control vs Patient Overlap")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% Euclidean Matrix
temp = zeros(400,400);
roi_load = eucd_perm(:,1);
idx_upper = logical(triu(ones(size(temp)),1));
temp(idx_upper) = roi_load;
eucd_perm2 = temp + temp';
figure; imagesc(-eucd_perm2); colormap(bluewhitered());
title("Control vs Patient Euclidean Distance")
[~,order] = sort(schaef_id);
eucd_perm2_order = eucd_perm2(order,order);


% NonVH vs VH
% Correlation Matrix
temp = zeros(400,400);
roi_load = corr_perm_nvhvh(:,1);
idx_upper = logical(triu(ones(size(temp)),1));
temp(idx_upper) = roi_load;
corr_perm_nvhvh2 = temp + temp';
figure; imagesc(corr_perm_nvhvh2); colormap(bluewhitered());
title("NonVH vs VH Correlation")

% Euclidean Matrix
temp = zeros(400,400);
roi_load = eucd_perm_nvhvh(:,1);
idx_upper = logical(triu(ones(size(temp)),1));
temp(idx_upper) = roi_load;
eucd_perm_nvhvh2 = temp + temp';
figure; imagesc(eucd_perm_nvhvh2); colormap(bluewhitered());
title("NonVH vs VH Euclidean Distance")

% Overlapping differences
overlap_conpd = corr_perm2==1 & eucd_perm2==1;
figure; imagesc(overlap_conpd); %colormap(bluewhitered());
title("Control vs Patient Overlap") 

for xx = 1:7
    for yy = 1:7
        overlap_conpd_schaef(xx,yy) = mean(mean(overlap_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(overlap_conpd_schaef); colormap(bluewhitered()); colorbar;
title("Control vs Patient Overlap")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

overlap_nvhvh = corr_perm_nvhvh2==1 & eucd_perm_nvhvh2==1;
figure; imagesc(overlap_nvhvh); %colormap(bluewhitered());
title("NonVH vs VH Overlap")

% Differences
% Control vs Patient
diff_conpd = corr_perm2 - eucd_perm2;
figure; imagesc(diff_conpd); colormap(bluewhitered());
title("Control vs Patient Differences")

for xx = 1:7
    for yy = 1:7
        diff_conpd_schaef(xx,yy) = mean(mean(diff_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(diff_conpd_schaef); colormap(bluewhitered()); colorbar;
title("Control vs Patient Differences")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% NonVH vs VH
diff_nvhvh = corr_perm_nvhvh2 - eucd_perm_nvhvh2;
figure; imagesc(diff_nvhvh); colormap(bluewhitered());
title("NonVH vs VH Differences")

for xx = 1:7
    for yy = 1:7
        diff_nvhvh_schaef(xx,yy) = mean(mean(diff_nvhvh(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(diff_nvhvh_schaef); colormap(bluewhitered()); colorbar;
title("NonVH vs VH Differences")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% ConPD vs NvhVH
overlap_diff = diff_conpd + diff_nvhvh;
figure; imagesc(overlap_diff); colormap(bluewhitered());
title("ConPD vs NvhVH Common Differences")

for xx = 1:7
    for yy = 1:7
        overlap_diff_schaef(xx,yy) = mean(mean(overlap_diff(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(overlap_diff_schaef); colormap(bluewhitered()); colorbar;
title("ConPD vs NvhVH Common Differences")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

overlap_diff_corr = overlap_diff==2;
figure; imagesc(overlap_diff_corr); colormap(bluewhitered());
title("ConPD vs NvhVH Correlation")

overlap_diff_eucd = overlap_diff==-2;
figure; imagesc(overlap_diff_eucd); colormap(bluewhitered());
title("ConPD vs NvhVH Euclidean Distance")

% Hallucinations in Euclidean Distance
eucd_diff_nvhvh = diff_nvhvh==-1;
eucd_diff_conpd = diff_conpd==-1;
eucd_diff_hallucination = eucd_diff_conpd - eucd_diff_nvhvh;
figure; imagesc(eucd_diff_hallucination); colormap(bluewhitered());
title("ConPD vs NvhVH Hallucination EUCD")

for xx = 1:7
    for yy = 1:7
        hallucination_all(xx,yy) = mean(mean(eucd_diff_hallucination(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(hallucination_all); colormap(bluewhitered()); colorbar;
title("ConPD vs NvhVH Hallucination EUCD")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

eucd_diff_hallucination2 = eucd_diff_hallucination==-1;
figure; imagesc(eucd_diff_hallucination2); colormap(bluewhitered());
title("Hallucination EUCD")

for xx = 1:7
    for yy = 1:7
        hallucination_only(xx,yy) = mean(mean(eucd_diff_hallucination2(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(hallucination_only); colormap(bluewhitered()); colorbar;
title("Hallucination EUCD")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% Hallucinations in Correlations
corr_nvhvh_only = diff_nvhvh==1;
corr_conpd_only = diff_conpd==1;
corr_diff_hallucination = corr_conpd_only - corr_nvhvh_only;
corr_hallucination = corr_diff_hallucination==-1;

% Hybrid Matrices
idx = logical(triu(ones(size(diff_conpd)),1));
v = diff_conpd(idx)';
temp = zeros(400,400);
roi_load = v;
temp(idx) = roi_load;
for ii = 1:7
    for jj = 1:7
        hybrid_diffconpd(schaef_id==ii,schaef_id==jj) = diff_conpd_schaef(ii,jj);
    end
end
hybrid_diffconpd(idx) = roi_load;
figure; imagesc(hybrid_diffconpd); colormap(bluewhitered());


% Net difference across 400 ROIs
regional_diff_conpd = mean(diff_conpd,2);
surf_schaef(regional_diff_conpd,'conpd');
regional_diff_nvhvh = mean(diff_nvhvh,2);
surf_schaef(regional_diff_nvhvh,'nvhvh');

% Hallucinations
mean_corr_hallucination = mean(corr_hallucination,2);
mean_eucd_hallucination = mean(eucd_diff_hallucination2,2);


surf_schaef(mean_corr_hallucination, 'corr');
surf_schaef(mean_eucd_hallucination, 'eucd');

% t-SNE - FC diff, hallucinations
hallucination_diff = eucd_diff_hallucination2 - corr_hallucination;
mean_hallucination_diff = mean(hallucination_diff,2);

surf_schaef(mean_hallucination_diff, 'eucd-corr');


%% Analysing Differences

% Average network differences
% Control vs Patient
% Correlation
for xx = 1:7
    for yy = 1:7
        mean_conpd_corr_schaef(xx,yy) = mean(mean(mean_diff_conpd_corr(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_corr_schaef); colormap(bluewhitered()); colorbar;
title("Control - Patient CORR")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% Euclidean Distance
for xx = 1:7
    for yy = 1:7
        mean_conpd_eucd_schaef(xx,yy) = mean(mean(mean_diff_conpd_eucd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_eucd_schaef); colormap(bluewhitered()); colorbar;
title("Control - Patient EUCD")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% NonVH vs VH
% Correlation
for xx = 1:7
    for yy = 1:7
        mean_nvhvh_corr_schaef(xx,yy) = mean(mean(mean_diff_nvhvh_corr(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_nvhvh_corr_schaef); colormap(bluewhitered()); colorbar;
title("NonVH - VH CORR")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% Euclidean Distance
for xx = 1:7
    for yy = 1:7
        mean_nvhvh_eucd_schaef(xx,yy) = mean(mean(mean_diff_nvhvh_eucd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_nvhvh_eucd_schaef); colormap(bluewhitered()); colorbar;
title("NonVH - VH EUCD")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% Control vs PD sig diff
% Correlation sig diff
roi_load = mean_diff_conpd_corr(corr_perm2==1);
corr_diff_conpd = zeros(400,400);
corr_diff_conpd(corr_perm2==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_conpd_corr_schaef(xx,yy) = mean(mean(corr_diff_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_corr_schaef); colormap(bluewhitered()); colorbar;
title("Control - Patient CORR")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

roi_load = mean_diff_conpd_eucd(corr_perm2==1);
corr_diff_conpd2 = zeros(400,400);
corr_diff_conpd2(corr_perm2==1) = roi_load;

% Euclidean Distance sig diff
roi_load = mean_diff_conpd_eucd(eucd_perm2==1);
eucd_diff_conpd = zeros(400,400);
eucd_diff_conpd(eucd_perm2==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_conpd_eucd_schaef(xx,yy) = mean(mean(eucd_diff_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_eucd_schaef); colormap(bluewhitered()); colorbar;
title("Control - Patient ")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

roi_load = mean_diff_conpd_corr(eucd_perm2==1);
eucd_diff_conpd2 = zeros(400,400);
eucd_diff_conpd2(eucd_perm2==1) = roi_load;

% Overlap
roi_load = mean_diff_conpd_corr(overlap_conpd==1);
corr_overlap_conpd = zeros(400,400);
corr_overlap_conpd(overlap_conpd==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_conpd_overlapcorr_schaef(xx,yy) = mean(mean(corr_overlap_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_overlapcorr_schaef); colormap(bluewhitered()); colorbar;
title("CORR = EUCD (CORR)")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

roi_load = mean_diff_conpd_eucd(overlap_conpd==1);
eucd_overlap_conpd = zeros(400,400);
eucd_overlap_conpd(overlap_conpd==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_conpd_overlapeucd_schaef(xx,yy) = mean(mean(eucd_overlap_conpd(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_conpd_overlapeucd_schaef); colormap(bluewhitered()); colorbar;
title("CORR = EUCD (EUCD)")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

% NonVH vs VH sig diff
% Correlation sig diff
roi_load = mean_diff_nvhvh_corr(corr_perm_nvhvh2==1);
corr_diff_nvhvh = zeros(400,400);
corr_diff_nvhvh(corr_perm_nvhvh2==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_nvhvh_corr_schaef(xx,yy) = mean(mean(corr_diff_nvhvh(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_nvhvh_corr_schaef); colormap(bluewhitered()); colorbar;
title("NonVH - VH (CORR)")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

corr_diff_nvhvh2 = zeros(400,400);
corr_diff_nvhvh2(corr_perm_nvhvh2==1) = mean_diff_nvhvh_eucd(corr_perm_nvhvh2==1);

% Euclidean distance sig diff
roi_load = mean_diff_nvhvh_eucd(eucd_perm_nvhvh2==1);
eucd_diff_nvhvh = zeros(400,400);
eucd_diff_nvhvh(eucd_perm_nvhvh2==1) = roi_load;
for xx = 1:7
    for yy = 1:7
        mean_nvhvh_eucd_schaef(xx,yy) = mean(mean(eucd_diff_nvhvh(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_nvhvh_eucd_schaef); colormap(bluewhitered()); colorbar;
title("NonVH - VH (EUCD)")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

eucd_diff_nvhvh2 = zeros(400,400);
eucd_diff_nvhvh2(eucd_perm_nvhvh2==1) = mean_diff_nvhvh_corr(eucd_perm_nvhvh2==1);

% ConPD vs NvhVH
% Distance
roi_load = mean_diff_nvhvh_eucd(eucd_diff_hallucination2==1);
eucd_hallucination = zeros(400,400);
eucd_hallucination(eucd_diff_hallucination2) = roi_load;
figure; imagesc(eucd_hallucination); colormap(bluewhitered());

for xx = 1:7
    for yy = 1:7
        mean_eucd_hallucination_schaef(xx,yy) = mean(mean(eucd_hallucination(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_eucd_hallucination_schaef); colormap(bluewhitered()); colorbar;
title("EUCD Hallucination")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});

corr_hallucination2 = zeros(400,400);
corr_hallucination2(corr_hallucination==1) = mean_diff_nvhvh_corr(corr_hallucination==1);
for xx = 1:7
    for yy = 1:7
        mean_corr_hallucination_schaef(xx,yy) = mean(mean(corr_hallucination2(schaef_id==xx,schaef_id==yy)));
    end
end
figure; imagesc(mean_corr_hallucination_schaef); colormap(bluewhitered()); colorbar;
title("CORR Hallucination")
xticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});
yticklabels({'VIS','SOM','DAN','VAN','LIM','FPN','DMN'});


% Unique differences for Correlation and tSNE (conpd)
unique_diff_conpd = corr_perm2 - eucd_perm2;

unique_pca_conpd = pca(unique_diff_conpd);
surf_schaef(unique_pca_conpd(:,1),'conpd');

mean_unique_conpd = mean(unique_diff_conpd,2);
surf_schaef(mean_unique_conpd, 'mean_conpd');

% Unique differences for Correlation and tSNE(nvhvh)
unique_diff_nvhvh = corr_perm_nvhvh2 - eucd_perm_nvhvh2;

unique_pca_nvhvh = pca(unique_diff_nvhvh);
surf_schaef(unique_pca_nvhvh(:,1),'nvhvh');

mean_unique_nvhvh = mean(unique_diff_nvhvh,2);
surf_schaef(mean_unique_nvhvh, 'mean_nvhvh');


%% Analying differences between correlation and t-SNE results (Eigendecomposition)

% Finding eigenvectors of each permutation matrix
[eig_vec_conpdcorr, ~] = eig(corr_perm2);
[eig_vec_conpdeucd, ~] = eig(eucd_perm2);

[eig_vec_nvhvhcorr, eig_val_nvhvhcorr] = eig(corr_perm_nvhvh2);
[eig_vec_nvhvheucd, ~] = eig(eucd_perm_nvhvh2);

% Sort eigenvectors
[d,ind] = sort(diag(eig_val_nvhvhcorr));
eig_val_nvhvhcorr2 = eig_val_nvhvhcorr(ind,ind);
eig_vec_nvhvhcorr2 = eig_vec_nvhvhcorr(:,ind);

% Correlate the correlation and eucd vectors
[coeff_conpd,pval_conpd] = corr(eig_vec_conpdcorr(:,1), eig_vec_conpdeucd(:,1),'Type','Pearson'); 
[coeff_nvhvh,pval_nvhvh] = corr(eig_vec_nvhvhcorr(:,1), eig_vec_nvhvheucd(:,1),'Type','Pearson');

% Plotting eigenvectors
surf_schaef(eig_vec_conpdcorr(:,1),'conpdcorr');
surf_schaef(eig_vec_conpdeucd(:,1),'conpdeucd');
surf_schaef(eig_vec_nvhvhcorr(:,1),'nvhvhcorr');
surf_schaef(eig_vec_nvhvheucd(:,1),'nvhvheucd');

% PCA
[coeff,score,~] = pca(corr_perm_nvhvh2);
surf_schaef(coeff(:,1),'pca1');


%% Comparison of Euclidean Distance with Margulies Principal Gradient

load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\grad_schaef.mat');

g1 = grad_schaef(:,1);
for ii = 1:length(g1)
    for jj = 1:length(g1)
        g1_diff(ii,jj) = abs(g1(ii) -  g1(jj));
    end
end

grad_diff_sig = g1_diff(eucd_perm_nvhvh2==1);
grad_diff_nsig = g1_diff(eucd_perm_nvhvh2==0);
[grad_sig,grad_pval] = perm_1d_delta(grad_diff_sig(:),grad_diff_nsig(:),5000);
[coeff_grad,pval_grad] = corr(grad_diff_sig, grad_diff_nsig, 'Type', 'Pearson');

mean_diff_sig = mean(grad_diff_sig);
mean_diff_nsig = mean(grad_diff_nsig);
mean_diff = mean_diff_sig - mean_diff_nsig;


grad_diff_sig = g1_diff(eucd_perm_nvhvh2==1);
eucd_diff_sig = mean_diff_eucd(eucd_perm_nvhvh2==1);
[coeff,pval] = corr(grad_diff_sig, eucd_diff_sig, 'Type', 'Pearson');
figure; scatter(eucd_diff_sig,grad_diff_sig,10,'filled');


%% Comparison of Euclidean Distance with Patient Gradients

load('D:\PhD\PD_Hallucinations\vh_gradient2.mat');
load('D:\PhD\PD_Hallucinations\nonvh_gradient2.mat');

% Make gradient roi diff matrix (roi x roi x subject)
% How far apart each region is on the unimodal-heteromodal gradient
% Smaller distance means closer together. Larger distance means further
% apart

% VH
for subject = 1:size(vh_gradient2,1)
    g1 = vh_gradient2(subject,:);
    g1_matrix = zeros(400,400);
    for ii = 1:length(g1)
        for jj = 1:length(g1)
            g1_matrix(ii,jj) = abs(g1(ii) - g1(jj));
        end
    end
    vh_grad_diff(:,:,subject) = g1_matrix;
end

% NonVH
for subject = 1:size(nonvh_gradient2,1)
    g1 = nonvh_gradient2(subject,:);
    g1_matrix = zeros(400,400);
    for ii = 1:length(g1)
        for jj = 1:length(g1)
            g1_matrix(ii,jj) = abs(g1(ii) - g1(jj));
        end
    end
    nonvh_grad_diff(:,:,subject) = g1_matrix;
end

% mean_vh_gradient = mean(vh_gradient2,1)';
% mean_nvh_gradient = mean(nonvh_gradient2,1)';
% mean_diff_gradient = mean_nvh_gradient - mean_vh_gradient;
% for ii = 1:length(mean_diff_gradient)
%     for jj = 1:length(mean_diff_gradient)
%         mean_grad_diff_matrix(ii,jj) = abs(mean_diff_gradient(ii) - mean_diff_gradient(jj));
%     end
% end

% Mean gradient distance matrix
mean_vh_gradDist = mean(vh_grad_diff,3);
mean_nvh_gradDist = mean(nonvh_grad_diff,3);
mean_diff_gradDist = mean_nvh_gradDist - mean_vh_gradDist;

% Mean eucd diff matrix
mean_vh_eucd = mean(vh_eucd,3);
mean_nvh_eucd = mean(nonvh_eucd,3);
mean_diff_eucd = mean_nvh_eucd - mean_vh_eucd;

eucd_sig = mean_diff_eucd(eucd_perm_nvhvh2==1);
grad_sig = mean_diff_gradDist(eucd_perm_nvhvh2==1);
[coeff_tsne, pval_tsne] = corr(eucd_sig,grad_sig,'Type','Pearson');
figure; scatter(eucd_sig,grad_sig,10,'filled');

mean_vh_corr = mean(vh_corr,3);
mean_nvh_corr = mean(nonvh_corr,3);
mean_diff_corr = mean_nvh_corr - mean_vh_corr;
corr_sig = mean_diff_corr(eucd_perm_nvhvh2==1);
[coeff_corr,pval_corr] = corr(corr_sig,grad_sig,'Type','Pearson');
figure; scatter(corr_sig,grad_sig,10,'filled');

% Average difference between gradient scores
mean_vh_grad = mean(vh_gradient2,1);
mean_nvh_grad = mean(nonvh_gradient2,1);
mean_diff_grad = mean_nvh_grad - mean_vh_grad;
% Average #edges different between groups t-SNE
mean_eucd_perm_nvhvh = mean(eucd_perm_nvhvh2,2);
% Correlation between gradient differences and number of different edges
[coeff,pval] = corr(mean_eucd_perm_nvhvh,mean_diff_grad','Type','Pearson'); % coeff = -0.0999, pval = 0.0458
figure; scatter(mean_eucd_perm_nvhvh,mean_diff_grad',10,'filled');

% Average #edges different between group correlation
mean_corr_perm_nvhvh = mean(corr_perm_nvhvh2,2);
[coeff,pval] = corr(mean_corr_perm_nvhvh,mean_diff_grad','Type','Pearson'); % coeff = -0.3241, pval = 3.0906e-11
figure; scatter(mean_corr_perm_nvhvh,mean_diff_grad',10,'filled');


%% Eigenvalue decomposition

% Eigenvalue decomposition for each subject correlation matrix
for ii = 1:length(storage)
    [eig_vec, eig_val] = eig(storage(ii).data);
    storage(ii).eig_vec = eig_vec;
    storage(ii).eig_val = eig_val; % To get as a vector real(diag(eig_val));
end

% Collate eigenvalues and check whether they are positive/negative
for ii = 1:length(storage)
    eig_val_all(:,ii) = real(diag(storage(ii).eig_val));
end
figure; imagesc(eig_val_all); colormap(bluewhitered()); colorbar;


%% Plotting

% Control vs PD correlation
test = mean(corr_diff_conpd,2);

% NonVH vs VH correlation
test = mean(corr_diff_nvhvh,2);

% Hallucination in Correlation only
test = mean(corr_hallucination,2);

% Control vs PD euclidean
test = mean(eucd_diff_conpd,2);

% NonVH vs VH euclidean
test = mean(eucd_diff_nvhvh,2);

% Hallucination in euclidean only
test = mean(eucd_diff_hallucination2,2);
