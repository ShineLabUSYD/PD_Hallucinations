%% Setting up data folder pathway and checking contents

% Specify the folder where the files lives.
myFolder = 'D:\Honours\Research\preprocessing_20211202_ICAAROMA\FunctionalConnectivity';
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

%% Calculate Mean Connectivity Matrix

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*ts_corr_schaef*.mat'); % Looking for the 400x400 matrices
theFiles = dir(filePattern); % Storing all filenames

% Loads in all the data into one variable
for k = 1:length(theFiles)
    % Set up pathway so that it is dynamic and changes with each iteration
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Loading in the subject matrix
    % Variable containing subject IDs to keep track for which matrix
    % belongs to which individual (thisArray)
    thisArray{k} = load(fullFileName); % Should have all the patient matrices
    % Variable for matrix to be used when calculating average matrix
    averageArray = load(fullFileName); 
    if k == 1
        sumArray = averageArray.ts_corr;
    else
        sumArray = sumArray + averageArray.ts_corr;
    end
end
% Calculating group-level (population) average matrix
meanArray = sumArray / length(theFiles);
% Visualise meanArray connectivity matrix
figure; imagesc(meanArray);colorbar;


%% Construct Population Average Gradient (Brainspace Toolbox Tutorial 1)

% Load Schaefer parcellation
labeling = load_parcellation('schaefer',400);
% Load brings in array as structure. Now we bring it into numeric array
labeling = labeling.schaefer_400;
% Load the conte69 hemisphere surfaces;
[surf_lh,surf_rh] = load_conte69();

% Constructing group-level average gradient
% Making GradientMap object Kernel = normalised angle, Approach = difusion mapping
avgGM = GradientMaps('kernel','na','approach','dm');
% Fitting group-level connectivity matrix to GradientMap
% Sparsity 90%, Alpha threshold = 0.5 (these values are also set as default
% parameters therefore no need to explicitly do this)
% avgGM = avgGM.fit(meanArray,'sparsity',90,'alpha',0.5);
avgGM = avgGM.fit(meanArray,'sparsity',90,'alpha',0.5);

gradients = avgGM.gradients{1};


%% Plotting 

% Plotting top 3 gradients for group-level
plot_hemispheres(avgGM.gradients{1}(:,1:3),{surf_lh,surf_rh}, ...
    'parcellation', labeling, ...
    'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});

% Checking which gradients should be kept for further analysis
scree_plot(avgGM.lambda{1}); % Higher Eigenvalues are more important, unless using laplacian eigenmapping in which case lower values are better

% Using colour to observe multidimensional interactions between the
% gradients
% Put these colours back on the cortical surface
gradient_in_euclidean(flipgradient(:,1:2),{surf_lh,surf_rh},labeling);
gradient_in_euclidean(avgGM.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling);


%% Individual Gradient Construction & Alignment for each Subject to Population Average Gradient (Brainspace Toolbox Tutorial 2)

% Extracting only the subject IDs from the full filename and making a
% vector containing all the subject IDs
for a = 1:length(theFiles)
    fullsubject = theFiles(a).name;
    subjectname = extractBefore(fullsubject,'-ts_corr');
    subjectname = strrep(subjectname,'-','_');
    subjects(a) = {subjectname};
    storage_hallucinations(a).name = subjectname;
end

% Using extracted subject IDs to label each matrix in thisArray
% thisArray will now have 60 structures with both a matrix and
% corresponding label
for index = 1:length(subjects)
    thisArray{index}.subID = subjects{index};
end

% Loop through subject ID, making an individual gradient using the
% corresponding matrix from thisArray
% Align each individual gradient to the group-level average gradient
for subject = 1:length(subjects)
    % Same parameters as when constructing group-level gradient. Also using
    % 'pa' = Procrustes Alignment to align individual gradients to the
    % group-level gradient
    indGM{subject} = GradientMaps('kernel','na','approach','dm','alignment','pa');
    % indGM{subject} = indGM{subject}.fit({meanArray,thisArray{subject}.ts_corr},'sparsity',90);
    % indGM{subject} = indGM{subject}.fit(thisArray{subject}.ts_corr,'reference',flipgradient,'sparsity',90,'alpha',0.5);
    indGM{subject} = indGM{subject}.fit(thisArray{subject}.ts_corr,'reference',avgGM.gradients{1},'sparsity',90,'alpha',0.5);
    % indGM will have 60 GradientMap objects, in the same order as the
    % subject IDs listed in subjects and thisArray
end


%% Collating all the new aligned gradients into a separate matrix

% Extracting the aligned gradient for each subject
% Get indGM{x}.aligned{2}(:,1) for each subject and place into a new matrix
% New matrix = 400x60
for alignG = 1:length(indGM)
    gradientMatrix(alignG,:) = (indGM{alignG}.aligned{1}(:,2)); 
end

% Saving the new matrix with aligned gradients as 'gradients.csv' file
writematrix(gradientMatrix,'C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\personal\principal_gradients_20220818.csv');
% Saving subject names to be copied into the principal gradient csv file
writecell(subjects,'C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\personal\subject_names.csv');


%% Setting up VH folder pathway

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


%% Making list of subject IDs with VH Status

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


%% Matching VH subject IDs and storing index 

% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsVH)
    for c = 1:length(thisArray)
        if strcmp(subjectsVH(b),thisArray{c}.subID)
            theFilesVH(b).Index = c;
            storage_hallucinations(c).Status = 'VH';
        end
    end
end


%% Extracting aligned gradients from indGM and calculating group-level average

% Add up all the gradients for VH subjects
a = 1;
for d = 1:length(theFilesVH)
    VHIndex = theFilesVH(d).Index;
    VHGradient = indGM{VHIndex}.aligned{1};
    vh_gradient2(a,:) = VHGradient(:,2);
    if d == 1
        sumVHGradient = VHGradient;
    else
        sumVHGradient = sumVHGradient + VHGradient;
    end
    a = a + 1;
end

% Calculate average gradient for VH subjects
meanVHGradient = sumVHGradient / length(theFilesVH);


%% Plotting average VH gradient onto cortical surface

VH_plot = plot_hemispheres(meanVHGradient(:,1),{surf_lh,surf_rh}, ...
    'parcellation', labeling, ...
    'labeltext',{'VH Gradient 1'}, ...
    'view', 'lmisap');
% VH_plot.colorlimits([-0.1611 0.0488])
% gradient_in_euclidean(meanVHGradient(:,1:2),{surf_lh,surf_rh},labeling);


%% Setting up Non-VH folder pathway

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


%% Making list of subject IDs with Non-VH Status

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


%% Matching Non-VH subject IDs and storing index 

% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsNVH)
    for c = 1:length(thisArray)
        if strcmp(subjectsNVH(b),thisArray{c}.subID)
            theFilesNVH(b).Index = c;
            storage_hallucinations(c).Status = 'nonVH';
        end
    end
end


%% Extracting aligned gradients from indGM and calculating group-level average

% Add up all the gradients for nonVH subjects
a = 1;
for d = 1:length(theFilesNVH)
    NVHIndex = theFilesNVH(d).Index;
    NVHGradient = indGM{NVHIndex}.aligned{1};
    nonvh_gradient2(a,:) = NVHGradient(:,2);
    if d == 1
        sumNVHGradient = NVHGradient;
    else
        sumNVHGradient = sumNVHGradient + NVHGradient;
    end
    a = a + 1;
end

% Calculate average gradient for nonVH subjects
meanNVHGradient = sumNVHGradient / length(theFilesNVH);


%% Plotting average Non-VH gradient onto cortical surface

nonVH_plot = plot_hemispheres(meanNVHGradient(:,1),{surf_lh,surf_rh}, ...
    'parcellation', labeling, ...
    'labeltext',{'Non-VH Gradient 1'}, ...
    'views','lmisap');
% nonVH_plot.colorlimits([-0.1611 0.0488])
% gradient_in_euclidean(meanNVHGradient(:,1:2),{surf_lh,surf_rh},labeling);


%% Creating patient gradient

% Add up all the gradients for Control subjects
for d = 1:length(theFilesNVH)
    NVHIndex = theFilesNVH(d).Index;
    NVHGradient = indGM{NVHIndex}.aligned{1};
    if d == 1
        sumPatientGradient = NVHGradient;
    else
        sumPatientGradient = sumPatientGradient + NVHGradient;
    end
end
for d = 1:length(theFilesVH)
    VHIndex = theFilesVH(d).Index;
    VHGradient = indGM{VHIndex}.aligned{1};
    sumPatientGradient = sumPatientGradient + VHGradient;
end

totalPatients = length(theFilesVH) + length(theFilesNVH);
meanPatientGradient = sumPatientGradient/totalPatients;

con_patient_diff = meanConGradient - meanPatientGradient;
nvh_vh_diff = meanNVHGradient - meanVHGradient;


%% Setting up Controls folder pathway

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


%% Making list of subject IDs with Control Status

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


%% Matching Control subject IDs and storing index 

% Match subject VH IDs with arrays and store the iteration number
for b = 1:length(subjectsCon)
    for c = 1:length(thisArray)
        if strcmp(subjectsCon(b),thisArray{c}.subID)
            theFilesCon(b).Index = c;
            storage_hallucinations(c).Status = 'control';
        end
    end
end


%% Extracting aligned gradients from indGM and calculating group-level average

% Add up all the gradients for Control subjects
a = 1;
for d = 1:length(theFilesCon)
    ConIndex = theFilesCon(d).Index;
    ConGradient = indGM{ConIndex}.aligned{1};
    control_gradient(a,:) = ConGradient(:,2);
    if d == 1
        sumConGradient = ConGradient;
    else
        sumConGradient = sumConGradient + ConGradient;
    end
    a = a + 1;
end

% Calculate average gradient for Control subjects
meanConGradient = sumConGradient / length(theFilesCon);


%% Plotting average Control gradient onto cortical surface

Control_plot = plot_hemispheres(meanConGradient(:,1),{surf_lh,surf_rh}, ...
    'parcellation', labeling, ...
    'labeltext',{'Control Gradient 1'}, ...
    'views','lmisap');
% Control_plot.colorlimits([-0.1611 0.0488])
% gradient_in_euclidean(meanConGradient(:,1:2),{surf_lh,surf_rh},labeling);


%% Gradient Checking

% Check correlation of 2nd gradient and network hierarchy
load('D:\PhD\PD_Hallucinations\network_hierarchy.mat');
load('C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\schaefer_parcellation\schaef_id.mat');
[~,order] = sort(schaef_id);

mean_control = mean(control_gradient,1)';
mean_control_order = mean_control(order);

figure; scatter(mean_control_order,hierarchy,10,'filled');lsline;
[coeff,pval] = corr(mean_control_order,hierarchy,'Type','Pearson');


%% Extracting gradients for each group

writematrix(meanVHGradient,'C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\personal\principal_meanVHgradient.csv');
writematrix(meanNVHGradient,'C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\personal\principal_meanNVHgradient.csv');
writematrix(meanConGradient,'C:\Users\JoshB\OneDrive\Documents\MATLAB_Analysis\MATLAB\Functions\personal\principal_meanCongradient.csv');
