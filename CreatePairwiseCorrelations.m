


%We put our files into a format
%filename root, frames to use, pix_to_um
%The filename format is
%[area/animal]_[gcamp type]_[repetition number]_["mocorr" if it's
%been motion corrected]
experiments = {
    
{'midbrain_h2bgcamp_fish81_mocorr',1:1000,  0.56}; %

};
%% Load results
experiment_number =23 ; %define in another script
experiment_info = experiments{experiment_number};
datafolder_name = 'Raw Data/';
processedfolder_name = 'Processed results/';
filename_root = experiment_info{1};
filename_parts = split(filename_root,'_');
molecule_type = filename_parts{2};
experiment_number = filename_parts{3};
part_to_use = experiment_info{2};
pix_to_um = experiment_info{3};

parts_to_title = sprintf('%s (%s)',experiment_number,molecule_type);

%Load the raw data, assumed to be X*Y*T grayscale tif file
%Also read_file() is from the CNMF library
filename = fullfile(datafolder_name,sprintf('%s.tif',filename_root));
M1 = double(read_file(filename));
M1 = M1(:,:,part_to_use);

%Has CNMF been ran for this file? If not, re-run it now
filename_results = fullfile(processedfolder_name,sprintf('%s_CNMF_results-p1.mat',filename_root));
if exist(filename_results,'file')
    load(filename_results)
else
    fprintf('Didnt find processed results file %s, processing now...\n',filename_results);
    processCNMF;
end

%% Using the CNMF built-in tools to extract the signal

d1 = size(M1,1);
d2 = size(M1,2);
T = size(M1,3);
%Some hardcoded parameters we got from Eftychios and Andreas
extractControl.baselineRatio = 0.25;
extractControl.thr = 0.95;
Yr = reshape(M1,d1*d2,T);

[ inferred, filtered, raw ] = signalExtraction(Yr,A,C,b,f,d1,d2,extractControl);
centroids=zeros(size(raw.df,1),2);
for raw_idx = 1:size(raw.df,1)
    pts = raw.CR{raw_idx};
    if isempty(pts)
        continue;
    end
    pts = pts';
    mean_pts = mean(pts,1);
    centroids(raw_idx,:) = [mean_pts(2),mean_pts(1)];
end

fprintf('signal extraction complete\n')
% As it takes a bit of time to run signalExtraction(), sometimes it was
% worth saving the intermedia calculatino
% save(sprintf('%s/%s_extractedComponents.mat',folder_name,filename_root),'inferred','raw','filtered')

%% Manually curate and save the ROIs
%This part requires a bit of explanation. Overall, with the three different
%variants of gcamp we've used GCaMP6f, h2bGCaMP and somaGCaMP,I've noticed
%that I get a lot of false positives. In many cases, the false positive
%have a similar temporal history (ie, they are background), so the
%downstream pairwise correlation suffers. To solve that, this chunk of code
%is how I manually curated the CNMF results: I would use their built-in
%classification method, classify_components, then refine those using the
%plot_components_gui.
%The workflow that I found the fastest was to copy the ones vector,
%indices_manually_curated_prefiltereed, into excel. Then, as I go through
%the plot_components_gui (using the 'keep' vector), I would turn any 1 into
%a zero if the CNMF result looked like a false positive to me. After runnning through
%the neurons, I would then copy the modified indices_manually_curated_prefiltereed
%vector out of Excel and into matlab using the variable explorer.
%Then the end of the if statement below has a few lines to then process and
%save the manually curated results
filename_manualcuration = fullfile(processedfolder_name,sprintf('%s_CNMF_results-manually_curated.mat',filename_root));
if exist(filename_manualcuration,'file')
    fprintf('Found manual curation file %s, processing now \n',filename_manualcuration);
    load(filename_manualcuration)
else
    fprintf('Didnt find manual curation file %s, process now!\n',filename_manualcuration);
    % If you want to see all the components
    % plot_components_GUI(M1,A,C,b,f,Cn,options);
    indices_manually_curated = zeros(size(raw.df,1),1);
    % For a quick calculation: ignore the current barf();
    % indices_manually_curated = ones(size(raw.df,1),1);
    
    %Use the built-in CNMF for a quick filtering of mouse data (fish was
    %manually curated)
    [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(M1,A,C,b,f,YrA,options);
    
    
    %Use the GUI tool to load only the computationally filtered ones
    %WHEN WORKING MANUALLY:
    plot_components_GUI(M1,A(:,keep),C(keep,:),b,f,Cn,options);
    
    %Then we create a vector to mask the keep, and manually prune any
    %incorrect ROIs
    indices_manually_curated_prefiltereed = ones(sum(keep),1);
    %^Load this vector, copy it into Excel, then step cell by cell to remove
    %any incorrect ROIs
    
    indices_keep = find(keep);
    
    %To manually filter, enter here
    % WHEN WORKING MANUALLY, I put this line in to crash the code so you
    % can stop to do the manual curation using the plot_components_gui
    % which has already been opened
    barf();
    
    %After finishing the manual review, run these lines:
    indices_keep_manually_filtered = indices_keep(logical(indices_manually_curated_prefiltereed));
    indices_manually_curated(indices_keep_manually_filtered)=1;
    save(filename_manualcuration,'indices_manually_curated');
end



indices_legit = find(indices_manually_curated>0);

indices_manually_curated = zeros(length(indices_manually_curated),1);
indices_manually_curated(indices_legit)=1;

indices_filtered = 1:size(raw.df,1);

indices_filtered = indices_filtered (logical(indices_manually_curated));
%The number of ROIs is the number that passed the filtering
nr = length(indices_filtered);

%Apply the manual filters to the extracted components
raw.df = raw.df(indices_filtered,:);
raw.dfof = raw.dfof(indices_filtered,:);
raw.CR = raw.CR(indices_filtered);
raw.bruteextraction = zeros(size(raw.dfof));

inferred.df = inferred.df(indices_filtered,:);
inferred.dfof = inferred.dfof(indices_filtered,:);
filtered.df = filtered.df(indices_filtered,:);
filtered.dfof = filtered.dfof(indices_filtered,:);
centroids = centroids(indices_filtered,:);

%% Using the mask to create a totally raw segmentation
% This is used to get the raw pixel values out of the tif file. We do this
% because I think even the raw.df variables in CNMF have been background
% subtracted, which is not what we want for the most naive of comparisons
% between molecule types

for r_idx = 1:length(raw.CR)
    xy = raw.CR{r_idx};
    total_signal = zeros(1,size(raw.bruteextraction,2));
    
    for p = 1:size(xy,2)
        total_signal = total_signal + squeeze(M1(xy(1,p),xy(2,p),:))';
    end
    raw.bruteextraction(r_idx,:) = total_signal;
    %     mask_img = zeros(size(M1,1),size(M1,2));
    %     indices1D = sub2ind(size(test_img),xy(1,:),xy(2,:));
    %     mask_img(indices1D) = 1;
    
end
figure;
subplot(1,2,1);
plot(raw.bruteextraction');
title('Traces of raw pixel data taken from neuron ROIs');
ylabel('Pixel intensity');
subplot(1,2,2);
plot(raw.df');
title('Traces of neuron ROIS from CNMF background subtraction');
ylabel('Pixel intensity');

%% Make a correlation matrix and a distance matrix

distances = zeros(nr,nr);
for y = 1:nr
    for x = 1:nr
        distances(y,x) = sqrt(sum((centroids(y,:)-centroids(x,:)).^2))*pix_to_um;
    end
end


corr_matrix_nobg = zeros(nr,nr);
corr_matrix_bgsub = zeros(nr,nr);
corr_matrix_bgsub_cainferred = zeros(nr,nr);
for y = 1:nr
    for x = 1:nr
        corr_matrix_bgsub_cainferred(y,x) = corr2(inferred.dfof(y,:),inferred.dfof(x,:));
        corr_matrix_bgsub(y,x) = corr2(raw.df(y,:),raw.df(x,:));
        corr_matrix_nobg(y,x) = corr2(raw.bruteextraction(y,:),raw.bruteextraction(x,:));
    end
end


figure;
subplot(1,3,1)
make2DHistogram(distances(:),corr_matrix_nobg(:))
title(sprintf('%s Corr Dist raw pixels of neurons',parts_to_title));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');

subplot(1,3,2)
make2DHistogram(distances(:),corr_matrix_bgsub(:))
title(sprintf('%s Corr Dist with background subtraction',parts_to_title));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');

subplot(1,3,3)
make2DHistogram(distances(:),corr_matrix_bgsub_cainferred(:))
title(sprintf('%s Corr Dist of inferred calcium conc. per neuron',parts_to_title));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');
%

outputfile = fullfile(processedfolder_name,sprintf('%s_corrhistograms.mat',filename_root));
save(outputfile,'distances','corr_matrix_nobg','corr_matrix_bgsub','corr_matrix_bgsub_cainferred','-v7.3');



