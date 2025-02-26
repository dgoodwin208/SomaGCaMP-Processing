rootdir = 'Processed results/';

%Mouse gcamp data
gcamp_data_files = {
    'mouse_gcamp_7136acsf_mocorr_corrhistograms.mat',...
    'mouse_gcamp_11815-2_mocorr_corrhistograms.mat',...
    'mouse_gcamp_11815-1_mocorr_corrhistograms.mat',...
    'mouse_gcamp_7136sch_mocorr_corrhistograms.mat'...
    'mouse_gcamp_4003_mocorr_corrhistograms.mat',...
    'mouse_gcamp_2565_mocorr_corrhistograms.mat',...
%     'mouse_gcamp_7584_mocorr_corrhistograms.mat',...
    };
somagcamp_data_files = {
    'mouse_somagcamp_032818_mocorr_corrhistograms.mat',...
    'mouse_somagcamp_032918-1_mocorr_corrhistograms.mat',...
    'mouse_somagcamp_032918-2_mocorr_corrhistograms.mat',...
    'mouse_somagcamp_8363_mocorr_corrhistograms.mat'};
 
 

MAX_HISTOGRAM_DISTANCE_UM = 300;
CORR_MIN = -.7; CORR_MAX = 1.;


%% Make combined somagcamp figures


distances_somagcamp = [];
corrvector_somagcamp = [];
N_somagcamp = 0;
for i = 1:length(somagcamp_data_files)
    load(fullfile(rootdir,somagcamp_data_files{i}))

    %CHOOOSE! 
    %corr_matrix_bgsub, correlations using CNMF's raw.df
    %corr_matrix_nobg, correlations using raw pixel values 
    %corr_matrix_bgsub_cainferred, correlations using CNMF's inferred.dfof
    corr_matrix = corr_matrix_nobg  ; %    corr_matrix_bgsub_cainferred
     
    N_somagcamp = N_somagcamp + size(corr_matrix,1);
    [distmatrix_masked,corrmatrix_masked] = deduplicateMatrix(distances,corr_matrix,1.);
    
    cropable_indices = distmatrix_masked>MAX_HISTOGRAM_DISTANCE_UM;
    distmatrix_masked(cropable_indices)=[];corrmatrix_masked(cropable_indices)=[];
    distances_somagcamp = [distances_somagcamp; distmatrix_masked(:)];
    corrvector_somagcamp = [corrvector_somagcamp; corrmatrix_masked(:)];
    

end


f = figure();
subplot(1,2,2);
X = [distances_somagcamp(:),corrvector_somagcamp(:)];
NPairs_somagcamp = size(X,1);
edges{2} = linspace(CORR_MIN,CORR_MAX,30);edges{1} = linspace(0,300,30);
n = hist3(X,'Edges',edges);
n1 = n';
xb = linspace(min(X(:,1)),max(X(:,1)),size(n,1));
yb = linspace(CORR_MIN,CORR_MAX,size(n,1));
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot) % heat map
colorbar;

max_somagcamp_count = max(n1(:));
% title('SomaGCaMP Correlation vs Spatial distance');
title(sprintf('SomaGCaMP \\rho vs distance (N=%i Fish, %i Neurons)',length(somagcamp_data_files),N_somagcamp));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');

%% Make combined gcamp figures
% figure;
markers = {'b+','m*','go','k.','md','kx','yd'};
distances_gcamp = [];
corrvector_gcamp = [];
N_gcamp = 0;
for i = 1:length(gcamp_data_files)
    load(fullfile(rootdir,gcamp_data_files{i}))
    %CHOOOSE! corr_matrix_bgsub, corr_matrix_nobg,corr_matrix_bgsub_cainferred
%     corr_matrix = corr_matrix_nobg; %corr_matrix_bgsub_cainferred;
    corr_matrix =  corr_matrix_nobg  ; %      corr_matrix_bgsub_cainferred
    N_gcamp = N_gcamp + size(corr_matrix,1);
    %Because the NxN matrices are symmetric, we don't want to report the
    %doubled results of neuron i vs neuron j AND neuron j vs neuron i.
    [distmatrix_masked,corrmatrix_masked] = deduplicateMatrix(distances,corr_matrix,1.);
    %Impose a maximum distance between the neurons above which we're not
    %interested in the correlation. (keeps the histogram neater)
    cropable_indices = distmatrix_masked>MAX_HISTOGRAM_DISTANCE_UM;
    distmatrix_masked(cropable_indices)=[];corrmatrix_masked(cropable_indices)=[];
    
    %Concatenate the results from multiple expeimrents
    distances_gcamp = [distances_gcamp; distmatrix_masked(:)];
    corrvector_gcamp = [corrvector_gcamp; corrmatrix_masked(:)];
    
    %This is pretty ugly and can be commented out, but this just overlays
    %the various experiments onto one plot before creating the 2D heatmap
%     plot(distmatrix_masked(:),corrmatrix_masked(:),markers{i},'MarkerSize',5);
%     ylabel('Pearson Correlation Coefficient');
%     xlabel('Spatial Distance (microns)');
% 
% 
%     title('GGaMP Correlation vs Spatial distance');
% 
%     hold on;
end
   

%To match the same number of cells for somagcamp and gcamp, we randomly
%choose a subset of the gcamp data
subset_indices = randperm(length(distances_gcamp(:)),NPairs_somagcamp);
figure(f);
subplot(1,2,1);
% X = [distances_gcamp(:),corrvector_gcamp(:)];
X = [distances_gcamp(subset_indices),corrvector_gcamp(subset_indices)];
edges{2} = linspace(CORR_MIN,CORR_MAX,30);edges{1} = linspace(0,300,30);
n = hist3(X,'Edges',edges);
n1 = n';

xb = linspace(min(X(:,1)),max(X(:,1)),size(n,1));
yb = linspace(CORR_MIN,CORR_MAX,size(n,1));
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
caxis([0 max_somagcamp_count]);
colormap(hot) % heat map
colorbar 

title(sprintf('GCaMP6f \\rho vs distance (N=%i Fish, %i Neurons)',length(gcamp_data_files),N_gcamp));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');





saveas(gcf,sprintf('combined_2DHistograms_CNMF.eps',rootdir),'epsc')

%%
[H, pValue, KSstatistic] = kstest_2s_2d([distances_somagcamp(:),corrvector_somagcamp(:)],[distances_gcamp(:),corrvector_gcamp(:)])

fprintf('SomaGCaMP: %f +- %f (SEM)\n', mean(corrvector_somagcamp), ...
   std(corrvector_somagcamp)/sqrt(length(corrvector_somagcamp)))

fprintf('GCaMP: %f +- %f (SEM)\n', mean(corrvector_gcamp), ...
   std(corrvector_gcamp)/sqrt(length(corrvector_gcamp)))


%% Draw a 



save('gcampcamp+somagcamp_vectors_cnmf.mat','distances_somagcamp',...
   'corrvector_somagcamp','distances_gcamp','corrvector_gcamp')