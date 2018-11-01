%% Load mocorr file and CNMF resultsresults

folder_name = '../population_analysis';
filename_root = 'midbrain_somagcamp_FISH38_mocorr';
filename = sprintf('%s/%s.tif',folder_name,filename_root);
M1 = double(read_file(filename));

% Get the resolution scale for converting pixel distances to microns
pix_to_um = 0.28; %For Fish 38 only: the 0.35x lightsheet zoom, perhaps a digital zoom?

parts_to_title = 'Fish38 (SomaGCaMP)';

filename_results = sprintf('%s/%s_CNMF_results-p1.mat',folder_name,filename_root);
if exist(filename_results,'file')
    load(filename_results)
else
   fprintf('Didnt find processed results file %s, processing now...\n',filename_results);
   processCNMF; 
end

%% Trying their official tools

d1 = size(M1,1);
d2 = size(M1,2);
T = size(M1,3);
extractControl.baselineRatio = 0.25;
extractControl.thr = 0.95;
Yr = reshape(M1,d1*d2,T);

[ inferred, filtered, raw ] = signalExtraction(Yr,A,C,b,f,d1,d2,extractControl);
centroids=zeros(size(raw.df,1),2);

Cn = correlation_image(M1);
for raw_idx = 1:size(raw.df,1)
    pts = raw.CR{raw_idx,1};
    if isempty(pts)
        continue;
    end
    pts = pts';
    mean_pts = mean(pts,1);
    centroids(raw_idx,:) = [mean_pts(2),mean_pts(1)];
end


fprintf('signal extraction complete\n')
save(sprintf('%s/%s_extractedComponents.mat',folder_name,filename_root),'inferred','raw','filtered')
%% Make the projection images 

figure; 
imagesc(mean(M1,3));
axis off; colormap gray
saveas(gcf,sprintf('%s/%s_meanprojection.tif',folder_name,filename_root))
close;
figure; 
imagesc(Cn,[0 1]);
axis off; colorbar;
saveas(gcf,sprintf('%s/%s_correlationImage.eps',folder_name,filename_root),'epsc')
close;

%% Manually curate and save the ROIs

% For the zebrafish figures, we carefully went through each ROI by human
% eye to confirm it was a legitimate neuron (judged by spatial footprint
% and calcium activity), and also to remove putative neurons that weren't in the tectum 
% filename_manualcuration = sprintf('%s/%s_CNMF_results-manually_curated.mat',folder_name,filename_root);
% if exist(filename_manualcuration,'file')
%     load(filename_manualcuration)
% else
%    fprintf('Didnt find manual curation file %s, process now!\n',filename_manualcuration);
%    plot_components_GUI(M1,A,C,b,f,Cn,options);
%    indices_manually_curated = zeros(size(raw.df,1),1);
%    barf();
%    %When ready:
% %    save(filename_manualcuration,'indices_manually_curated');
% end

%For all other cases, you can just accept all pututative neurons
indices_manually_curated = ones(size(raw.df,1),1);

indices_legit = find(indices_manually_curated>0);
num_to_use = length(indices_legit);

indices_of_the_legits_to_keep = randperm(length(indices_legit),num_to_use);
indices_legit = indices_legit(indices_of_the_legits_to_keep);
indices_manually_curated = zeros(length(indices_manually_curated),1);
indices_manually_curated(indices_legit)=1;


%% Plot the the centroids of the ROIs


figure; 
roikeep = ones(size(centroids,1),1);

indices_filtered = 1:size(raw.df,1);
% indices_filtered = indices_filtered(logical(roikeep));
indices_filtered = indices_filtered (logical(indices_manually_curated));
%The number of ROIs is the number that passed the filtering
nr = length(indices_filtered);
%Plot the correlation image
imagesc(Cn,[0 1]); colormap gray;

hold on;
for c_idx = 1:nr
    plot(centroids(indices_filtered(c_idx),1),centroids(indices_filtered(c_idx),2),'ro')
end
hold off;
title(sprintf('%s centroids of activity rois,N=%i',filename,sum(indices_manually_curated)));
saveas(gcf,sprintf('%s/%s_ROIfig_manuallycurated.eps',folder_name,filename_root),'epsc')
%% Make a correlation matrix and a distance matrix
corr_matrix = zeros(nr,nr);
for y = 1:nr
    for x = 1:nr
%         corr_matrix(y,x) = corr2(inferred.dfof(indices_filtered(y),:),inferred.dfof(indices_filtered(x),:));
        corr_matrix(y,x) = corr2(raw.df(indices_filtered(y),:),raw.df(indices_filtered(x),:));
    end
end

distances = zeros(nr,nr);
for y = 1:nr
    for x = 1:nr
        distances(y,x) = sqrt(sum((centroids(indices_filtered(y),:)-centroids(indices_filtered(x),:)).^2))*pix_to_um;
    end
end

figure;
subplot(1,2,1)
plot(distances(:),corr_matrix(:),'.','MarkerSize',5);
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');


title(sprintf('%s Correlation vs Spatial distance',parts_to_title));
subplot(1,2,2);
X = [distances(:),corr_matrix(:)];
n = hist3(X,[30,30]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(X(:,1)),max(X(:,1)),size(n,1)+1);
yb = linspace(min(X(:,2)),max(X(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot) % heat map
title(sprintf('%s Correlation vs Spatial distance',parts_to_title));
ylabel('Pearson Correlation Coefficient');
xlabel('Spatial Distance (microns)');
% grid on

filename_corrVectors = sprintf('%s/%s_correlationScatterPlot_N=%i.mat',folder_name,filename_root,num_to_use);
save(filename_corrVectors,'distances','corr_matrix');    



%% Show all ROIs
stdProjc = std(M1,0,3);

figure;
plot_contours(A(:,logical(indices_manually_curated)),stdProjc,options,1);
axis off; colorbar;
saveas(gcf,sprintf('%s/%s_correlationImageWROIS.eps',folder_name,filename_root),'epsc')
close;


