M1 = load3DTif_uint16(''); %Load the motion corrected file here

sizY = size(M1);
patch_size = [64,64];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [16,16];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 10;                                            % number of components to be found
% tau = 10; worked decently for zoomed in         % std of gaussian kernel (half size of neuron) 
tau = 5; %worked decently for whole fish data % std of gaussian kernel (half size of neuron) 
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',2,...                                % downsample in space
    'tsub',1,...                                % downsample in time
    'merge_thr',0.95,...                         % merging threshold
    'gSig',tau,... 
    'gnb',2,...                                 % number of background components
    'spatial_method','regularized'...
    );

%% run CNMF algorithm on patches and combine
Cn = correlation_image(M1);
tic;
[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(M1,K,patches,tau,p,options);
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(M1,A,C,b,f,YrA,options);
Coor = plot_contours(A,Cn,options,1); close;
traces = C + YrA;
ROIvars.fitness = compute_event_exceptionality(traces,0);
ROIvars.fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);
toc
%% set our own thresholds

keep = (ROIvars.rval_space > 0.5) | (ROIvars.fitness < - 40) | (ROIvars.fitness_delta < - 15) ;
% keep_std = traces_std>.005;
% keep = keep & keep_std;
%% show the selected components

throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')

%%
A_keep = A(:,keep);
C_keep = C(keep,:);

% plot_components_GUI(M1,A_keep,C_keep,b,f,Cn,options);

%%
% results_file = strrep(filename,'.tif',sprintf('_CNMF_results-p%i.mat',p));

save(filename_results,'A_keep','A','C','p','ROIvars','YrA','traces','C_keep','b','f','Cn','options','-v7.3');

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

% inferred.dfof = inferred.dfof(:,part_to_use);
% raw.df = raw.df(:,part_to_use);
fprintf('signal extraction complete\n')
save(sprintf('%s/%s_extractedComponents.mat',folder_name,filename_root),'inferred','raw','filtered')

