clear all;
file_folder = '../population_analysis';
filename = 'midbrain_gcamp_FISH41.tif';
Y = read_file(fullfile(file_folder,filename));

%% motion correction
Y = double(Y);
[d1,d2,T] = size(Y);
template_in = median(Y,3);
d3=1;
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',50,'max_shift',40,'us_fac',50,'iter',3);
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'d3',d3,...
            'bin_width',50,'max_shift',10,'us_fac',20,'iter',3, ...
            'grid_size',[48,48,1],'overlap_pre',[16,16,1],'mot_uf',4 ...
                );
%Running the piecewise rigid alignment
[M1,shifts1,template1] = normcorre_batch(Y,options_nonrigid,template_in);
%Alternatively, one could run a faster rigid alignment 
% [M1,shifts1,template1] = normcorre_batch(Y,options_rigid,template_in);

mocorr_file = strrep(filename,'.tif','_mocorr.tif');

save3DTif_uint16(M1,fullfile(file_folder,mocorr_file));
