function [body,wing1,wing2,all,frames] = par_seg(sp,kcam,segcls,varargin)
% run parallel code for segmenting the wings and body. 
% input: sp - cell array of the sparse movies
% kcam  - camera number
% segcls - initilized segmentation class 

segcls.flyallImage(sp,kcam); % create a temporary field of the sparse movie in full binary image configurations [800 1280]
segcls.BodyCM_estimate(sp,kcam); % calculate the location of body CM and perform a polynomic estimation 
segcls.transIm_SegWings(sp,kcam); % loop on frames to run translte image and segment wing.

% save only the relevant output
body = segcls.body(kcam);
wing1 = segcls.wing1(kcam);
wing2 = segcls.wing2(kcam);
all = segcls.all(kcam);
frames = segcls.st_enfr;

end