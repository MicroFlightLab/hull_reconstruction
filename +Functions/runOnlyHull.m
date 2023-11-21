function [hull_recon] = runOnlyHull(im4hull,parts2run,hullRec,seedname)

hullRec.parts2run = parts2run;%{'all','all','all','all'}
hullRec.im4hull.sprs.all = im4hull;%{seg.all{1}(frm).indIm,seg.all{2}(frm).indIm,seg.all{3}(frm).indIm,seg.all{4}(frm).indIm};

hullRec.FindSeed(seedname);
hullRec.hull_params();


% build an initial volume from "all" * 4 images

createvol = 1;
[ind0_frame,framevolume] = hullRec.createVol(hullRec.voxelSize4search);

[~,ind0_frame] =  hullRec.hull_reconstruction(hullRec.parts2run,createvol,ind0_frame);
createvol = 2;
[~,ind0_frame] =  hullRec.hull_reconstruction(hullRec.parts2run,createvol,ind0_frame);
createvol = 0;
 hull_recon = hullRec.hull_reconstruction(hullRec.parts2run,createvol,ind0_frame);

end



 