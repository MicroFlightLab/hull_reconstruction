function [ind0_frame,realcoord] = par_inivol(hullRec,parts2run,ind0_frame,varargin)

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

try
   createvol = 1;
%     % build an initial volume from "all" * 4 images

   [~,ind0_frame] =  hullRec.hull_reconstruction(parts2run(1,1:4),createvol,ind0_frame);
   realcoord = hullRec.real_coord;
catch
    ind0_frame = {[nan nan nan]};
    

    
end
end