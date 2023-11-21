function [body,Twowings,realC,body4plot] = par_hullRec(hullRec,parts2run,ind0_frame,varargin)
parser = inputParser;
addParameter(parser,'onlyhull',0); %  defines the size of the sphere around the CM - leaves only the heat and tail to choose refine the outcome of PCA
parse(parser, varargin{:});
body4plot = [];

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

try
   createvol = 1;
%     % build an initial volume from "all" * 4 images

   [~,ind0_frame] =  hullRec.hull_reconstruction(parts2run(1,1:4),createvol,ind0_frame);
   createvol = 2;
   [~,ind0_frame] =  hullRec.hull_reconstruction(parts2run(1,1:4),createvol,ind0_frame);

%     
    
    
    % use the already built volume to reconstruct the rest of parts2run.
    % (body * 4, wing1 all*3....)
    indname = 1; createvol =0;
    for k = 2:1:size(parts2run,1)
        if strcmp(parts2run{k,end},parts2run{k - 1,end}) == 0
            indname = 1;
        end
        
        hull_recon = hullRec.hull_reconstruction(parts2run(k,1:size(parts2run,2) - 1),createvol,ind0_frame);
        hullRec.framehull.(parts2run{k,end}){indname} = uint16(hull_recon);
        indname = indname + 1;
        
    end
    if parser.Results.onlyhull == 0
        hullRec.substrBody('uglywing',1);  
        body4plot = hullRec.framehull.body4plot;   
        Twowings = hullRec.framehull.Twowings; 
    else
        Twowings = {unique([cell2mat(hullRec.framehull.wing1');cell2mat(hullRec.framehull.wing2')],'rows')}; 
    end  
    body = hullRec.framehull.body3;
    realC = {hullRec.real_coord};
    hullRec.framehull=[];
    hullRec.im4hull.sprs=[];
catch
    body = {[nan nan nan]};
    Twowings = {[nan nan nan]};
    realC = {[nan nan nan]};
    body4plot = {[nan nan nan]};
    
end
end

