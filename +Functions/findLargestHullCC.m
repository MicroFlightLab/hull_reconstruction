function [largestCC, CC, idx] = findLargestHullCC (hull, conn, sizevec)
% finds the largest connected component in the given 3D hull (x,y,z)
% and returns it in largestCC.
% otherCC contains all the other voxels that were not included in largestCC
% CC contains the connected components analysis returned by bwconncomp
% idx is the indices of the connected components in CC sorted by descending
% size
%
% sizevec is the size of the 3D volume. Dafault value is [512 512 512] 
% conn is the connectivity used for finding connected components in 3D


% -----------
% DEFINITIONS
% -----------
xmin = min(hull(:,1));
ymin = min(hull(:,2));
zmin = min(hull(:,3));
% 


if (~exist('conn','var'))
    conn = 6 ;
end


VOL=Functions.ptcloud2binMat(hull+[sign(xmin)*xmin+1,sign(ymin)*ymin+1,sign(zmin)*zmin+1]);% BUILD 3D VOLUME FROM HULL
CC = bwconncomp(VOL,conn);
% ----------------------------
% ANALYZE CONNECTED COMPONENTS
% ----------------------------

if (CC.NumObjects>1)
    
    SizeCC=cell2mat(cellfun(@(x) length(x),CC.PixelIdxList,'uni',0));
    CC.NumObjects=sum(SizeCC>50);
    CC.PixelIdxList=CC.PixelIdxList(SizeCC>50);
    SizeCC=SizeCC(SizeCC>50);
    
    if isempty(SizeCC)==1
     largestCC = hull ;
     idx = 1 ;
     return
    end
        
    if (CC.NumObjects>1)
    [~, idx] = sort(SizeCC,'descend') ;
    else
        idx=1;
    end
        
        
    
    [idx1, idx2, idx3] = ind2sub(size(VOL), CC.PixelIdxList{idx(1)}) ;
    largestCC = double([idx1 idx2 idx3 ]) ;
    largestCC = largestCC-double([sign(xmin)*xmin+1,sign(ymin)*ymin+1,sign(zmin)*zmin+1]) ; % xxx
 
else
    % if there is only one connected component
    largestCC = hull ;
    idx = 1 ;
end



return
