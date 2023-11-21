
function [hullReal,rtmat] = hullRec2lab(inds,Rotation_Matrix,RotMat_vol,realC)
% get hull in lab axis
rtmat =  Rotation_Matrix*RotMat_vol';
hullReal = (rtmat * ([realC{1}(inds(:,1))',realC{2}(inds(:,2))',realC{3}(inds(:,3))'])')';
end


