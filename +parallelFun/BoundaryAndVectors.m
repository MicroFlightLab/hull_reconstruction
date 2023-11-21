function [outputCell] = BoundaryAndVectors(hull,realC,winghull)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

try
hull.Xaxis(realC); % calculate X body
hull.splitwings(winghull); % define wings

wingnameCell = {'rightwing','leftwing'};

for k = 1:1:length(wingnameCell)
% calculate tip and span---------------------
[cm_coordR,wingCMR] = hull.EstTip_calcSpanChord(wingnameCell{k},'name2save','span');
hull.calculateTip(wingnameCell{k},hull.(wingnameCell{k}).hull3d,wingCMR);
[LE,TE ]= parallelFun.par_boundary(hull,realC,wingnameCell{k}); % LE TE
bound.(wingnameCell{k}) = {LE,TE};
end


outputCell = {uint16(hull.rightwing.hull3d),uint16(hull.leftwing.hull3d),(hull.rotmat_EWtoL*(hull.rightwing.vectors.span'))'...
    ,hull.rightwing.coords.tip,(hull.rotmat_EWtoL*(hull.rightwing.vectors.chord'))',(hull.rotmat_EWtoL*(hull.leftwing.vectors.span'))'...
    ,hull.leftwing.coords.tip,(hull.rotmat_EWtoL*(hull.leftwing.vectors.chord'))',hull.body.vectors.X,(bound.rightwing{1}),(bound.rightwing{2}),...
    (bound.leftwing{1}),(bound.leftwing{2}),double(hull.body.coords.CM_real),uint16(hull.rightwing.vectors.wingCM),uint16(hull.leftwing.vectors.wingCM)...
    ,uint16(hull.rightwing.vectors.alltip),uint16(hull.leftwing.vectors.alltip)};
catch
    
    outputCell = {{nan nan nan},{nan nan nan}...
        [nan nan nan],[nan nan nan],[nan nan nan],[nan nan nan],[nan nan nan]...
        ,[nan nan nan],[nan nan nan],{nan nan nan},{nan nan nan},{nan nan nan},{nan nan nan},...
        [nan nan nan],{nan nan nan},{nan nan nan},{nan nan nan},{nan nan nan}};

end
end

