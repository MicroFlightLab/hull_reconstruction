function [LE,TE] = par_boundary(hull,realC,wingname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
try
    chord_dir = (hull.(wingname).vectors.chord')'; % rotate chord to EW axis
    Tip = hull.(wingname).coords.tip;
    wing = ((double(hull.(wingname).hull3d')))';
    [LE ,TE] = hull.SplitWin_LETE(chord_dir,Tip,wing);
    [LE_srt,TE_srt] = hull.CreateBoundary(realC,LE,TE,wingname);
    hull.boundOptions(LE_srt,TE_srt,wingname);
    LE = uint16(hull.(wingname).hull.LE);
    TE = uint16(hull.(wingname).hull.TE);
    
catch
    LE = {[nan nan nan]};
    TE = {[nan nan nan]};
end
end

  