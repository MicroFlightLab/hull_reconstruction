function [LE ,TE] = wing2LETE(hull3d,hull,wingname,fr)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
frm = find(fr == hull.frames);

    wing = hull3d.(wingname).hull.hull3d{frm};
    chrd = hull.(wingname).vectors.chord(frm,:);
    chord_dir = (hull.rotmat_EWtoL'* (chrd'))'; % rotate chord to EW axis
    Tip = hull.(wingname).coords.tip(frm,:);
    [LE ,TE] = hull.SplitWin_LETE(chord_dir,Tip,double(wing));
end

