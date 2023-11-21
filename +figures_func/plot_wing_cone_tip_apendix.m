function plot_wing_cone_tip_apendix(hull,hull3d,fr,wingname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
frm = find(fr == hull.frames);
ind2mm = 50e-6*1000;


hull.(wingname).hull3d = double(hull3d.(wingname).hull.hull3d{frm});
hull.body.hull3d = double(hull3d.body.hull{frm});

hull.body.vectors.X = hull.body.vectors.X(frm,:);
hull.EstTip_calcSpanChord('rightwing','name2save','span','plotcone',1);
wingsname = {'rightwing','leftwing'};
otherwing = wingsname{(1 - strcmp(wingsname,wingname)*1) == 1};
wingplt = (hull.rotmat_EWtoL*(double(hull3d.(otherwing).hull.hull3d{frm})'))';
hold on;plot3(ind2mm*wingplt(:,1),ind2mm*wingplt(:,2),...
    ind2mm*wingplt(:,3),'b.','markersize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');

tip = hull.rotmat_EWtoL*hull.(wingname).coords.tip(frm,:)';
hold on;plot3(ind2mm*tip(1),ind2mm*tip(2),...
    ind2mm*tip(3),'k.','markersize',20,'MarkerFaceColor','y','MarkerEdgeColor','k');
box on;grid on;
xlabel('mm');ylabel('mm');
zlabel('mm');
axis equal
end

