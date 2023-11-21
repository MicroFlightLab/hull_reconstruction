function plot_fly_LETE_strk_fig2(hull,hull3d,fr,figsdir,figurepath,view_v,addSP)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
ind2mm = 50e-6*1000;
frm = find(fr == hull.frames);
wingcell = {'rightwing','leftwing'};
figure;

colormat = {'r','m';'b','c'};
hull.body.hull3d = hull3d.body.hull{frm};
hull.body.vectors.Xtmp =hull.body.vectors.X;
hull.body.vectors.X =hull.body.vectors.X(frm,:);
bod = ind2mm*(hull.rotmat_EWtoL *( double([hull3d.body.hull{frm}])'))';
CMbody = mean(bod);
realC = hull.real_coord{frm};
[head_tailCM,headVoxels,tailVoxels] =hull.Xaxis(realC,'plot',0,'save2hull',0); % calculate X body
tailVoxels = ind2mm*(hull.rotmat_EWtoL * [double(tailVoxels)]')';
headVoxels = ind2mm*(hull.rotmat_EWtoL * [double(headVoxels)]')';


a = ismember(bod,[tailVoxels;headVoxels],'rows');
middle = bod(a == 0,:) - CMbody;

for kwing = 1:1:2
    wingname = wingcell{kwing};
    wing = double(hull3d.(wingname).hull.hull3d{frm});
    chrd = hull.(wingname).vectors.chord(frm,:);
    chord_dir = (hull.rotmat_EWtoL'* (chrd'))'; % rotate chord to EW axis
    Tip = hull.(wingname).coords.tip(frm,:);
    CMwing = mean(ind2mm*(hull.rotmat_EWtoL*wing')') - CMbody
    [LE ,TE] = hull.SplitWin_LETE(chord_dir,Tip,wing);
    
    hullRealLE = ind2mm*(hull.rotmat_EWtoL *LE')'-CMbody;
    hullRealTE = ind2mm*(hull.rotmat_EWtoL *TE')'-CMbody;
    hull.(wingname).hull3d = hull3d.(wingname).hull.hull3d{frm};
    [cm_coord,wingCM,TE_chrd,LE_chrd]  = hull.EstTip_calcSpanChord(wingname);
    RealLE = ind2mm*(hull.rotmat_EWtoL*LE_chrd')'-CMbody;
    RealTE = ind2mm*(hull.rotmat_EWtoL*TE_chrd')'-CMbody;
    
    hold on;plot3(hullRealLE(:,1),hullRealLE(:,2),hullRealLE(:,3),'.','color',colormat{kwing,1});
    hold on;plot3(hullRealTE(:,1),hullRealTE(:,2),hullRealTE(:,3),'.','color',colormat{kwing,2});
    hold on;plot3(RealLE(:,1),RealLE(:,2),RealLE(:,3),'.k');
    hold on;plot3(RealTE(:,1),RealTE(:,2),RealTE(:,3),'.','color',[0.5 0.5 0.5]);  
%     hold on;quiver3(CMwing(1),CMwing(2),CMwing(3),chrd(1),chrd(2),chrd(3),1,'color',colormat{kwing,1},'linewidth',2);
end



tailmm = tailVoxels-CMbody;
headmm = headVoxels-CMbody;
cmtail = mean(tailmm);
hold on;plot3(tailmm(:,1),tailmm(:,2),tailmm(:,3),'.','color',[0,0.9,0],'markersize',5);hold on;
plot3(headmm(:,1),headmm(:,2),headmm(:,3),'.','color',[0,0.9,0],'markersize',5);
plot3(middle(:,1),middle(:,2),middle(:,3),'.','color',[0.3 0.67 0.04],'markersize',5);

% hold on;a0 = quiver3(cmtail(1),cmtail(2),cmtail(3),hull.body.vectors.X(1),hull.body.vectors.X(2),hull.body.vectors.X(3),3,'color','k','linewidth',3);
xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]');
axis equal 
grid on
% set(gcf,'renderer','painters');

if addSP == 1
% add wing stroke, get X Y Z wing stroke and plot their plane
pln =hull.body.vectors.strkPlan(frm,:);
ybody = hull.body.vectors.Y(frm,:);
% hold on;a1 = quiver3(0,0,0,ybody(1),ybody(2),ybody(3),1,'linewidth',3,'color','k','MaxHeadSize',0.8);
Zbody = hull.body.vectors.Z(frm,:);
% hold on;a2 = quiver3(0,0,0,Zbody(1),Zbody(2),Zbody(3),1,'linewidth',3,'color','k','MaxHeadSize',0.8);

Ystrk = ybody;
Xstrk = cross(Ystrk,pln);
sqr = [Xstrk-Ystrk;Xstrk+Ystrk;-Xstrk+Ystrk;-Xstrk-Ystrk];
X = sqr(:,1);
Y = sqr(:,2);
Z = -1/pln(3)*(pln(1)*X + pln(2)*Y);


ptch = patch('XData',1.5*X,'YData',1.5*Y,'ZData',1.5*Z) ;
ptch.FaceColor = [0.5 0 1];
ptch.FaceAlpha = 0.1;
hold on;quiver3(0,0,0,pln(1),pln(2),pln(3),2,'linewidth',2,'color',[0.8 0 1]);
end
view(view_v)
box on
% axis([-2.9447    2.1175   -2.1201    2.9134   -1.2303    1.9975])
set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 10 10 8]);
legend off

flname = ['\bodyAx and sp'];
saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig']);
exportgraphics(gcf,[figurepath,figsdir,'\pdf\',flname,'.pdf']);
exportgraphics(gcf,[figurepath,figsdir,'\eps\',flname,'.eps']);

plot2svg([figurepath,figsdir,'\svg\',flname,'.svg']);


end

