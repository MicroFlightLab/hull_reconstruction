function plot_boundary_strkpln_fig2(ploters,hull,hull3d,fr,figurepath,figsdir,flname,view_v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ploters.labAx_ploter(hull,fr,hull3d,'bound',1,'LETE',1,'normal',0,'cam',0,'chord',0,'strk',0,'move2CM',1);
frm = find(fr == hull.frames);


% add wing stroke, get X Y Z wing stroke and plot their plane
pln =hull.body.vectors.strkPlan(frm,:);
ybody = hull.body.vectors.Y(frm,:);
Ystrk = ybody;
Xstrk = cross(Ystrk,pln);
sqr = [Xstrk-Ystrk;Xstrk+Ystrk;-Xstrk+Ystrk;-Xstrk-Ystrk];
X = sqr(:,1);
Y = sqr(:,2);
Z = -1/pln(3)*(pln(1)*X + pln(2)*Y);
ptch = patch('XData',1.5*X,'YData',1.5*Y,'ZData',1.5*Z) ;
ptch.FaceColor = [0.5 0 1];
ptch.FaceAlpha = 0.1;
% hold on;quiver3(0,0,0,pln(1),pln(2),pln(3),2,'linewidth',2,'color',[0.8 0 1]);

strkvec = [0,0,0;pln*2];
plot3(strkvec(:,1),strkvec(:,2),strkvec(:,3),'linewidth',2,'color',[0.8 0 1]);

box on
legend off
view(view_v)
title('')
legend off
axis([-2.5    1   -1.2    1.5   -1.2303    1.5])

set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 10 10 8]);
set(gcf,'renderer','painters');
exportgraphics(gcf,[figurepath,figsdir,'\eps\',flname,'.eps']);

% print(gcf,[figurepath,figsdir,'\svg\',flname],'-dsvg')
% saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig'])
% exportgraphics(gcf,[figurepath,figsdir,'\pdf\',flname,'.pdf']);

end

