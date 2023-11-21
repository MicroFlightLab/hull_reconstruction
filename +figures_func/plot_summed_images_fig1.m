function plot_summed_images_fig1(seg,sp,flname,kcam,figurepath,sparsefile,fr,ax)

seg.BodyCM_estimate(sp,kcam);
repPix = seg.sumimages(fr,sp,kcam,'plotSummedIm',1);

% hold on;plot(seg.bodCM.xy(:,2),seg.bodCM.xy(:,1),'ow','markersize',10);
hold on;plot(seg.bodCM.xy(:,2),seg.bodCM.xy(:,1),'+w','markersize',5,'linewidth',0.5);

colormap('jet'); hc =colorbar('fontsize',8)
hc.Ticks = [0:20:60,72];

set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 10 5 4]);
axis(ax);
print(gcf,[figurepath,sparsefile,'\svg\',flname],'-dsvg');
saveas(gcf,[figurepath,sparsefile,'\fig\',flname,'.fig'])
exportgraphics(gcf,[figurepath,sparsefile,'\pdf\',flname,'.pdf']);

end

