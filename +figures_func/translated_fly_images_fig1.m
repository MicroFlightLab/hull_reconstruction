function ax = translated_fly_images_fig1(kcam,flname,seg,sp,fr,figurepath,sparsefile)
% plot the translated fly images.

seg.flyallImage(sp,kcam); % create a temporary field of the sparse movie in full binary image configurations [800 1280]
seg.BodyCM_estimate(sp,kcam); % calculate the location of body CM and perform a polynomic estimation
seg.trans_frame(fr,sp,kcam,'pltAllFrames',1,'ofst',10,'colormax',[0,73]);axis off % translate the batch images to the location of current frame
colormap('jet'); hc =colorbar('fontsize',8) 
hc.Ticks = [0:20:60,72];
% hold on;plot(seg.bodCM.xy(:,2),seg.bodCM.xy(:,1),'*w');
set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 10 5 4]);
ax = axis;
print(gcf,[figurepath,sparsefile,'\svg\',flname],'-dsvg');
saveas(gcf,[figurepath,sparsefile,'\fig\',flname,'.fig'])
exportgraphics(gcf,[figurepath,sparsefile,'\pdf\',flname,'.pdf']);


end

