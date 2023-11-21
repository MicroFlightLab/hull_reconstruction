function plot_sparse_im_fig1(fr,sp,figurepath,sparsefile,ofset)
%% plot sparse images

for kcam = 1:1:4
    ttl =sprintf('Camera %d frame %d',kcam,fr)
    flname = sprintf('Camera_%d',kcam);
  
imcam{kcam} = Functions.ImfromSp(sp{1}.metaData.frameSize,sp{kcam}.frames(fr).indIm);
imcam_bg = (double(sp{1}.metaData.bg).*(1-(imcam{kcam}>0)*1) + imcam{kcam})/256/256;

  axlims = [min(sp{kcam}.frames(fr).indIm(:,2)),max(sp{kcam}.frames(fr).indIm(:,2))...
                ,min(sp{kcam}.frames(fr).indIm(:,1)),max(sp{kcam}.frames(fr).indIm(:,1))];
            szim = [axlims(2) - axlims(1),axlims(4) - axlims(3)];
            I2 = imcrop(imcam_bg,[axlims(1) - ofset,axlims(3)- ofset,szim(1) + 2*ofset,szim(2) + 2*ofset]);

figure;
imshow(I2,[])  
% title(ttl)
% imwrite(imcam_bg,[figurepath,sparsefile,flname])
print(gcf,[figurepath,sparsefile,'\svg\',flname],'-dsvg')
saveas(gcf,[figurepath,sparsefile,'\fig\',flname,'.fig'])
imwrite(I2,[figurepath,sparsefile,'\png\',flname,'.png'])

end
end

