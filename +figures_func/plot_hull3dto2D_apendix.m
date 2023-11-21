function plot_hull3dto2D_apendix(fr,hull,loaders,figsdir,flname,figurepath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
frm = find(fr == hull.frames);

easy = loaders.easywand();
seg = loaders.loadSegfile(1);
parts2run = {'all','all','all','all'}

hullRec = hullrec_class(easy,'ofst',5,'VxlSize4search',20e-5,'TseviMethod',0);

im4hull = {seg.all{1}(frm).indIm,seg.all{2}(frm).indIm,seg.all{3}(frm).indIm,seg.all{4}(frm).indIm};

figure;
[hull_recon] = Functions.runOnlyHull(im4hull,parts2run,hullRec,'all');
[hullRealR,rtmat] = Functions.hullRec2lab(hull_recon,hullRec.cameras.all.Rotation_Matrix,hullRec.cameras.all.RotMat_vol,hullRec.real_3d);


print(gcf,[figurepath,figsdir,'\svg\',flname],'-dsvg')
saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig'])

figure;
plot3(hullRealR(:,1),hullRealR(:,2),hullRealR(:,3),'.b');title('Ristroph');axis equal;grid on;box on;xlabel('mm');ylabel('mm');zlabel('mm')

end

