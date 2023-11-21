function plot_all_all_wakk_fig2(seg,fr,figsdir,figurepath,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
parser = inputParser;
addParameter(parser,'camvec',[1,2,3,4]);
parse(parser, varargin{:});
figure;
ofset = 30
for kcam = parser.Results.camvec
    allIm = Functions.ImfromSp(seg.frameSize,seg.all{kcam}(fr).indIm);
    bodyIm = Functions.ImfromSp(seg.frameSize,seg.body{kcam}(fr).indIm);
    allIm = allIm.*(1-bodyIm);
    
    rgbimage_tmp = cat(3, bodyIm*(0), bodyIm*(1), bodyIm*(0));
    rgbimage = rgbimage_tmp + allIm;
    subplot(2,4,kcam);h = imshow(rgbimage);hold on;
    axlims = [min(seg.all{kcam}(fr).indIm(:,2))-ofset,max(seg.all{kcam}(fr).indIm(:,2))+ofset...
        ,min(seg.all{kcam}(fr).indIm(:,1))-ofset,max(seg.all{kcam}(fr).indIm(:,1))+ofset];
    axis(axlims)
    ttl = sprintf('Camera %d',kcam);
    title(ttl,'fontsize',10,'fontweight','normal')
end


for kcam = parser.Results.camvec
    ttl = sprintf('Camera %d',kcam);
    rgbimage_tmp_wing1 = zeros([seg.frameSize,3]);   
    rgbimage_tmp_wing2 = zeros([seg.frameSize,3]);
    wingIm1 = zeros(seg.frameSize);
    wingIm2 = zeros(seg.frameSize);
      bodyIm = zeros(seg.frameSize);
    rgbimage_tmp  = zeros([seg.frameSize,3]);
  
    if kcam == 1
        wingIm1 = Functions.ImfromSp(seg.frameSize,seg.wing1{kcam}(fr).indIm);
        wingIm2 = Functions.ImfromSp(seg.frameSize,seg.wing2{kcam}(fr).indIm);
        rgbimage_tmp_wing1 = cat(3, wingIm1*(1), wingIm1*(0), wingIm1*(0));
        rgbimage_tmp_wing2 = cat(3, wingIm2*(0), wingIm2*(0), wingIm2*(1));
    bodyIm = Functions.ImfromSp(seg.frameSize,seg.body{kcam}(fr).indIm);
    rgbimage_tmp = cat(3, bodyIm*(0), bodyIm*(1), bodyIm*(0));

    end
    
    allIm = Functions.ImfromSp(seg.frameSize,seg.all{kcam}(fr).indIm);
    
    allIm = allIm.*(1-(bodyIm + wingIm2 + wingIm1));
    
    
    rgbimage = rgbimage_tmp + allIm + rgbimage_tmp_wing1 + rgbimage_tmp_wing2;
    subplot(2,4,kcam + 4);h = imshow(rgbimage);hold on;
    axlims = [min(seg.all{kcam}(fr).indIm(:,2))-ofset,max(seg.all{kcam}(fr).indIm(:,2))+ofset...
        ,min(seg.all{kcam}(fr).indIm(:,1))-ofset,max(seg.all{kcam}(fr).indIm(:,1))+ofset];
    axis(axlims)
    title(ttl,'fontsize',10,'fontweight','normal')
end
flname = ['\hull 2d all all all_body'];
print(gcf,[figurepath,figsdir,'\svg\',flname],'-dsvg')
saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig'])
set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 5 12 10]);
exportgraphics(gcf,[figurepath,figsdir,'\pdf\',flname,'.pdf']);

end

