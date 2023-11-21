function plot2DBound_fig2(sp,ploters,fr,hull,hull3d,figurepath,figsdir,flname,kcam,varargin)

parser = inputParser;
            addParameter(parser,'savefig',1); %
            addParameter(parser,'only2d',0); %
            addParameter(parser,'nobg',0); %

            parse(parser, varargin{:})




frm = find(fr == hull.frames);


realC = hull.real_coord{frm};

allIm = Functions.ImfromSp(sp{kcam}.metaData.frameSize,sp{kcam}.frames(fr).indIm);
rgbimage = zeros(hull.cameras.all.size_image);
wingcell = {'rightwing','leftwing'};

realC = hull.real_coord{frm};

rgbimage = zeros(hull.cameras.all.size_image);SE = strel('disk',3);
for kwing =1:1:2
    wingname = wingcell{kwing};
    hull.body.hull3d = hull3d.body.body4plot{frm};
    hull.(wingname).hull3d = hull3d.(wingname).hull.hull3d{frm};
    [LE ,TE] = figures_func.wing2LETE(hull3d,hull,wingname,fr);
    
    TwoD = Functions.Create2DCoords(realC,hull.cameras.all.DLT_coefs(:,kcam),{LE...
        ,TE,hull.(wingname).hull3d,hull.(wingname).coords.tip(frm,:),...
        mean(hull.(wingname).hull3d),hull.body.hull3d},1,hull.cameras.all.size_image(1),hull.cameras.all.RotMat_vol);
    [wingIm,DownIm,UpIm,spn] = figures_func.Split2Dbound4fig(hull,TwoD,SE,kcam);
     spnW{kwing} = spn;
    if kcam == hull.cameras.all.ZaxCam;
        TwoD{6}(:,2) = 801 - TwoD{6}(:,2);
        TwoD{5}(:,2) = 801 - TwoD{5}(:,2);
         spnW{kwing} = [spn(1), -spn(2)];
    end
   
    XYspnW{kwing} = TwoD{5};
    [body] = Functions.ImfromSp(hull.cameras.all.size_image,fliplr(TwoD{6}));
    
    bounWing_dil = bwperim(wingIm);
    bounWing_dil = imdilate(bounWing_dil,hull.wingana.bound_dilate);
    bounWing{1}=bounWing_dil.*DownIm;
    bounWing{2}=bounWing_dil.*UpIm;
    
    [ro co ] = find(bounWing{1});
    [ro2 co2 ] = find(bounWing{2});
        [rob cob ] = find(body);

    cm_2d = [mean(rob),mean(cob)];


    coord_B{kwing} = [ro co ];
    coord_B2{kwing} = [ro2 co2 ];
    
    matcol = {[1,0,0;1,0,1;0 1 0;0,1,0],[0.1,0.1,1;0,1,1;0 1 0;0,1,0]};
    
    
    k = 1;
    rgbimage_tmp = cat(3, DownIm*(matcol{kwing}(k,1)), DownIm*(matcol{kwing}(k,2)), DownIm*(matcol{kwing}(k,3)))...
        +cat(3, UpIm*(matcol{kwing}(k+1,1)), UpIm*(matcol{kwing}(k+1,2)), UpIm*(matcol{kwing}(k+1,3)))...
        +cat(3, body*(matcol{kwing}(k+2,1)), body*(matcol{kwing}(k+2,2)), body*(matcol{kwing}(k+2,3)));
    rgbimage = rgbimage + rgbimage_tmp;
    
end
axlims = ploters.plotSpIm(sp,fr,kcam,'inscm',body);
alpha  = 0.6;

if parser.Results.savefig == 0
% hold on;scatter(coord_B{1}(:,2),coord_B{1}(:,1),'o','filled','MarkerEdgeColor',matcol{1}(1,:),'MarkerFaceColor',matcol{1}(1,:));
% hold on;scatter(coord_B2{1}(:,2),coord_B2{1}(:,1),'o','filled','MarkerEdgeColor',matcol{1}(2,:),'MarkerFaceColor',matcol{1}(2,:));
% hold on;scatter(coord_B{2}(:,2),coord_B{2}(:,1),'o','filled','MarkerEdgeColor',matcol{2}(1,:),'MarkerFaceColor',matcol{2}(1,:));
% hold on;scatter(coord_B2{2}(:,2),coord_B2{2}(:,1),'o','filled','MarkerEdgeColor',matcol{2}(2,:),'MarkerFaceColor',matcol{2}(2,:));
if parser.Results.nobg == 0    
bgnofly = double(sp{kcam}.metaData.bg)/double(max(sp{kcam}.metaData.bg(:))).*(1-allIm>0*1);
img = rgbimage+bgnofly+(double(allIm)/max(double(allIm(:))));
else
    img = rgbimage;
end
if parser.Results.only2d == 0
h = imshow(img);hold on;
set(h, 'AlphaData', alpha);
end



if parser.Results.only2d == 0
hold on;quiver(XYspnW{1}(1),XYspnW{1}(2),spnW{1}(1),spnW{1}(2),1200,'k','LineWidth',3);
hold on;quiver(XYspnW{2}(1),XYspnW{2}(2),spnW{2}(1),spnW{2}(2),1200,'k','LineWidth',3);
end
end
legend off

if parser.Results.savefig == 1
    bgnofly = double(sp{kcam}.metaData.bg)/double(max(sp{kcam}.metaData.bg(:))).*(1-allIm>0*1);
h = imshow(rgbimage+bgnofly+(double(allIm)/max(double(allIm(:)))));hold on;

set(h, 'AlphaData', alpha);

hold on;scatter(coord_B{1}(:,2),coord_B{1}(:,1),'o','filled','MarkerEdgeColor','k','MarkerFaceColor',matcol{1}(1,:));
hold on;scatter(coord_B2{1}(:,2),coord_B2{1}(:,1),'o','filled','MarkerEdgeColor','k','MarkerFaceColor',matcol{1}(2,:));
hold on;scatter(coord_B{2}(:,2),coord_B{2}(:,1),'o','filled','MarkerEdgeColor','k','MarkerFaceColor',matcol{2}(1,:));
hold on;scatter(coord_B2{2}(:,2),coord_B2{2}(:,1),'o','filled','MarkerEdgeColor','k','MarkerFaceColor',matcol{2}(2,:));
hold on;quiver(XYspnW{1}(1),XYspnW{1}(2),spnW{1}(1),spnW{1}(2),1200,'k','LineWidth',3);
hold on;quiver(XYspnW{2}(1),XYspnW{2}(2),spnW{2}(1),spnW{2}(2),1200,'k','LineWidth',3);


    
set(gcf,'inverthardcopy','off','color','w','paperpositionmode','auto','units','centimeters'...
    ,'position',[10 10 10 8]);
set(gcf,'renderer','painters');

print(gcf,[figurepath,figsdir,'\svg\',flname],'-dsvg');
saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig']);
end

end

