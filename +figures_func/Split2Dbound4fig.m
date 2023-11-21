function [wingIm,DownIm,UpIm,span2D] = Split2Dbound4fig(hull,TwoD,SE,kcam)
            % Generate the boundary of the projected hull and split it.
            % make sure there are no holes in the projected image by
            % closeing it (SE)
            % TwoD contains the projection of LE,TE,wing_cut,Tip,Wing
            % CM,body.
            ypix = 0;
            mult = -1;
            if kcam == hull.cameras.all.ZaxCam
                ypix = 801
                mult = 1;
            end
            span2D=diff([TwoD{5};TwoD{4}])/norm([TwoD{5};TwoD{4}]); % calculate the 2D span
            slp=1./(span2D/span2D(2));
            b=TwoD{5}(2)-TwoD{5}(1)*slp(1);
            wingAll=[TwoD{1};TwoD{2}];
            if isinf(slp(1))
                yperp_wall=(ones(1,length(wingAll(:,1)))*TwoD{5}(2))';
            else
                yperp_wall=slp(1)*wingAll(:,1)+b;
            end
            % split the image according to span
            Down=wingAll(wingAll(:,2)<yperp_wall,:);
            Down(:,2) = ypix - Down(:,2)*mult;
            
            Up=wingAll(wingAll(:,2)>=yperp_wall,:);
            Up(:,2) = ypix - Up(:,2)*mult;
            TwoD{3}(:,2) = ypix - TwoD{3}(:,2)*mult;
            [wingIm] = Functions.ImfromSp(hull.cameras.all.size_image,fliplr(TwoD{3}));
            wingIm=imclose(wingIm,SE);
            
            
            [DownIm] = Functions.ImfromSp(hull.cameras.all.size_image,fliplr(Down));
            DownIm=imclose(DownIm,SE);
            
            [UpIm] = Functions.ImfromSp(hull.cameras.all.size_image,fliplr(Up));
            UpIm=imclose(UpIm,SE);
end
        
