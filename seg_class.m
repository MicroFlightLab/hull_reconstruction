classdef seg_class<handle
    % seg_class: segmenting the wing and body from the video to later use in hull reconstruction. 
    % input: sp - a cell array of the video without the background
    % 	to segment the fly the user need to run the following functions: 
    %       segcls.flyallImage(sp,kcam);  create a temporary field of the sparse movie in full binary image configurations [800 1280]
    %       segcls.BodyCM_estimate(sp,kcam);  calculate the location of body CM and perform a polynomic estimation 
    %       segcls.transIm_SegWings(sp,kcam);  loop on frames to run translate  image and segment wing.
    % the last 2 functions can run in parallel mode. ( a function doing it
    % is in \hull_reconstruction\+Functions\par_seg ) 
    % when runniing in a non parallel mode: each function has a ploter used for debugging and a progress bar that the user can
    % activate (ploter details inside the functions, progress bar in class initilization). 
    %   plotter for the results: plotSeg(fr,kcam) - plot a segmented frame
   
    
    properties
        st_enfr
        bodCM
        legTH
        delta
        showProgBar
        fit_scope
        cm_poly_degree
        bodyTH
        body
        wing1
        wing2
        all
        frameSize
        tmpImage
        numofwings
        camvec
        openbody
        flythin
    end
    
    methods
        function obj = seg_class(sp,varargin)
            parser = inputParser;
            addParameter(parser,'st_en_fr',0); % start and end frames. if 0 - run all. user input [start end] 
            addParameter(parser,'delta',36); % frames before and after current frame. used to find initial body CM
            addParameter(parser,'showProgBar',0); % show progress bar - for non parallel runs
            addParameter(parser,'fit_scope',100); % scope used for body CM polynom calculation 
            addParameter(parser,'cm_poly_degree',2); % body CM polynom degree 
            addParameter(parser,'bodyTH',0); % threashold used to find body, if 0 - size of scope
            addParameter(parser,'legTH',50); % threashold used to find the wings as the 2 largest blobs (min ammount of pixels)
            addParameter(parser,'frameSize',[800 1280]); % frame size [y x]  
            addParameter(parser,'camvec',[1,2,3,4]); % camera vector  
            addParameter(parser,'openbody',2); % camera vector       
            addParameter(parser,'flythin',1); % make fly thinner by one pixel           
            parse(parser, varargin{:});
            

            
            if length(parser.Results.st_en_fr) == 2
                obj.st_enfr = parser.Results.st_en_fr;
            else
                obj.st_enfr = [1+ parser.Results.delta,length(sp{1}.frames) - parser.Results.delta];
            end
            obj.delta = parser.Results.delta; % ammount of frames to use before and after current frame
            obj.showProgBar = parser.Results.showProgBar; % show/dont show progress bar
            obj.fit_scope = parser.Results.fit_scope; % size of scope to fit CM ploynom
            obj.cm_poly_degree = parser.Results.cm_poly_degree; % degree of CM ploynom
            obj.legTH = parser.Results.legTH;  
            obj.frameSize = parser.Results.frameSize;
            obj.openbody = 2;
            if parser.Results.bodyTH == 0 
                obj.bodyTH = parser.Results.fit_scope;
            else
              obj.bodyTH = parser.Results.bodyTH;
            end
            obj.tmpImage = cell(1,obj.st_enfr(2) - obj.st_enfr(1));
            obj.numofwings = zeros(4,obj.st_enfr(2) - obj.st_enfr(1));
            obj.camvec = parser.Results.camvec;
            % initilize body, wing1 and wing2 cell array 
            obj.body  = cell(1,length(obj.camvec));
            obj.wing1  = cell(1,length(obj.camvec));
            obj.wing2  = cell(1,length(obj.camvec));
            obj.all = cell(1,length(obj.camvec));
            obj.flythin = parser.Results.flythin;
        end
        
        function flyallImage(obj,sp,kcam)
            % create a temporary field of the sparse movie in full binary image configurations [800 1280]
            % apply thin on the image - erasing one pixel 
            for fr = obj.st_enfr(1) : obj.st_enfr(2) + 1
                try
                if isempty(obj.all{kcam}) || isempty(obj.all{kcam}) || fr > length(obj.all{kcam}) || isempty(obj.all{kcam}(fr).indIm)
                   
                    currframe = Functions.ImfromSp([800 1280],sp{kcam}.frames(fr).indIm);
                    currframe = logical(currframe);
                    if obj.flythin == 1
                    currframe = bwmorph(currframe,'thin');
                    end
                    [ro co] = find(currframe);
                    obj.all{kcam}(fr).indIm  = uint16([ro co ones(length(co),1)]);
                    obj.tmpImage{fr} = currframe;
                end
                catch
                    continue
                end
            end
        end
        
        function BodyCM_estimate(obj,sp,kcam,varargin)
            % first pass - estimating body CM by summing 2*delta images and
            % threasholding
            fprintf('BodyCM_estimate \n');
            parser = inputParser;
            addParameter(parser,'plotBod',0);           
            parse(parser, varargin{:});
            
            if obj.showProgBar == 1
                estttl = sprintf('Estimating body CM - camera %d',kcam);
            f = waitbar(0,estttl,'Name',estttl,...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            setappdata(f,'canceling',0);
            end
            
            c = 1;
            obj.bodCM.xy = zeros(obj.st_enfr(2) - obj.st_enfr(1),2);
            for fr = obj.st_enfr(1) : obj.st_enfr(2)
                % loop on frames. for each frame sum +-delta images. find
                % the body by thresholding
                try
                if obj.showProgBar  == 1 && getappdata(f,'canceling')
                    break
                end

                sumim = sumimages(obj,fr,sp,kcam);
                body = sumim(sumim(:,3) == 2*obj.delta+1,:) ;
                
                if parser.Results.plotBod == 1 % plot the estimated body
                    bodimage = Functions.ImfromSp(obj.frameSize,body);
                    imshow(bodimage,[]);
                end
                
                obj.bodCM.xy(c,:) = [mean(body(:,1)),mean((body(:,2)))] ; % xy body location
                if mod(c,100) == 0 && obj.showProgBar == 1
                    waitbar(c/(obj.st_enfr(2)-obj.st_enfr(1)),f,estttl);
                end
                c = c+1;
                catch
                    continue
                end
            end
            obj.bodCM.xy = obj.bodCM.xy(1:c - 1,:); % if initial frame smaller than delta, cut zeros
            obj.st_enfr = [fr - c + 2,obj.st_enfr(1) + c]; % update start and end frame
            
            if obj.showProgBar == 1
            delete(f)
            end
        end
          
        function repPix = sumimages(obj,fr,sp,kcam,varargin)
            % sum images delta before and after fr - used to estimate body
            % CM. fr - frame
            %     sp - sparse images cell
            %     kcam - number of camera
            %     plotSummedIm - plot the summed image
            parser = inputParser;
            addParameter(parser,'plotSummedIm',0);
            parse(parser, varargin{:})
            
            sumim = cell2mat(struct2cell(sp{kcam}.frames((fr-obj.delta):(fr+obj.delta)))');
            repPix = Functions.countRepRow(sumim(:,[1,2]));
            
            if parser.Results.plotSummedIm == 1
                sumImShow = Functions.ImfromSp(obj.frameSize,repPix);
                imshow(sumImShow,[]);
                [ro co] = find(sumImShow>0);
                ofst = 20;
%                 axis([max(min(co)-ofst,0),min(max(co)+ofst,obj.frameSize(2)),...
%                     max(min(ro)-ofst,0),min(max(ro)+ofst,obj.frameSize(1))])
            end
        end
        
        function [drvec,frvec_poly] = calc_dr(obj,fr,varargin)
            % estimate a second (defult) degree polynom for the CM
            % in a scope of 100 (defult) around current frame.
            % calculate the ammount of pixels each image in the scope
            % should move (dr) to be in the same position as current frame
            parser = inputParser;
            addParameter(parser,'pltpolyfit',0);
            parse(parser, varargin{:});
            try
            frvec = [fr - obj.fit_scope : fr + obj.fit_scope];
            polyfr = 1:sum(frvec > obj.st_enfr(1) & frvec < obj.st_enfr(2));
            frvec_poly = frvec(frvec >= obj.st_enfr(1) & frvec < obj.st_enfr(2));
            fr4cm = frvec_poly - obj.st_enfr(1) + 1;
            fr4cm = fr4cm(fr4cm >= 1 & fr4cm <= size(obj.bodCM.xy,1));
            frvec_poly = frvec_poly(1:length(fr4cm));
            
            % estimate polynom
            px = polyfit(frvec_poly, obj.bodCM.xy(fr4cm,1), obj.cm_poly_degree) ;
            py = polyfit(frvec_poly, obj.bodCM.xy(fr4cm,2), obj.cm_poly_degree) ;
            
            if parser.Results.pltpolyfit == 1
                x1 = polyval(px,frvec_poly);
                y1 = polyval(py,frvec_poly);
                hold on
                plot(polyfr,obj.bodCM.xy(polyfr,1));
                plot(polyfr,obj.bodCM.xy(polyfr,2));
                plot(polyfr,x1);
                plot(polyfr,y1);
                figure;plot(obj.bodCM.xy(fr4cm,1),obj.bodCM.xy(fr4cm,2));hold on;plot(x1,y1);
            end
            % calculate CM
            cmvec  = [ polyval(px,frvec_poly)' polyval(py,frvec_poly)' ] ;
            cmcurr = cmvec(polyfr(frvec_poly == fr),:) ;
            % calculate dr
            drvec  = round(cmvec - repmat(cmcurr,size(cmvec,1),1)) ;
            catch
                wakk = 2
            end
        end
        
        function trans_frame(obj,fr,sp,kcam,varargin)
            % substract the fly boundaty (1 extra pixel from sparssing) and
            % translate the fly images in the frames around input fr to the 
            % same location (the translation is calculated using calc_dr function  
            % find body by threasholding (TH can be changed in class
            % initilization: bodyTH)
            
            parser = inputParser;
            addParameter(parser,'pltAllFrames',0);
            addParameter(parser,'colormax',[0,73]);
            addParameter(parser,'ofst',20);
            parse(parser, varargin{:});
            
            [drvec,frvec] = calc_dr(obj,fr); % calc translation
            cnt = 1;
            frvec_delta = fr - obj.delta:fr + obj.delta;
            [frvec in]= intersect(frvec,frvec_delta);
            moved_frame = cell(length(frvec),1);
            drvec = drvec(in,:);
            % move frames to the location of fr
            for frloop = frvec
                drframe = drvec(cnt,:);

                currframe = obj.all{kcam}(frloop).indIm;
                
                moved_frame{cnt} = double(currframe(:,1:2)) - (drframe);
                cnt = cnt + 1;
            end
            
            moved_frame_mat = cell2mat(moved_frame);
            [Result,maph]  = Functions.countRepRow(moved_frame_mat);
            body = uint16(Result(Result(:,3)> obj.bodyTH,:)); % find body by threasholding 
            
            obj.body{kcam}(fr).indIm = uint16(body);
            obj.wing1{kcam}(fr).indIm = uint16(Result(Result(:,3)< obj.bodyTH,:)); 
            bodCM(fr,:) = mean(body(:,1:2));
            if parser.Results.pltAllFrames == 1
                ofst = parser.Results.ofst;
                figure;
                img = Functions.ImfromSp(obj.frameSize,Result);
                imagesc(img,parser.Results.colormax);%colorbar
                [ro co] = find(img>0);axis equal;%hold on;plot(bodCM(:,2),bodCM(:,1),'*r')
                axis([max(min(co)-ofst,0),min(max(co)+ofst,obj.frameSize(2)),...
                    max(min(ro)-ofst,0),min(max(ro)+ofst,obj.frameSize(1))])
                

                hold on;plot(bodCM(fr,2),bodCM(fr,1),'+w','markersize',5,'linewidth',1)
%                 figure;
%                  img = Functions.ImfromSp(obj.frameSize,body);
%                 imshow(img,[]);
%                 
%                 currframe = obj.all{kcam}(fr).indIm;
%                 img2 = Functions.ImfromSp(obj.frameSize,currframe);
%                 imagesc((img2>0)*1+(img>0)*1)
            end
        end
            
        function transIm_SegWings(obj,sp,kcam,varargin)
            % loop on frames to run translte image and segment wing.

            parser = inputParser;
            addParameter(parser,'frmprgbar',500);
            parse(parser, varargin{:});
            
            fprintf('translate image and segment wings \n')
            if obj.showProgBar == 1
                estttl = sprintf('translating fly image and segmenting wing - cam%d',kcam);
                f = waitbar(0,'translating fly image and segmenting wing','Name',estttl,...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                setappdata(f,'canceling',0);
            end
            for fr = obj.st_enfr(1):obj.st_enfr(2)
                try
                    if mod(fr,parser.Results.frmprgbar) == 0 && obj.showProgBar == 1
                        waitbar((fr - obj.st_enfr(1))/(obj.st_enfr(2)-obj.st_enfr(1)),f,estttl);
                    end
                    if obj.showProgBar  == 1 && getappdata(f,'canceling')
                        break
                    end
                    % this is the interesting part of the loop: 
                    obj.trans_frame(fr,sp,kcam); % translate the batch images to the location of current frame
                    obj.segmentwings(fr,kcam);  % segment wings
                    
                catch
                    continue
                end
            end
            if obj.showProgBar == 1
                delete(f)
            end
        end

        function segmentwings(obj,fr,kcam,varargin)
            % segment wings after body segmentation. substract the image
            % with a dilated body image. define a blob as wing only if its
            % larger than legTH
            parser = inputParser;
            addParameter(parser,'plotwings',0);
            parse(parser, varargin{:});
            
            SE = strel('disk',2);
            wingname = {'wing1','wing2'};
            imgBody = Functions.ImfromSp(obj.frameSize,obj.body{kcam}(fr).indIm);
            imgall = obj.tmpImage{fr};%Functions.ImfromSp(obj.frameSize,obj.all{kcam}(fr).indIm);%
            imgBody = imdilate(imgBody,SE);
            
            imgdiff = ((imgall>0)*1-(imgBody>0)*1)>0;
            imgdiff = imclose(imgdiff,SE);
            propsIm = bwconncomp(imgdiff);
            [val,ind] = sort(cellfun(@length, propsIm.PixelIdxList),'descend');
            
            ind(val<obj.legTH) = [];
            if length(ind) == 1
                ind = [ind ind];
            end
            obj.numofwings(kcam,fr) = length(ind);
            for k = 1:1:length(wingname)
                [ro co] = ind2sub(obj.frameSize, propsIm.PixelIdxList{ind(k)});
                obj.(wingname{k}){kcam}(fr).indIm = uint16([ro co ones(length(ro),1)]);
            end
            if parser.Results.plotwings == 1
                imgallWing1 = Functions.ImfromSp(obj.frameSize,obj.wing1{kcam}(fr).indIm);
                imgallWing2 = Functions.ImfromSp(obj.frameSize,obj.wing2{kcam}(fr).indIm);
                imagesc(imgallWing1+2*imgallWing2);
            end
        end
        
        function prep4save(obj,body,wing1,wing2,image,frames)
            % arange output from parallel run in prepro class
            for k = obj.camvec
                obj.body{k} = body{k};
                obj.wing1{k} = wing1{k};
                obj.wing2{k} = wing2{k};
                obj.all{k} = image{k};
                obj.st_enfr = frames(1,:);
            end
            obj.tmpImage = []
            
        end
        
        
        function plotSeg(obj,fr,kcam,varargin)
            parser = inputParser;
            addParameter(parser,'body',1);
            addParameter(parser,'wings',1);
            addParameter(parser,'ofset',0);
            parse(parser, varargin{:});
            ofset = parser.Results.ofset;
            % plot the segmented fly
            allIm = Functions.ImfromSp(obj.frameSize,obj.all{kcam}(fr).indIm);
            imshow(allIm,[]);
            hold on
            if parser.Results.body == 1
            scatter(obj.body{kcam}(fr).indIm(:,2),obj.body{kcam}(fr).indIm(:,1),'g');
            end
            if parser.Results.wings == 1
            scatter(obj.wing1{kcam}(fr).indIm(:,2),obj.wing1{kcam}(fr).indIm(:,1),'r');
            scatter(obj.wing2{kcam}(fr).indIm(:,2),obj.wing2{kcam}(fr).indIm(:,1),'b');
            end
            axlims = [min(obj.all{kcam}(fr).indIm(:,2))-ofset,max(obj.all{kcam}(fr).indIm(:,2))+ofset...
                ,min(obj.all{kcam}(fr).indIm(:,1))-ofset,max(obj.all{kcam}(fr).indIm(:,1))+ofset];
            axis(axlims)
        end
        
        
            
    end
    
    
end

