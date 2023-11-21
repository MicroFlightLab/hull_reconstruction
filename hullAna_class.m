classdef hullAna_class<handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        body
        rightwing
        leftwing
        bodyana
        wingana
        cameras
        real_coord
        rotmat_EWtoL
        frames
        sprs
        video
    end
    
    methods
        function obj = hullAna_class(hullRec,varargin)
            parser = inputParser;
            addParameter(parser,'bodyfrac',0.3); %  defines the size of the sphere around the CM - leaves only the heat and tail to choose refine the outcome of PCA
            addParameter(parser,'distfromCM',5); %
            addParameter(parser,'percWing_distfromBod',[0.4 0.65]);
            addParameter(parser,'DeltaTubeMax',10);
            addParameter(parser,'DeltaTubeMin',0.6);
            addParameter(parser,'per2cutWing',0.1);
            addParameter(parser,'angleTH',140);
            addParameter(parser,'deltaChordSlice',1.1);
            addParameter(parser,'image2D_close',3);
            addParameter(parser,'bound_dilate',2);
            addParameter(parser,'wingSize',43);
            addParameter(parser,'per2cutWing4CM',0.3);
            addParameter(parser,'wingfract_tip',0.05);
            addParameter(parser,'part4TE',2);
            addParameter(parser,'part4LE',6);
            addParameter(parser,'numofsec',4);
            addParameter(parser,'coneang',15);
            parse(parser, varargin{:});
            
            obj.wingana.coneang = parser.Results.coneang;
            obj.wingana.part4TE = parser.Results.part4TE;
            obj.wingana.part4LE = parser.Results.part4LE;
            obj.wingana.numofsec = parser.Results.numofsec;
            obj.wingana.per2cutWing4CM = parser.Results.per2cutWing4CM;
            obj.wingana.image2D_close =strel('disk',parser.Results.image2D_close);
            obj.wingana.bound_dilate =strel('disk',parser.Results.bound_dilate);
            obj.wingana.per2cutWing = parser.Results.per2cutWing;
            obj.wingana.angleTH = parser.Results.angleTH;
            obj.bodyana.frac = parser.Results.bodyfrac;
            obj.bodyana.distfromCM = parser.Results.distfromCM;
            obj.wingana.wingfract_tip = parser.Results.wingfract_tip;
            obj.wingana.DeltaTubeMax = parser.Results.DeltaTubeMax;
            obj.wingana.DeltaTubeMin = parser.Results.DeltaTubeMin;
            obj.wingana.deltaChordSlice = parser.Results.deltaChordSlice;
            obj.wingana.wingSize = parser.Results.wingSize;
            obj.wingana.percWing_distfromBod = parser.Results.percWing_distfromBod;
            
            obj.cameras = hullRec.cameras;
            obj.rotmat_EWtoL = hullRec.cameras.all.Rotation_Matrix*hullRec.cameras.all.RotMat_vol';
        end
        
        function [head_tailCM,headVoxels,tailVoxels] = Xaxis(obj,realC,varargin)
            parser = inputParser;
            addParameter(parser,'plot',0); %
            addParameter(parser,'prop','body'); %
            addParameter(parser,'savename','X'); %
            addParameter(parser,'save2hull',1)
            parse(parser, varargin{:})
            
            prop = parser.Results.prop;
            savename = parser.Results.savename;
            
            hullReal = (obj.rotmat_EWtoL * double(obj.(prop).hull3d'))';
            
            CM = mean(hullReal);
            
            p = pca(double(hullReal)) ; % find AHat from pca
            AHat  = sign(p(3,1))*p(:,1)' ;
            Nvox = ceil(size(hullReal,1)*obj.bodyana.frac);
            
            % find distances of body voxels from bodyCM
            distFromRcm = ( sum ((hullReal - CM).^2,2) ) .^ 0.5 ;
            [~, sortedInd] = sort(distFromRcm,'descend') ;
            
            selectedInd    = sortedInd(1:Nvox) ;
            
            % separate the voxels that are at the head and tail of the fly. use AHat.
            % if AHat is not avilable then an esitmation for the axis could be found by
            % taking the most distant pair of voxels.
            
            mat1 = double(hullReal(selectedInd,:)) - repmat(CM, Nvox, 1) ;
            mat2 = repmat (AHat, Nvox, 1) ;
            dotprod = sum( mat1 .* mat2 , 2) ;
            
            % find the voxels with distance around the maximum/minimum
            maxdist = max(dotprod) ;
            mindist = min(dotprod) ;
            deltahead = (maxdist-mindist)/obj.bodyana.distfromCM; % threshold for distance from cm
            deltarail = (maxdist-mindist)/obj.bodyana.distfromCM; % threshold for distance from cm
            headFlag = (dotprod>(maxdist-deltahead)) ;
            tailFlag = (dotprod<(mindist+deltarail)) ;
            
            
            % check if the head voxels are in one or more connected component
            % if there's more than one, take the largest
            headVoxels = obj.(prop).hull3d(selectedInd(headFlag),:) ;
            headVoxels =  Functions.findLargestHullCC (int32(headVoxels)) ;
            
            tailVoxels = obj.(prop).hull3d(selectedInd(tailFlag),:) ;
            tailVoxels =  Functions.findLargestHullCC(int32(tailVoxels)) ;
            
            head_tailCM = [mean(headVoxels);mean(tailVoxels)] ;
            
            
            
            if strcmp(prop,'body') == 0
                bodLab = (obj.rotmat_EWtoL * double(obj.('body').hull3d)')';
                cmbody = mean(bodLab);
                
                tlHdDistCM = sqrt(sum((cmbody - head_tailCM).^2,2));
                if tlHdDistCM(1) < tlHdDistCM(2)
                    head_tailCM = flipud(head_tailCM);
                end
                rtmat = eye(3);
            end
            
            
            Xax = (obj.rotmat_EWtoL *(head_tailCM(1,:) - head_tailCM(2,:))')' ;
            if parser.Results.save2hull == 1
                obj.(prop).vectors.(savename) = ((Xax / norm(Xax))')' ;
                
                
                if strcmp(prop,'body') == 1
                    headreal = Functions.hullRec2lab(headVoxels,obj.cameras.all.Rotation_Matrix,...
                        obj.cameras.all.RotMat_vol,realC);
                    
                    tailreal = Functions.hullRec2lab(tailVoxels,obj.cameras.all.Rotation_Matrix,...
                        obj.cameras.all.RotMat_vol,realC);
                    obj.(prop).coords.CM_real = (mean(headreal)+mean(tailreal))/2;
                end
                
                
                %             obj.(prop).coords.CM = ((head_tailCM(1,:) + head_tailCM(2,:))/2)';
                
            end
            
            if parser.Results.plot == 1
                %                 pcshow(hullReal,'r');
                %                 hold on;quiver3(CM(1),CM(2),CM(3),obj.(prop).vectors.(savename)(1),obj.(prop).vectors.(savename)(2),obj.(prop).vectors.(savename)(3),50);
                cmtail = mean(tailVoxels);
                plot3(tailVoxels(:,1),tailVoxels(:,2),tailVoxels(:,3),'.g');hold on;
                plot3(headVoxels(:,1),headVoxels(:,2),headVoxels(:,3),'.g');
                hold on;quiver3(cmtail(1),cmtail(2),cmtail(3),obj.(prop).vectors.(savename)(1),obj.(prop).vectors.(savename)(2),obj.(prop).vectors.(savename)(3),50,'color','k','linewidth',3);
                
            end
        end
        
        function splitwings(obj,Twowings,varargin)
            parser = inputParser;
            addParameter(parser,'name2save','span'); %
            addParameter(parser,'plotkmeans',0); %
            parse(parser, varargin{:})
            
            %-----------identify wing1 and wing2. if its one blob split it----
            Two_wings = double(Twowings);
            [idx,C] = kmeans(Two_wings,4,'Distance','sqeuclidean');
            wings_coor_cell={Two_wings(idx==1,:),Two_wings(idx==2,:),Two_wings(idx==3,:),Two_wings(idx==4,:)};
            
            ind_del=0;ind2del=[];
            wings_coor_cell_split=[];
            for k=1:1:4
                [mat3d] = Functions.ptcloud2binMat([wings_coor_cell{k}]);
                CC = bwconncomp(mat3d); % count if the two wings have 1 or 2 clusters to identify the case
                if CC.NumObjects>1
                    ind_del=ind_del+1;ind2del(ind_del)=k;
                    for kobj=1:1:CC.NumObjects
                        blob=CC.PixelIdxList{kobj};
                        [wings_all1,wings_all2 , wings_all3]=ind2sub([CC.ImageSize],blob);
                        wings_coor_cell_split{1+length(wings_coor_cell_split)}=[wings_all1,wings_all2 , wings_all3];
                    end
                end
            end
            if ind_del>0
                wings_coor_cell(ind2del)=[];
                wings_coor_cell(end+1:end+length(wings_coor_cell_split))=wings_coor_cell_split;
            end
            max_dist_wings=zeros(1,length(wings_coor_cell));
            for k=1:1:length(wings_coor_cell)
                max_dist_wings(k)=Functions.find_max_min_dist(wings_coor_cell{k},mean(obj.body.hull3d));
            end
            [a,indDist]=sort(max_dist_wings);
            % Calculate the distance between blobs that were found by kmeans
            % and identify the 2 blobs (each one is a wing)
            notFarBlobs=cell2mat(wings_coor_cell(indDist(1:end-2))');
            dist_1=sqrt(sum((notFarBlobs-mean(wings_coor_cell{indDist(end)})).^2,2));
            dist_2=sqrt(sum((notFarBlobs-mean(wings_coor_cell{indDist(end-1)})).^2,2));
            
            closer= dist_1<dist_2 ;
            wing1=[notFarBlobs(closer,:);wings_coor_cell{indDist(end)}];
            wing2=[notFarBlobs((1-closer)>0,:);wings_coor_cell{indDist(end-1)}];
            
            obj.rightwing.hull3d = Functions.findLargestHullCC(wing1);
            obj.leftwing.hull3d = Functions.findLargestHullCC(wing2);
            
            
            if isempty(parser.Results.plotkmeans) == 1
                km1 = (obj.rotmat_EWtoL*notFarBlobs(closer,:)')';
                km2 = (obj.rotmat_EWtoL*wings_coor_cell{indDist(end)}')';
                km3 = (obj.rotmat_EWtoL*notFarBlobs((1-closer)>0,:)')';
                km4 = (obj.rotmat_EWtoL*wings_coor_cell{indDist(end-1)}')';
                bod = (obj.rotmat_EWtoL*double(parser.Results.plotkmeans)')';
                plot3(km1(:,1),km1(:,2),km1(:,3),'.');hold on;
                plot3(km2(:,1),km2(:,2),km2(:,3),'.');
                plot3(km3(:,1),km3(:,2),km3(:,3),'.');
                plot3(km4(:,1),km4(:,2),km4(:,3),'.');
                plot3(bod(:,1),bod(:,2),bod(:,3),'g.');
            end
            
            
        end
        
        function [spanVec,wingCM,wingLargestCC,forCM] = first_est_span(obj,prop)

            % take the farthest "fraction" of the voxels and calculate
            % their mean position - this is the estimated wing tip.
            % recalculate the span using this tip
            % find distances of coords from wingCM
            wingSize = obj.wingana.wingSize; % 41
            percWing_distfromBod = obj.wingana.percWing_distfromBod;%[ 0.4 0.65];

            wing = (double(obj.(prop).hull3d)')';
            distfromBod = wingSize * percWing_distfromBod; %
            cBody = mean(obj.body.hull3d) ;
            wingLargestCC_initial = Functions.findLargestHullCC(wing);
            [farPoint, ~, list] = Functions.farthestPoint(wingLargestCC_initial, cBody, wingSize) ;
            wingLargestCC = Functions.findLargestHullCC (wingLargestCC_initial(list,:));

            [~, ~, farfrombod] = Functions.farthestPoint(wingLargestCC, cBody, distfromBod(2)) ; % list2 used for finding wing cm later
            [~, ~, closetobod] = Functions.farthestPoint(wingLargestCC, cBody, distfromBod(1)) ; % list2 used for finding wing cm later
            % cut a percentage from the wing's base and keep the largest
            % hull

            % etimate wing center of mass from a strip located aproximatly
            % at the center of the wing
            wingCoordsForCM_1 = wingLargestCC(closetobod,:) ;
            wingCoordsForCM_2 = wingLargestCC(farfrombod,:) ;

            forCM = (1 - ismember(wingCoordsForCM_2,wingCoordsForCM_1,'rows')) == 1;
            forCM = wingCoordsForCM_2(forCM,:);
            wingCM = mean(forCM) ;

            % calculate span vector with respect to wing cm
            % later will also calculate leading edge vector
            spanHat = (double(farPoint) - wingCM)' ;
            spanVec = spanHat / norm(spanHat) ;
        end

        function [spanHat,cm_coord4tip,coords_cone] = est_tip_recalc_span(obj,coords,wingCM,spanVec,forCM,prop,name2save)
            coneang = obj.wingana.coneang;
           
            Nvox = size(coords,1) ;
            % find tip location
            % 1st screening
            % -------------
            % consider only voxels in a cone of obj.wingana.coneang degrees around the line starting
            % from wingCM along the direction of spanVec
            mat1 = double(coords) - repmat(wingCM, Nvox, 1) ;
            mat1 = mat1 ./ repmat( vecnorm(mat1,2,2), 1,3) ;
            mat2 = repmat(spanVec', Nvox, 1) ;
            dotprod = dot(mat1, mat2, 2) ;
            
            ind1 = dotprod > cosd(coneang) ; 
            
            % 2nd screening
            % -------------
            % take the farthest "fraction" of the voxels and calculate their mean position

            coords_cone = double(coords(ind1,:)) ;
            [spanHat,cm_coord4tip] =span_from_tip_CM(obj,coords_cone,forCM,prop,name2save);

%             Nvox   = size(coords_cone,1) ;
%             
%             dst  = vecnorm(coords_cone - repmat(wingCM,Nvox,1),2,2) ;
%             
%             [~, sortedInd] = sort(dst,'descend') ;
%             Nvox1 = ceil(Nvox*fraction) ;
%             selectedInd    = sortedInd(1:Nvox1) ;
%             
%             cm_coord4tip = coords_cone(selectedInd,:) ;
%             
%             % recalculate span
%             if (isempty(selectedInd)) == 0
%                 spanHat = (mean(cm_coord4tip,1) - wingCM) ;
%                 spanHat = (spanHat / norm(spanHat))'';
%             end
        end

        function [TE_chrd,LE_chrd] = est_chord(obj,forCM,prop)
            % estimate the chord, take the wing strip arounde CM wing and
            % kmeans to 2 groups. define as LE the group that is above the projection of wingCM on Xbody
       
            [idx,C] = kmeans(double(forCM),2,'Distance','sqeuclidean');
            grp_chrd={forCM(idx==1,:),forCM(idx==2,:)};
            mean_chrds_grpLoc=dot(obj.body.vectors.X,mean(forCM));
            dist_headChord=mean_chrds_grpLoc - [ dot(mean(grp_chrd{1}),obj.body.vectors.X), dot(mean(grp_chrd{2},1),obj.body.vectors.X)];

            [~,ind_chor]=min(dist_headChord);
            TE_chrd=grp_chrd{ind_chor(1)};
            LE_chrd=grp_chrd{(ind_chor(1)==2)*1+(ind_chor(1)==1)*2};

            obj.(prop).vectors.chord=(mean(TE_chrd)-mean(LE_chrd))/norm(mean(TE_chrd)-mean(LE_chrd));


        end

        function [forCM_fixed] = relocateCM(obj,forCM,wingLargestCC,spanHat,prop,cm_coord4tip)
            % Perform dot product of wingCM with the span. devide wingCM to
            % 4 stripes, calculate the projection of each strip on the chord
            % use as reference the stripe that is closest to the initial
            % guess of CM. then, devide the whole wing to 8 stripes. 
            % check the difference in the amount of voxels
            % and the length of the stripe projected on the chord 
            % choose a stripe that is similiar to the referance but located
            % the furthest from the tip
            
            crd = obj.(prop).vectors.chord;
            % count number of voxels for each strip:
            [amnt,rw] = getamount_stripe(obj,forCM,spanHat,4);
            [vl in_v ] = max(amnt);
            % use the one with maximal aount of voxels to calculate the dot 
            % product on chord and to get its size, this is our referance
            maxst = forCM(rw{in_v},:);
            distn_strp = (dot(maxst',repmat(crd',1,size(maxst,1))));
            maxpt_ref = max(distn_strp);
            minpt_ref = min(distn_strp);
            mean_CM_ref = mean(forCM(rw{in_v},:));
            
            % devide the whole wing to 8 stripes, count the number of
            % voxels for each strip
            [amnt_allwing,rw_all_wing] = getamount_stripe(obj,wingLargestCC,spanHat,8);
            
            % find the referance stripe and calculate the size on cord of
            % the others
            diff_meanCM_old= inf;
            for k = 1:1:length(rw_all_wing)
                strip = wingLargestCC(rw_all_wing{k},:);
                diff_meanCM = sum((mean(strip) - mean_CM_ref).^2);
                if diff_meanCM<diff_meanCM_old
                    rws_ref = size(strip,1);
                    diff_meanCM_old = diff_meanCM;
                end
                
                distn_strp = (dot(strip',repmat(crd',1,size(strip,1))));
                maxpt(k) = max(distn_strp);
                minpt(k) = min(distn_strp);
                rws_all(k) = length(rw_all_wing{k});
            end
            
            % use the stripe that is located the furthest from the tip and
            % is about the same size of the referance strip
            maxdiff = maxpt - maxpt_ref;
            mindiff = minpt - minpt_ref;
            
            strips_wng = find(abs(maxdiff)< 2 & abs(mindiff)< 2 & (rws_ref-rws_all)<100);
            cms = wingLargestCC(cell2mat(rw_all_wing(strips_wng)),:);
            
            
            for k = 1:1:length(strips_wng)
                strip = wingLargestCC(rw_all_wing{strips_wng(k)},:);
                mean_dist_cm(k) = sum(mean(cm_coord4tip,1) - mean(strip)).^2;
            end
            if length(strips_wng) == 0
                forCM_fixed = maxst;
            else
                
            [v ij] = max(mean_dist_cm);
            maxind = max(length(ij));
            forCM_fixed = wingLargestCC(rw_all_wing{strips_wng(maxind)},:);
            end
        end

        function [spanHat,cm_coord4tip,new_wCM ]= span_from_tip_CM(obj,coords_cone,forCM_fixed,prop,name2save)
             fraction = obj.wingana.wingfract_tip;
            new_wCM = mean(forCM_fixed);
            Nvox   = size(coords_cone,1) ;
            
            dst  = vecnorm(coords_cone - repmat(new_wCM,Nvox,1),2,2) ;
            
            [~, sortedInd] = sort(dst,'descend') ;
            Nvox1 = ceil(Nvox*fraction) ;
            selectedInd    = sortedInd(1:Nvox1) ;
            cm_coord4tip = coords_cone(selectedInd,:) ;
            
            % recalculate span
            if (isempty(selectedInd)) == 0
                spanHat = (mean(cm_coord4tip,1) - new_wCM) ;
                spanHat=(spanHat / norm(spanHat))';
                obj.(prop).vectors.(name2save) = spanHat';
            end
            obj.(prop).vectors.wingCM = forCM_fixed;
            obj.(prop).vectors.alltip = cm_coord4tip;
        end


        function [cm_coord,new_wCM,TE_chrd,LE_chrd] = EstTip_calcSpanChord(obj,prop,varargin)
            parser = inputParser;
            addParameter(parser,'name2save','span'); %
            addParameter(parser,'plotChord',0); %plotcone
            addParameter(parser,'plotcone',0);
            parse(parser, varargin{:})
            name2save = parser.Results.name2save;
            

            % take the farthest "fraction" of the voxels and calculate
            % their mean position - this is the estimated wing tip.
            % recalculate the span using this tip
            % find distances of coords from wingCM
            [spanVec,wingCM,coords,forCM] = first_est_span(obj,prop);

            % find tip location
            % 1st screening
            % -------------
            % consider only voxels in a cone of obj.wingana.coneang degrees around the line starting
            % from wingCM along the direction of spanVec
            % 2nd screening
            % -------------
            % take the farthest "fraction" of the voxels and calculate their mean position

            [spanHat,cm_coord4tip,coords_cone] = est_tip_recalc_span(obj,coords,wingCM,spanVec,forCM,prop,name2save);


            % estimate the chord, take the wing strip arounde CM wing and
            % kmeans to 2 groups. define as LE the group that is above the projection of wingCM on Xbody
            [TE_chrd,LE_chrd] = est_chord(obj,forCM,prop);


            % devide the wing to strips, choose the one furtherst from the
            % tip as the new CM
            [cm_coord] = relocateCM(obj,forCM,coords,spanHat',prop,cm_coord4tip);
            [spanHat,cm_coord4tip,new_wCM ] = span_from_tip_CM(obj,coords_cone,forCM,prop,name2save);
            [TE_chrd,LE_chrd] = est_chord(obj,cm_coord,prop);


            if parser.Results.plotcone == 1
                ind2mm = 50e-6*1000;
                coords = (obj.rotmat_EWtoL*(coords'))';
                coords_cone = (obj.rotmat_EWtoL*(coords_cone'))';
                cm_coord = (obj.rotmat_EWtoL*(cm_coord'))';
                bod = (obj.rotmat_EWtoL*double(obj.body.hull3d)')';
                
                plot3(ind2mm*coords(:,1),ind2mm*coords(:,2),ind2mm*coords(:,3),'.r','markersize',5);hold on;
                plot3(ind2mm*coords_cone(:,1),ind2mm*coords_cone(:,2),ind2mm*coords_cone(:,3),'marker','.','color',[0.93,0.69,0.13],'markersize',12)
                plot3(ind2mm*cm_coord(:,1),ind2mm*cm_coord(:,2),ind2mm*cm_coord(:,3),'.k','markersize',6)
                plot3(ind2mm*bod(:,1),ind2mm*bod(:,2),ind2mm*bod(:,3),'.g','markersize',7);axis equal
                
                
            end
            
            if parser.Results.plotChord == 1
                meanslice = mean(forCM);
                plot_3d(coords,'b')


                hold on;quiver3(meanslice(1),meanslice(2),meanslice(3),obj.(prop).vectors.chord(1),obj.(prop).vectors.chord(2),obj.(prop).vectors.chord(3),15)
                hold on;quiver3(meanslice(1),meanslice(2),meanslice(3),obj.(prop).vectors.(name2save)(1),obj.(prop).vectors.(name2save)(2),obj.(prop).vectors.(name2save)(3),30)
                
                
                plot_3d(forCM,'g');hold on
                plot_3d(TE_chrd,'m');hold on;
                plot_3d(LE_chrd,'y')
            end
        end
        
        function [amnt,rw] = getamount_stripe(obj,vxls,spanHat,amnt)
            % perform dot product of voxels (vxls) on a vector (spanHat),
            % return amnt - the amount of voxels in each stripe. and rw -
            % a cell array with the rows of each stripe
            distn = (dot(vxls',repmat(spanHat',1,size(vxls,1))));
            maxpt = max(distn);
            minpt = min(distn);
            sz_st = linspace(minpt,maxpt,amnt);
            
            for k = 1:1:length(sz_st)-1
                rws = find(distn>sz_st(k) & distn<sz_st(k + 1));
                amnt(k) = length(rws);
                rw{k} = rws;
            end
        end
        
        function wingTip = calculateTip(obj,prop,coords,wingCM,varargin)
            % use the recent span vector to find final tip position:
            % the tip is the farthest wing voxel in the direction of the
            % new span (delta tube stuff, no need to explain this here).
            
            parser = inputParser;
            addParameter(parser,'spaname','span'); %
            addParameter(parser,'spanframe',0); %
            addParameter(parser,'tipop',0);
            parse(parser, varargin{:})
            % recalculates new tip as the farthest voxel in the new span
            % direction
            
            DeltMin = obj.wingana.DeltaTubeMin;
            DeltMax = obj.wingana.DeltaTubeMax;
            spnam = parser.Results.spaname;
            if parser.Results.spanframe ~= 0
                spn = (obj.rotmat_EWtoL'*obj.(prop).vectors.(spnam)(parser.Results.spanframe,:)')';
            else

            spn = obj.(prop).vectors.(spnam);
            end
            % mat1 is vector from wing cm (not normalized).
            Nvox = size(coords,1);
            mat1 = double(coords) - repmat(wingCM, Nvox, 1) ;
            mat2 = repmat(spn, Nvox, 1) ;
            

            % cm_coord
            mat3 = mat1 - dot(mat1, mat2, 2) .* mat2 ; % "dot(mat1, mat2, 2)" should be column vector
            
            % calc norm of mat3 rows
            perps = vecnorm(mat3,2,2);
            
            % find the voxels from cm_coord that are inside a tube around
            % the span
            DELTA_TUBE = DeltMin;
            ind1 = perps <= DELTA_TUBE ;
            
            while sum(ind1) == 0
                DELTA_TUBE = DELTA_TUBE+0.05;
                ind1 = perps <= DELTA_TUBE ;
                if DELTA_TUBE >= DeltMax
                    break
                end
            end
            
            %
            mat1_tube = mat1(ind1,:) ;
            dst = dot(mat1_tube',repmat(spn,size(mat1_tube,1),1)');
%             dst = vecnorm(mat1_tube,2,2) ; % distances of "tube" voxels from wing cm.
            [~, maxind] = max(dst) ; % look for farthest voxel in tube
            wingTip_tube = coords(ind1,:);
            wingTip = wingTip_tube(maxind,:) ; % finally update wingTip
            
            if DELTA_TUBE == DeltMax
                wingTip = [999 999 999];
            end
            if parser.Results.tipop == 0
            obj.(prop).coords.tip = wingTip;
            end
            plt = 0;
            if plt ==1 
figure;plot_3d(coords,'r');hold on;plot_3d(wingTip_tube,'b');hold on;plot_3d(wingTip,'k')
meanc = mean(coords);
quiver3(meanc(1),meanc(2),meanc(3),spn(1),spn(2),spn(3),30)
figure;plot_3d(mat1,'r');hold on;plot_3d(mat1_tube,'b');hold on;plot_3d(wingTip,'k')


            end
        end
        
        function hull3d = prep4save_spnXTip(obj,op,body,realCfrm,framevec)
            
            % outputCell = {rw,lw,sp...
            %     ,tip,chordr,spl...
            %     ,tipl,chordl,X,b1,b2,...
            %     b1l,b2l};
            
            
            hull3d.body.hull = body;
            obj.real_coord = realCfrm;
            hull3d.frames =framevec;
            obj.frames = framevec;
            for k = 1:1:length(op)
                if sum(cellfun(@isempty,(op{k})))>0
                    inds = cellfun(@isempty,(op{k}));
                    op{k}(inds) = {[nan nan nan]};
                end
                hull3d.rightwing.hull.hull3d(k) = op{k}(1);
                hull3d.leftwing.hull.hull3d(k) = op{k}(2);
                obj.body.hull3d = [];
                
                obj.rightwing.vectors.span(k,:) = op{k}{3};
                obj.rightwing.coords.tip(k,:) = op{k}{4};
                obj.rightwing.vectors.chord(k,:) = op{k}{5};
                
                obj.leftwing.vectors.span(k,:) = op{k}{6};
                obj.leftwing.coords.tip(k,:) = op{k}{7};
                obj.leftwing.vectors.chord(k,:) = op{k}{8};
                obj.body.vectors.X(k,:) = op{k}{9};
                obj.body.coords.CM_real(k,:) = op{k}{14};
                
                hull3d.rightwing.hull.LE(k) = op{k}(10);
                hull3d.rightwing.hull.TE(k) = op{k}(11);
                
                hull3d.leftwing.hull.LE(k) = op{k}(12);
                hull3d.leftwing.hull.TE(k) = op{k}(13);
                
                hull3d.leftwing.hull.LE(k) = op{k}(12);
                hull3d.leftwing.hull.TE(k) = op{k}(13);
                
                hull3d.rightwing.hull.wingCM(k) = op{k}(15);
                hull3d.leftwing.hull.wingCM(k) = op{k}(16);
 
                
                hull3d.rightwing.hull.tip(k) = op{k}(17);
                hull3d.leftwing.hull.tip(k) = op{k}(18);
            end
            
        end
        
        function labAx_ploter(obj,fr,hull3d,varargin)
            parser = inputParser;
            addParameter(parser,'span',1); % number of camera pointing in Z lab axis
            addParameter(parser,'spanprop','span'); % number of camera pointing in Z lab axis
            addParameter(parser,'strkpln',1); % number of camera pointing in Z lab axis
            addParameter(parser,'cam',1); % number of camera pointing in Z lab axis
            addParameter(parser,'leg',1); % number of camera pointing in Z lab axis
            addParameter(parser,'bodyax',1); % number of camera pointing in Z lab axis
            addParameter(parser,'flipcam',1); % number of camera pointing in Z lab axis
            addParameter(parser,'bound',1); % number of camera pointing in Z lab axis
            addParameter(parser,'LETE',1); % number of camera pointing in Z lab axis
            addParameter(parser,'chord',1); % number of camera pointing in Z lab axis
            addParameter(parser,'normal',1); % number of camera pointing in Z lab axis
            addParameter(parser,'tip',1); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            % plot hull in lab axes
            
            frm = (fr == obj.frames);
            body_plt = hull3d.body.body4plot{frm};
            
            if parser.Results.bound == 1
                rightwing_plt = [hull3d.rightwing.hull.LE{frm};hull3d.rightwing.hull.TE{frm}];
                leftwing_plt = [hull3d.leftwing.hull.LE{frm};hull3d.leftwing.hull.TE{frm}];
            else
                rightwing_plt = hull3d.rightwing.hull.hull3d{frm};
                leftwing_plt = hull3d.leftwing.hull.hull3d{frm};
                
            end
            
            % transform indices to lab axis
            
            spanprop = parser.Results.spanprop;
            realC = obj.real_coord{frm}';
            [hullRealbody rotmat] = Functions.hullRec2lab(body_plt,obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            
            hullrightwing = Functions.hullRec2lab(rightwing_plt,obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            hullleftwing = Functions.hullRec2lab(leftwing_plt,obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            if parser.Results.LETE == 1
                hullrightwing_LE = Functions.hullRec2lab(hull3d.rightwing.hull.LE{frm},obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
                hullrightwing_TE = Functions.hullRec2lab(hull3d.rightwing.hull.TE{frm},obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
                hullleftwing_LE = Functions.hullRec2lab(hull3d.leftwing.hull.LE{frm},obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
                hullleftwing_TE = Functions.hullRec2lab(hull3d.leftwing.hull.TE{frm},obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            end
            [tipR rotmat] = Functions.hullRec2lab(obj.rightwing.coords.tip(frm,:),obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            [tipL rotmat] = Functions.hullRec2lab(obj.leftwing.coords.tip(frm,:),obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            strk = obj.body.vectors.strkPlan(frm,:);
            bodyCM = mean(hullRealbody);
            rwingCM = mean(hullrightwing);
            lwingCM = mean(hullleftwing);
            % plot wing and body
            lg(1) = plot3(hullRealbody(:,1),hullRealbody(:,2),hullRealbody(:,3),'.g');hold on
            if parser.Results.LETE == 0
                lg(2) = plot3(hullrightwing(:,1),hullrightwing(:,2),hullrightwing(:,3),'.r');
                lg(3) = plot3(hullleftwing(:,1),hullleftwing(:,2),hullleftwing(:,3),'.b');hold on
            end
            if parser.Results.cam == 1
                for k= 1:1:length(obj.cameras.cams_arr)
                    dircam = obj.cameras.all.Rotation_Matrix*(parser.Results.flipcam * obj.cameras.cams_arr(k).camera_dir);
                    lg(4) = quiver3(bodyCM(1),bodyCM(2),bodyCM(3),dircam(1),dircam(2),dircam(3),10,'color',[0.5 0.5 0.5],'linewidth',2);hold on
                end
            end
            rspn = obj.rightwing.vectors.(spanprop)(frm,:)';
            lspan = obj.leftwing.vectors.(spanprop)(frm,:)';
            rchord = obj.rightwing.vectors.chord(frm,:)';
            lchord = obj.leftwing.vectors.chord(frm,:)';
            strkpln = obj.body.vectors.strkPlan(frm,:)';
            %             rnrml = obj.rightwing.vectors.nrml(frm,:)';
            %             lnrml = obj.leftwing.vectors.nrml(frm,:)';
            if parser.Results.bodyax == 1
                hold on;lg(5) = quiver3(bodyCM(1),bodyCM(2),bodyCM(3),obj.body.vectors.X(frm,1),obj.body.vectors.X(frm,2),obj.body.vectors.X(frm,3),0.005,'color','k','linewidth',2);
            end
            if parser.Results.span == 1
                hold on;lg(6) = quiver3(rwingCM(1),rwingCM(2),rwingCM(3),rspn(1),rspn(2),rspn(3),0.002,'color','r','linewidth',2);
                hold on;lg(7) = quiver3(lwingCM(1),lwingCM(2),lwingCM(3),lspan(1),lspan(2),lspan(3),0.002,'color','b','linewidth',2);
            end
            if parser.Results.strkpln == 1
                quiver3(bodyCM(1),bodyCM(2),bodyCM(3),strkpln(1),strkpln(2),strkpln(3),0.005,'color','k','linewidth',2);
            end
            
            if parser.Results.tip == 1
                lg(8) = plot3(tipR(:,1),tipR(:,2),tipR(:,3),'sk','MarkerFaceColor','k','markersize',10);hold on
                plot3(tipL(:,1),tipL(:,2),tipL(:,3),'sk','MarkerFaceColor','k','markersize',10);hold on
            end
            
            if parser.Results.chord == 1
                hold on;lg(9) = quiver3(rwingCM(1),rwingCM(2),rwingCM(3),rchord(1),rchord(2),rchord(3),0.001,'color','r','linewidth',2);
                hold on;quiver3(lwingCM(1),lwingCM(2),lwingCM(3),lchord(1),lchord(2),lchord(3),0.001,'color','b','linewidth',2);
            end
            if parser.Results.normal == 1
                hold on;lg(10) = quiver3(rwingCM(1),rwingCM(2),rwingCM(3),rnrml(1),rnrml(2),rnrml(3),0.001,'color','r','linewidth',2);
                hold on;quiver3(lwingCM(1),lwingCM(2),lwingCM(3),lnrml(1),lnrml(2),lnrml(3),0.001,'color','b','linewidth',2);
            end
            if parser.Results.LETE == 1
                lg(2) = plot3(hullrightwing_LE(:,1),hullrightwing_LE(:,2),hullrightwing_LE(:,3),'.r');hold on
                plot3(hullrightwing_TE(:,1),hullrightwing_TE(:,2),hullrightwing_TE(:,3),'.m');
                lg(3) = plot3(hullleftwing_LE(:,1),hullleftwing_LE(:,2),hullleftwing_LE(:,3),'.b');hold on
                plot3(hullleftwing_TE(:,1),hullleftwing_TE(:,2),hullleftwing_TE(:,3),'.c');
            end
            
            
            
            xlabel('mm');ylabel('mm');zlabel('mm');
            ttl = sprintf('frame %d',fr);
            title(ttl);
            axis equal
            grid on
            box on
            if parser.Results.leg == 1
                legend(lg,'body','right wing','left wing','cameras','Xbody','span','span','tip','chord');
            end
        end
        
        function hull3d = calcNormChord_bound(obj,hull3d,frm)
            % calculate the wing's chord and normal from the boundary
            wingnameC = {'rightwing','leftwing'};
            for k = 1:1:2
                try
                    % rotate the boundary to lab axis and fit a plane to
                    % calculate the normal
                    boundWing = (obj.rotmat_EWtoL*([double(hull3d.(wingnameC{k}).hull.LE{frm});double(hull3d.(wingnameC{k}).hull.TE{frm})]'))';
                    nrml = Functions.affine_fit(boundWing);
                    obj.(wingnameC{k}).vectors.nrml(frm,:) = nrml / norm(nrml);
                    chord = cross(obj.(wingnameC{k}).vectors.span(frm,:), obj.(wingnameC{k}).vectors.nrml(frm,:));
                    dirchord_bound = dot(chord,obj.body.vectors.strkPlan(frm,:));
                    dirchord = dot(obj.(wingnameC{k}).vectors.chord(frm,:),obj.body.vectors.strkPlan(frm,:));
                    obj.(wingnameC{k}).vectors.chord_bound(frm,:) = chord/norm(chord);






                    if dirchord_bound < 0
                        obj.(wingnameC{k}).vectors.nrml(frm,:) = -obj.(wingnameC{k}).vectors.nrml(frm,:);
                        obj.(wingnameC{k}).vectors.chord_bound(frm,:) = -obj.(wingnameC{k}).vectors.chord_bound(frm,:);

%                         TEemp = hull3d.rightwing.hull.LE{frm};
%                         hull3d.rightwing.hull.LE{frm} = hull3d.rightwing.hull.TE{frm};
%                         hull3d.rightwing.hull.TE{frm} = TEemp;
                    end
                    if dirchord < 0
                        obj.(wingnameC{k}).vectors.chord(frm,:) = -obj.(wingnameC{k}).vectors.chord(frm,:);
                    end

                    if sum(isnan(obj.(wingnameC{k}).vectors.chord_bound(frm,:)))==3
                    obj.(wingnameC{k}).vectors.chord_bound(frm,:) = obj.(wingnameC{k}).vectors.chord(frm,:);
                    obj.(wingnameC{k}).vectors.nrml(frm,:) = cross(obj.(wingnameC{k}).vectors.span(frm,:),  obj.(wingnameC{k}).vectors.chord(frm,:));

                    end

                    LEemp = (obj.rotmat_EWtoL*double(hull3d.(wingnameC{k}).hull.LE{frm})')';
                    meanLE = mean(LEemp,1);
                    TEemp = (obj.rotmat_EWtoL*double(hull3d.(wingnameC{k}).hull.TE{frm})')';
                    meanTE = mean(TEemp,1);

                    dotLEChord = dot(meanLE,obj.(wingnameC{k}).vectors.chord_bound(frm,:));
                    dotTEChord = dot(meanTE,obj.(wingnameC{k}).vectors.chord_bound(frm,:));
                    if dotLEChord < dotTEChord



                        LEemp2 = (double(hull3d.(wingnameC{k}).hull.LE{frm})')';
                        hull3d.(wingnameC{k}).hull.LE{frm} = hull3d.(wingnameC{k}).hull.TE{frm};
                        hull3d.(wingnameC{k}).hull.TE{frm} = LEemp2;
                    end




                catch
                    obj.(wingnameC{k}).vectors.chord_bound(frm,:) = obj.(wingnameC{k}).vectors.chord(frm,:);
                    obj.(wingnameC{k}).vectors.nrml(frm,:) = cross(obj.(wingnameC{k}).vectors.span(frm,:),  obj.(wingnameC{k}).vectors.chord(frm,:));
                end
            end
        end
%          if dirchord_bound < 0
%                         obj.(wingnameC{k}).vectors.nrml(frm,:) = -obj.(wingnameC{k}).vectors.nrml(frm,:);
%                         obj.(wingnameC{k}).vectors.chord_bound(frm,:) = -obj.(wingnameC{k}).vectors.chord_bound(frm,:);
%                         
%                         TEemp = (hull3d.rightwing.hull.LE{frm}')';
%                         hull3d.rightwing.hull.LE{frm} = hull3d.rightwing.hull.TE{frm};
%                         hull3d.rightwing.hull.TE{frm} = TEemp;
% 
%                     end
        function idx4StrkPln = ChooseSpan(obj,varargin)
            parser = inputParser;
            addParameter(parser,'plot',0); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            % choose the frames used to calculate the Y axis. In those
            % frames the wings are the farthest apart
            
            % projection of each wing span on body axis
            dotspanAx_wing1=dot(obj.rightwing.vectors.span', obj.body.vectors.X');
            dotspanAx_wing2=dot(obj.leftwing.vectors.span', obj.body.vectors.X');
            dotspanAx_wing2 = filloutliers(dotspanAx_wing2,'pchip');
            dotspanAx_wing1 = filloutliers(dotspanAx_wing1,'pchip');
            
            % Calculate and choose the greatest angles (between both wings)
            distSpans=real(acosd(dot(obj.rightwing.vectors.span',obj.leftwing.vectors.span')));
            angBodSp=real(acosd(dot([obj.rightwing.vectors.span'],obj.body.vectors.X')));
            
            % calculate the mean projection of both spans on body axis. if
            % the sign changes it means the wings are located aproximatly 90 degrees to the body axis
            mean_strks=mean(([dotspanAx_wing1;dotspanAx_wing2]));
            changeSgn=[mean_strks<0;mean_strks>=0];
            [FrbckStrk] = FindUp_downStrk(obj,changeSgn,mean_strks,0);
            [idx4StrkPln] = Choose_GrAng_wing1_wing2(obj,distSpans,FrbckStrk,obj.wingana.angleTH,10);
            % make sure there is only one wing flap between each cosen
            % points
            while sum(diff(idx4StrkPln)<65)>0
                for k= 1:1:length(idx4StrkPln)-1
                    if idx4StrkPln(k+1)-idx4StrkPln(k)<65
                        idx4StrkPln(k+1) = [];
                        break
                    end
                end
            end
            % leave only indices that are aproximatly 90 deg from body axis
            idx4StrkPln=unique(idx4StrkPln);
            idx4StrkPln(idx4StrkPln<1)=[];
            idx4StrkPln(abs(angBodSp(idx4StrkPln)-90)>20) = [];
            if parser.Results.plot == 1
                figure;
                plot(mean_strks);hold on;grid on;title('Projection of mean span on body axis');xlabel('frames');ylabel('angle [deg]')
                plot(idx4StrkPln,mean_strks(idx4StrkPln),'*r')
                plot(idx4StrkPln,mean_strks(idx4StrkPln),'*k')
                
                figure;
                plot(distSpans);hold on;grid on;title('angle between spans');xlabel('frames');ylabel('angle [deg]')
                plot(idx4StrkPln,distSpans(idx4StrkPln),'*r')
                plot(idx4StrkPln,distSpans(idx4StrkPln),'*k')
            end
            
        end
        
        function obj = BodyAxesYZ(obj)
            
            idx4StrkPln = ChooseSpan(obj);
            spn_wing1 = obj.rightwing.vectors.span;
            spn_wing2 = obj.leftwing.vectors.span;
            
            % make the span direction of both wings the same - we dont know
            % which is right/left
            bodAxCrosSpan = [cross([obj.body.vectors.X;obj.body.vectors.X],...
                [spn_wing1;spn_wing2])];
            FlipDirWing = dot([repmat([0,0,1],size(bodAxCrosSpan,1),1)'],bodAxCrosSpan');
            [ind] = ([FlipDirWing(1:size(FlipDirWing,2)/2)',FlipDirWing(size(FlipDirWing,2)/2+1:end)'])<0;
            spn_wing1(ind(:,1),:) = -spn_wing1(ind(:,1),:);
            spn_wing2(ind(:,2),:) = -spn_wing2(ind(:,2),:);
            
            for k= 1:1:length(idx4StrkPln)
                
                % define a window of C around idx4StrkPln(k) (the point in
                % ehich the wings are at maximum angle from each other.
                inifr = idx4StrkPln(k)-10;
                if inifr < 1
                    inifr = 1;
                end
                maxfr = idx4StrkPln(k)+10;
                if maxfr > size(obj.rightwing.coords.tip,1) || maxfr > size(obj.leftwing.coords.tip,1)
                    szLR = [size(obj.rightwing.coords.tip,1),size(obj.leftwing.coords.tip,1)];
                    maxfr =max(szLR)  ;
                end
                %-----------------
                
                NormSpan_rw{k} = obj.rightwing.vectors.span(inifr:maxfr,:)./vecnorm(obj.rightwing.vectors.span(inifr:maxfr,:)',1)';
                Zvec = NormSpan_rw{k}(:,3).*repmat([0,0,1],size(NormSpan_rw{k},1),1);
                
                % find 2 angles that represent the location of the span,
                % use it to find outliers.
                vec4ang = NormSpan_rw{k} - Zvec;
                angXY{k} = atan2(vec4ang(:,2),vec4ang(:,1));
                angZ{k} = acos(dot(NormSpan_rw{k}',repmat([0,0,1],size(NormSpan_rw{k},1),1)')');
                NormSpan_lw{k} = obj.leftwing.vectors.span(inifr:maxfr,:)./vecnorm(obj.leftwing.vectors.span(inifr:maxfr,:)',1)';
                Zvec = NormSpan_lw{k}(:,3).*repmat([0,0,1],size(NormSpan_lw{k},1),1);
                
                vec4ang = NormSpan_lw{k} - Zvec;
                angXY_lw{k} = atan2(vec4ang(:,2),vec4ang(:,1));
                angZ_lw{k} = acos(dot(NormSpan_lw{k}',repmat([0,0,1],size(NormSpan_lw{k},1),1)')');
                frms{k} = [inifr:maxfr];
                
                angles_spn_rw = [(angXY{k}'),(angZ{k}')];
                angles_spn_lw = [(angXY_lw{k}'),(angZ_lw{k}')];
                
                [ro_r,col] = find(isoutlier(angles_spn_rw));
                [ro_l,col] = find(isoutlier(angles_spn_lw));
                
                [Crw,ia ]  = setdiff([inifr:maxfr] , ro_r);
                [Clw,ia ]  = setdiff([inifr:maxfr] , ro_l);
                
                
                %                 calculate Ybody as the mean of both spans (when the angle
                %                 between them is maximal).
                sptmp =   [mean(spn_wing1(Crw,:),1)+mean(spn_wing2(Clw,:),1)]/2;
                sptmp = sptmp/norm(sptmp);
                Ybody(k,1:3) =sptmp;
                
                % find plane of wing tips (n) - its aproximatly the stroke plane. Check the the
                % direction of the calculated Y axis is the same as [stroke plane X Xbody]
                TpWng = [double([obj.rightwing.coords.tip]) - obj.body.coords.CM;double([obj.leftwing.coords.tip]) - obj.body.coords.CM];
                TpWng(find((sum(isnan(TpWng),2) == 3)),:) = [];
                n = Functions.affine_fit(TpWng);
                CheckDirY = cross(n,obj.body.vectors.X(k,:));
                CheckDirY = CheckDirY/norm(CheckDirY);
                
                if dot(CheckDirY,Ybody(k,1:3))<0
                    Ybody(k,1:3) = -  Ybody(k,1:3);
                end
                
                if k>1 && dot(Ybody(k,1:3),Ybody(k-1,1:3))<0
                    Ybody(k,1:3) = -  Ybody(k,1:3);
                end
            end
            
            % interpulate the points in between the maximal angle to get
            % a continius Y axis
            Ybody_inter = [interp1(idx4StrkPln,Ybody(:,1),1:1:size(obj.body.vectors.X,1),'spline')', interp1(idx4StrkPln,Ybody(:,2),1:1:size(obj.body.vectors.X,1),'spline')'...
                ,interp1(idx4StrkPln,Ybody(:,3),1:1:size(obj.body.vectors.X,1),'spline')'];
            Ybody_inter = Ybody_inter - obj.body.vectors.X .* dot(obj.body.vectors.X',Ybody_inter')'; % make Ybody perpendicular to Xbody (for interpulated points)
            obj.body.vectors.Y = Ybody_inter./vecnorm(Ybody_inter,2,2);
            Zbody = cross(obj.body.vectors.X,obj.body.vectors.Y); % calculate Z body (Xbody X Ybody)
            obj.body.vectors.Z = Zbody./vecnorm(Zbody,2,2);
            obj.body.idxInterp = idx4StrkPln;
            
            % calculate stroke plane as a 45 deg rotation of Xbody around Y body
            for fr = 1:1:size(obj.body.vectors.X,1)
                strkPlan = Functions.rodrigues_rot(obj.body.vectors.X(fr,:),obj.body.vectors.Y(fr,:),-45*pi/180);
                if dot(strkPlan,[0,0,1])<0
                    strkPlan = Functions.rodrigues_rot(obj.body.vectors.X(fr,:),obj.body.vectors.Y(fr,:),45*pi/180);
                end
                obj.body.vectors.strkPlan(fr,1:3) = strkPlan;
            end
            
        end
        
        function obj = BodyAxesYZV2(obj)
            
            idx4StrkPln = ChooseSpan(obj);
            spn_wing1 = obj.rightwing.vectors.span;
            spn_wing2 = obj.leftwing.vectors.span;
            
            % make the span direction of both wings the same - we dont know
            % which is right/left
            bodAxCrosSpan = [cross([obj.body.vectors.X;obj.body.vectors.X],...
                [spn_wing1;spn_wing2])];
            FlipDirWing = dot([repmat([0,0,1],size(bodAxCrosSpan,1),1)'],bodAxCrosSpan');
            [ind] = ([FlipDirWing(1:size(FlipDirWing,2)/2)',FlipDirWing(size(FlipDirWing,2)/2+1:end)'])<0;
            spn_wing1(ind(:,1),:) = -spn_wing1(ind(:,1),:);
            spn_wing1(ind(:,2),:) = -spn_wing2(ind(:,2),:);
            
            for k= 1:1:length(idx4StrkPln)
                
                % define a window of C around idx4StrkPln(k) (the point in
                % ehich the wings are at maximum angle from each other.
                inifr = idx4StrkPln(k)-10;
                if inifr < 1
                    inifr = 1;
                end
                maxfr = idx4StrkPln(k)+10;
                if maxfr > size(obj.rightwing.coords.tip,1) || maxfr > size(obj.leftwing.coords.tip,1)
                    szLR = [size(obj.rightwing.coords.tip,1),size(obj.leftwing.coords.tip,1)];
                    maxfr =max(szLR)  ;
                end
                %-----------------
                
                NormSpan_rw{k} = spn_wing1(inifr:maxfr,:)./vecnorm(spn_wing1(inifr:maxfr,:)',1)';
                Zvec = NormSpan_rw{k}(:,3).*repmat([0,0,1],size(NormSpan_rw{k},1),1);
                
                % find 2 angles that represent the location of the span,
                % use it to find outliers.
                vec4ang = NormSpan_rw{k} - Zvec;
                angXY{k} = atan2(vec4ang(:,2),vec4ang(:,1));
                angZ{k} = acos(dot(NormSpan_rw{k}',repmat([0,0,1],size(NormSpan_rw{k},1),1)')');
                NormSpan_lw{k} = spn_wing2(inifr:maxfr,:)./vecnorm(spn_wing2(inifr:maxfr,:)',1)';
                Zvec = NormSpan_lw{k}(:,3).*repmat([0,0,1],size(NormSpan_lw{k},1),1);
                
                vec4ang = NormSpan_lw{k} - Zvec;
                angXY_lw{k} = atan2(vec4ang(:,2),vec4ang(:,1));
                angZ_lw{k} = acos(dot(NormSpan_lw{k}',repmat([0,0,1],size(NormSpan_lw{k},1),1)')');
                frms_cell{k} = [inifr:maxfr];
                
                angles_spn_rw = [(angXY{k}'),(angZ{k}')];
                angles_spn_lw = [(angXY_lw{k}'),(angZ_lw{k}')];
                
                [ro_r,col] = find(isoutlier(angles_spn_rw));
                [ro_l,col] = find(isoutlier(angles_spn_lw));
                frms = [inifr:maxfr];
                Clw = frms;
                Clw(ro_l) = [];
                Crw = frms;
                Crw(ro_r) = [];
                
                tips = [obj.rightwing.coords.tip(Crw,:);obj.leftwing.coords.tip(Clw,:)];
                tips_lab = (obj.rotmat_EWtoL*tips')';
                row2del = find([sum(isnan(tips_lab),2) > 0]);
                tips_lab(row2del,:) = [];
                vec_norm_tips = Functions.affine_fit(tips_lab)';
                
                
                Ybody(k,1:3) =cross(vec_norm_tips,obj.body.vectors.X(idx4StrkPln(k),:));
                Ybody(k,1:3) = Ybody(k,1:3)/norm(Ybody(k,1:3));
                
                
                %                 [Crw,ia ]  = setdiff([inifr:maxfr] , ro_r);
                %                 [Clw,ia ]  = setdiff([inifr:maxfr] , ro_l);
                %
                %
                % %                 calculate Ybody as the mean of both spans (when the angle
                % %                 between them is maximal).
                %                 sptmp =   [mean(spn_wing1(Crw,:),1)+mean(spn_wing2(Clw,:),1)]/2;
                %                 sptmp = sptmp/norm(sptmp);
                %                 Ybody(k,1:3) =sptmp;
                %
                % find plane of wing tips (n) - its aproximatly the stroke plane. Check the the
                % direction of the calculated Y axis is the same as [stroke plane X Xbody]
                TpWng = [double([obj.rightwing.coords.tip]) - obj.body.coords.CM;double([obj.leftwing.coords.tip]) - obj.body.coords.CM];
                TpWng(find((sum(isnan(TpWng),2) == 3)),:) = [];
                
                TpWng = (obj.rotmat_EWtoL*TpWng')';
                n = Functions.affine_fit(TpWng);
                CheckDirY = cross(n,obj.body.vectors.X(k,:));
                CheckDirY = CheckDirY/norm(CheckDirY);
                
                if dot(CheckDirY,Ybody(k,1:3))<0
                    Ybody(k,1:3) = -  Ybody(k,1:3);
                end
                
                if k>1 && dot(Ybody(k,1:3),Ybody(k-1,1:3))<0
                    Ybody(k,1:3) = -  Ybody(k,1:3);
                end
            end
            
            % interpulate the points in between the maximal angle to get
            % a continous Y axis
            Ybody_inter = [interp1(idx4StrkPln,Ybody(:,1),1:1:size(obj.body.vectors.X,1),'spline')', interp1(idx4StrkPln,Ybody(:,2),1:1:size(obj.body.vectors.X,1),'spline')'...
                ,interp1(idx4StrkPln,Ybody(:,3),1:1:size(obj.body.vectors.X,1),'spline')'];
            Ybody_inter = Ybody_inter - obj.body.vectors.X .* dot(obj.body.vectors.X',Ybody_inter')'; % make Ybody perpendicular to Xbody (for interpulated points)
            obj.body.vectors.Y = Ybody_inter./vecnorm(Ybody_inter,2,2);
            Zbody = cross(obj.body.vectors.X,obj.body.vectors.Y); % calculate Z body (Xbody X Ybody)
            obj.body.vectors.Z = Zbody./vecnorm(Zbody,2,2);
            obj.body.idxInterp = idx4StrkPln;
            
            % calculate stroke plane as a 45 deg rotation of Xbody around Y body
            for fr = 1:1:size(obj.body.vectors.X,1)
                strkPlan = Functions.rodrigues_rot(obj.body.vectors.X(fr,:),obj.body.vectors.Y(fr,:),-45*pi/180);
                if dot(strkPlan,[0,0,1])<0
                    strkPlan = Functions.rodrigues_rot(obj.body.vectors.X(fr,:),obj.body.vectors.Y(fr,:),45*pi/180);
                end
                obj.body.vectors.strkPlan(fr,1:3) = strkPlan;
            end
            
        end
        
        
        function hull3d = RightLeftWing(obj,hull3d)
            % define right and left wing according to Y body
            Right = find(dot(obj.leftwing.vectors.span',obj.body.vectors.Y')'<dot(obj.rightwing.vectors.span',obj.body.vectors.Y')');
            if isempty(Right) == 0
                tmpRspn = obj.rightwing.vectors.span;
                tmpRspn(Right,:) = obj.leftwing.vectors.span(Right,:);
                obj.leftwing.vectors.span(Right,:) = obj.rightwing.vectors.span(Right,:);
                obj.rightwing.vectors.span = tmpRspn;
                
                tmpRcrd = obj.rightwing.vectors.chord;
                tmpRcrd(Right,:) = obj.leftwing.vectors.chord(Right,:);
                obj.leftwing.vectors.chord(Right,:) = obj.rightwing.vectors.chord(Right,:);
                obj.rightwing.vectors.chord = tmpRcrd;
                
                tmpRtip = obj.rightwing.coords.tip;
                tmpRtip(Right,:) = obj.leftwing.coords.tip(Right,:);
                obj.leftwing.coords.tip(Right,:) = obj.rightwing.coords.tip(Right,:);
                obj.rightwing.coords.tip = tmpRtip;
                
                tmpRhull = hull3d.rightwing.hull.hull3d;
                tmpRhull(Right) = hull3d.leftwing.hull.hull3d(Right);
                hull3d.leftwing.hull.hull3d(Right) = hull3d.rightwing.hull.hull3d(Right);
                hull3d.rightwing.hull.hull3d = tmpRhull;
                
                tmpRhull = hull3d.rightwing.hull.LE;
                tmpRhull(Right) = hull3d.leftwing.hull.LE(Right);
                hull3d.leftwing.hull.LE(Right) = hull3d.rightwing.hull.LE(Right);
                hull3d.rightwing.hull.LE = tmpRhull;
                
                tmpRhull = hull3d.rightwing.hull.TE;
                tmpRhull(Right) = hull3d.leftwing.hull.TE(Right);
                hull3d.leftwing.hull.TE(Right) = hull3d.rightwing.hull.TE(Right);
                hull3d.rightwing.hull.TE = tmpRhull;
                
                tmpRspn = hull3d.rightwing.hull.wingCM;
                tmpRspn(Right) = hull3d.leftwing.hull.wingCM(Right);
                hull3d.leftwing.hull.wingCM(Right) = hull3d.rightwing.hull.wingCM(Right);
                hull3d.rightwing.hull.wingCM = tmpRspn;

                tmpRspn = hull3d.rightwing.hull.tip;
                tmpRspn(Right) = hull3d.leftwing.hull.tip(Right);
                hull3d.leftwing.hull.tip(Right) = hull3d.rightwing.hull.tip(Right);
                hull3d.rightwing.hull.tip = tmpRspn;
                
            end
        end
        
        function hull3d = flipLETE(obj,hull3d)
            wingnamecell = {'rightwing','leftwing'};
            for kwing = 1:1:2
                %                 chord2flip = dot(obj.body.vectors.strkPlan',obj.(wingnamecell{kwing}).vectors.chord');
                chord2flip = cross(obj.(wingnamecell{kwing}).vectors.span',obj.(wingnamecell{kwing}).vectors.chord');
                chord2flip = dot(obj.body.vectors.strkPlan',chord2flip);
                chord2flip = chord2flip<0;
                obj.(wingnamecell{kwing}).vectors.chord(chord2flip,:) = -obj.(wingnamecell{kwing}).vectors.chord(chord2flip,:);
                
                LEtmp = hull3d.(wingnamecell{kwing}).hull.LE ;
                hull3d.(wingnamecell{kwing}).hull.LE(chord2flip) = hull3d.(wingnamecell{kwing}).hull.TE(chord2flip);
                hull3d.(wingnamecell{kwing}).hull.TE(chord2flip) = LEtmp(chord2flip);
            end
        end
        
        function [LE TE] = SplitWin_LETE(obj,chord_dir,Tip,wing,varargin)
            parser = inputParser;
            addParameter(parser,'plot',0); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            % split wing to LE and TE according to chord and tip
            
            % project tip on chord, split to LE and TE if above and below
            % tip projection on chord
            z_tipCoords_wing=dot(chord_dir,double(Tip));
            z_wingCoords_wing=dot(repmat(chord_dir',1,size(wing,1)),wing');
            LE_indx_wing=z_wingCoords_wing>z_tipCoords_wing;
            LE=double(wing(LE_indx_wing==1,:));
            TE=double(wing(LE_indx_wing==0,:));
            
            if parser.Results.plot == 1
                pcshow(wing,'r')
                meanRwin = mean(wing);
                hold on; quiver3(meanRwin(1),meanRwin(2),meanRwin(3),chord_dir(1),chord_dir(2),chord_dir(3),10);
                pcshow(LE,'m');hold on
                pcshow(TE,'b');
            end
            
        end
        
        function [LE_srt,TE_srt,ind,LE_srt2d,TE_srt2d] = CreateBoundary(obj,realC,LE,TE,wingname)
            kcamind = 1;
            for kcam = obj.cameras.all.camvec
                TwoD = Functions.Create2DCoords(realC,obj.cameras.all.DLT_coefs(:,kcam),{LE...
                    ,TE,obj.(wingname).hull3d,obj.(wingname).coords.tip,mean(obj.(wingname).hull3d),obj.body.hull3d},1,obj.cameras.all.size_image(1),obj.cameras.all.RotMat_vol);
                [TwoDwingIm,TwoDDownIm,TwoDUpIm] = Split2Dbound(obj,TwoD,obj.wingana.image2D_close);
                % project 3d LE and TE and define the intersecting
                % pixels as LE and TE accordingaly
                BodOnW=0;LETE_same=0;
                
                if size(intersect(TwoD{3},TwoD{6},'rows'),1)>size(unique(TwoD{3},'rows'),1)*3/4
                    % if more than 75% of the wing intersects the body use
                    % the LE hull (not boundary)
                    wingLE = {LE,LE};
                    wingTE = {TE,TE};BodOnW=1;
                else
                    [wingLE,wingTE,wingLE2d,wingTE2d] = intersect2D(obj,TwoDwingIm,LE,TE,TwoD{1},TwoD{2},TwoDUpIm,TwoDDownIm);
                end
                % Each wingLE/wingTE has 2 options, each option is constructed
                % using one of the 2D boundary. (each boundary can
                % intersect with the LE as well as the TE. usually the
                % right option will contain more pixels)
                idxs_LETE = [cell2mat(cellfun(@(x) length(x), wingLE, 'UniformOutput', false));...
                    cell2mat(cellfun(@(x) length(x), wingTE, 'UniformOutput', false))];
                [~,indM]=max([idxs_LETE(1,1)+idxs_LETE(2,2),idxs_LETE(1,2)+idxs_LETE(2,1)]);
                if indM==1
                    op_indLE = 1;op_indTE=2;
                else
                    op_indLE = 2;op_indTE=1;
                end
                % calculate the difference between amount of pixels for
                % each option. keep only the indices marked as LE or TE
                % in the 3D hull. (eventually having 3 hulls, one for
                % each camera)
                
                valDiff(kcamind)=abs(idxs_LETE(1,op_indLE)+idxs_LETE(2,op_indTE)-(idxs_LETE(1,op_indTE)+idxs_LETE(2,op_indLE)));
                LE_srt{kcamind}= [wingLE(op_indLE) wingLE(op_indTE)];
                TE_srt{kcamind}= [wingTE(op_indTE) wingTE(op_indLE)];
                LE_srt2d{kcamind}= [wingLE2d(op_indLE) wingLE2d(op_indTE)];
                TE_srt2d{kcamind}= [wingTE2d(op_indTE) wingTE2d(op_indLE)];
                [~,ind]=sort(valDiff(:),'descend');
                kcamind = kcamind + 1;
            end
        end
        
        function boundOptions(obj,LE_srt,TE_srt,wingname)
            % define LE and TE by the ammount of voxels in each edge.
            % choose the combination with the maximal voxels
            for opt = 1:1:2
                LE_optio = [(LE_srt{1}{opt});(LE_srt{2}{opt});(LE_srt{3}{opt});(LE_srt{4}{opt})];
                occ = Functions.countRepRow(LE_optio);
                LE_op{opt} = occ(occ(:,4) > 2,1:3);
                TE_optio = [(TE_srt{1}{opt});(TE_srt{2}{opt});(TE_srt{3}{opt});(TE_srt{4}{opt})];
                occ = Functions.countRepRow(TE_optio);
                TE_op{opt} = occ(occ(:,4) > 2,1:3);
            end
            
            lentot = cellfun(@length, LE_op) + cellfun(@length, TE_op);
            [~,ind] = max(lentot);
            obj.(wingname).hull.LE = [LE_op{ind}];
            obj.(wingname).hull.TE = [TE_op{ind}];
            
        end
        
        function [wingLE,wingTE,wingLE2d,wingTE2d,bounWing_dil] = intersect2D(obj,TwoDwingIm,LE_wing,TE_wing,TwoDwingLE,TwoDwingTE,Bup,Bdown)
            % Create the boundary of the wing from the image and dilate it
            LETE_same=0;
            [wing_coordsy wing_coordsx]=find(TwoDwingIm);
            CropedIm=TwoDwingIm(min(wing_coordsy):max(wing_coordsy),min(wing_coordsx):max(wing_coordsx));
            CropedIm=imclose(CropedIm,obj.wingana.image2D_close);
            bounWing_dil = bwperim(CropedIm);
            
            bounWing_dil = imdilate(bounWing_dil,obj.wingana.bound_dilate);
            
            
            bounWing{1}=bounWing_dil.*Bdown(min(wing_coordsy):max(wing_coordsy),min(wing_coordsx):max(wing_coordsx));
            bounWing{2}=bounWing_dil.*Bup(min(wing_coordsy):max(wing_coordsy),min(wing_coordsx):max(wing_coordsx));
            
            
            for k=1:1:2
                % count the amount of intersecting pixels for each boundary
                % with the projection of LE and TE
                [boundRo,boundCo]=find(bounWing{k});
                Boundary_wing{k}=[boundCo+min(wing_coordsx)-1,boundRo+min(wing_coordsy)-1];
                
                ibLE=(ismember(TwoDwingLE,Boundary_wing{k},'rows'));
                ibTE=(ismember(TwoDwingTE,Boundary_wing{k},'rows'));
                wingLE{k}=(LE_wing(ibLE,:)')';
                wingTE{k}=(TE_wing(ibTE,:)')';
                
                wingLE2d{k}=(TwoDwingLE(ibLE,:)')';
                wingTE2d{k}=(TwoDwingTE(ibTE,:)')';
                
                
                if abs(sum(ismember(TwoDwingLE,TwoDwingTE,'rows')))>size(TwoDwingLE,1)*4/5
                    LETE_same = 1;
                end
            end
        end
        
        function bodyAngles(obj)
            % calculate pitch, roll and yaw. the body angles.
            % pitch ---------
            obj.body.angles.pitch =(90 - acos(dot(obj.body.vectors.X',repmat([0,0,1],size(obj.body.vectors.X,1) ,1)'))*180/pi)';
            
            % yaw ---------
            obj.body.angles.yaw = (atan2(obj.body.vectors.X(:,2), obj.body.vectors.X(:,1)) * 180/pi)';
            flipyaw = obj.body.angles.yaw<0;
            if isempty(flipyaw) == 0
                obj.body.angles.yaw(flipyaw) = obj.body.angles.yaw(flipyaw) + 360;
            end
            
            % roll--------
            roll_rotation_mat =  Functions.eulerRotationMatrix(obj.body.angles.yaw*pi/180, obj.body.angles.pitch*pi/180,0);
            Y = obj.body.vectors.Y;
            Yvec = permute(Y,[3,2,1]);
            Yvect = pagetranspose(Yvec);
            yflyZerPitch = pagemtimes(roll_rotation_mat,Yvect);
            obj.body.angles.roll = atan2(squeeze(yflyZerPitch(3,:,:)),squeeze(yflyZerPitch(2,:,:)))*180/pi;
        end
        
        function wingAngles(obj,varargin)
            parser = inputParser;
            addParameter(parser,'spanProp','span'); % number of camera pointing in Z lab axis
            addParameter(parser,'chordProp','chord'); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            
            spanProp = parser.Results.spanProp;
            chordProp = parser.Results.chordProp;
            wingname = {'rightwing','leftwing'};
            strkpln = obj.body.vectors.strkPlan;
            for k = 1:1:2
                spn = obj.(wingname{k}).vectors.(spanProp);
                
                [spnOnPln]=Functions.projection(spn,strkpln);
                [XbodOnPln]=Functions.projection(obj.body.vectors.X,strkpln);
                
                % calculate phi, when calculating using dot product we get
                % the small angle. sometimes we need an angles that is gt
                % 180. to calculate it, rotate the vector (eg, rotate xbody
                % around the stroke plane and check if we got the span)
                % and check if we get the second vector
                ypln = cross(XbodOnPln,strkpln);
                yphi = dot(spnOnPln',ypln');
                xphi = dot(spnOnPln',XbodOnPln');
                if k == 2
                    yphi = -yphi;
                end
                
                phitmp = mod(atan2(yphi,xphi),2*pi);
                
                obj.(wingname{k}).angles.phi = phitmp'*180/pi;
                obj.(wingname{k}).angles.theta = 90 - real(acosd(dot(spn',obj.body.vectors.strkPlan')))';
                calcPsi(obj,spn,obj.(wingname{k}).vectors.(chordProp),strkpln,wingname{k})
            end
        end
        
        
        function calcPsi(obj,spn,chord,strkpln,wingname,varargin)
            parser = inputParser;
            addParameter(parser,'psiname','psi');
            parse(parser, varargin{:})
            psisave = parser.Results.psiname;
            if strfind(wingname,'rightwing')
                Surf_sp=cross(strkpln,spn);
                signy = 1;
            else
                Surf_sp=cross(spn,strkpln);
                signy = -1;
            end
            Surf_sp = Surf_sp./vecnorm(Surf_sp')';
            yax = cross(spn,Surf_sp);
            yax = yax./vecnorm(yax')';
            ypsi = signy*dot(chord',yax');
            xpsi = dot(chord',Surf_sp');
            angtmp = mod(atan2(ypsi,xpsi),2*pi);
            if strfind(psisave,'sec')
                name2save = split(psisave);
                obj.(wingname).angles.(name2save{1}).(name2save{2})= angtmp'*180/pi;
            else
                obj.(wingname).angles.(psisave)=angtmp*180/pi;
            end
            
        end
        
        function cooTE = cutwingBound(obj,vmax,vmin,coall,spanEW,part,flgTE)
            % take only part of the wing (part: number of parts to take).
            % sevide the wing according to the projection on the span.
            delroot = (vmax - vmin)/part;
            vecofsec = linspace(vmin + delroot,vmax,part);
            spanmat = repmat(spanEW',size(coall,1),1);
            projOnSpnall = dot(coall',spanmat');
            if flgTE == 1
                quart_coor =  (projOnSpnall > (vecofsec(end-1)) & projOnSpnall < (vecofsec(end)))';
            else
                quart_coor =  (projOnSpnall > (vecofsec(1)) & projOnSpnall < (vecofsec(end)))';
            end
            cooTE = coall(quart_coor==1,:);
        end
        
        function [coords] =  calcSecChord(obj,LETEname,wingname,hull3d,varargin)
            parser = inputParser;
            addParameter(parser,'plot',0);
            addParameter(parser,'tipprop','projtip');
            addParameter(parser,'chordprop','chord_bound');

            addParameter(parser,'wingsize',30);
            parse(parser, varargin{:})
            signLE = 1;
            chordprop =  parser.Results.chordprop;
            tipprop = parser.Results.tipprop;
            if strcmp(LETEname,'TE')
                signLE = -1;
            end
            for fr = obj.frames
                frm = (fr == obj.frames);
                try
                    bnd =[double(hull3d.(wingname).hull.(LETEname){frm})];
                    
                    % rotate span, stroke plane and chord to EW axes
                    spanEW = obj.rotmat_EWtoL'*obj.(wingname).vectors.span(frm,:)';
                    strkplnEW = obj.rotmat_EWtoL'*obj.body.vectors.strkPlan(frm,:)';
                    chordEW = obj.rotmat_EWtoL'*obj.(wingname).vectors.(chordprop)(frm,:)';
                    if dot(obj.(wingname).coords.tip(frm,:) - mean(bnd),spanEW')< 0
                        obj.(wingname).coords.tip(frm,:) = obj.calculateTip(wingname,double(hull3d.(wingname).hull.hull3d{frm}),mean(double(hull3d.(wingname).hull.wingCM{frm})),'spanframe',find(frm),'tipop',1);
                    end
                    % project LE/TE on span, find min and max value and devide the
                    % wing to sections (numberof_sec)
                    vmax = dot(obj.(wingname).coords.tip(frm,:)',spanEW');
                    vmin = vmax - parser.Results.wingsize;
                    spanmat = repmat(spanEW',size(bnd,1),1);
                    projOnSpn = dot(bnd',spanmat');
                    vmax = max(projOnSpn);
                    vmin = min(projOnSpn);
                    % start from the second section because the root is noisy
                    delroot = (vmax - vmin)/obj.wingana.numofsec;
                    vecofsec = linspace(vmin + delroot,vmax,obj.wingana.numofsec);
                    
                    if strcmp(LETEname,'LE')
                        % use part of the TE and most of the LE as a plane to
                        % calculate the local chord
                        cooTE = cutwingBound(obj,vmax,vmin,double(hull3d.(wingname).hull.TE{frm}),spanEW,obj.wingana.part4TE,1);
                        cooLE = cutwingBound(obj,vmax,vmin,double(hull3d.(wingname).hull.LE{frm}),spanEW,obj.wingana.part4LE,2);
                        coo = [cooTE;cooLE];
                        mean_win = mean([coo;obj.(wingname).coords.tip(frm,:)]);
                        obj.(wingname).vectors.LETEPlane(frm,:) = Functions.affine_fit([coo;obj.(wingname).coords.tip(frm,:)]);
                        projspan = Functions.projection(spanEW',obj.(wingname).vectors.LETEPlane(frm,:))';
                        
                        tpo = obj.(wingname).coords.tip(frm,:);
                        v = tpo-mean_win;
                        d = obj.(wingname).vectors.LETEPlane(frm,:)*v';
                        projectedtip = tpo - d*obj.(wingname).vectors.LETEPlane(frm,:);
                        
                        
                        obj.(wingname).vectors.spanforpsi(frm,:) = -projspan';
                        obj.(wingname).vectors.projtip(frm,:) = projectedtip;
                        
                        if parser.Results.plot == 1
                            allbound = [hull3d.(wingname).hull.LE{frm};hull3d.(wingname).hull.TE{frm}];
                            meanb = mean(allbound);
                            used4span = [coo;obj.(wingname).coords.tip(frm,:)];
                            figure;plot3(allbound(:,1),allbound(:,2),allbound(:,3),'.b','markersize',5);hold on;
                            plot3(used4span(:,1),used4span(:,2),used4span(:,3),'.r','markersize',5)
                            tp = obj.(wingname).vectors.projtip(frm,:);
                            tpori = obj.(wingname).coords.tip(frm,:);
                            plot3(tpori(:,1),tpori(:,2),tpori(:,3),'.k','markersize',10)
                            plot3(tp(:,1),tp(:,2),tp(:,3),'.m','markersize',10)
                            hold on;quiver3(meanb(1),meanb(2),meanb(3),spanEW(1),spanEW(2),spanEW(3),50)
                            axis equal
                            grid on ; box on ;
                            grid minor
                            
                        end
                    else
                        obj.(wingname).vectors.spanforpsi(frm,:) = -obj.rotmat_EWtoL' *obj.(wingname).vectors.spanforpsi(frm,:)';
                    end
                    tip = obj.(wingname).vectors.(tipprop)(frm,:);
                    
                    for k = 1:1:length(vecofsec)-1
                        secname = sprintf('sec%d',k);
                        quart_coor =  (projOnSpn > (vecofsec(k)) & projOnSpn < (vecofsec(k+1)))';
                        coords{k} = bnd(quart_coor==1,:);
                        secCM = mean(coords{k},1);
                        
                        % for each section calculate the direction of vector between the span and the section's CM
                        vec1 = (dot(secCM - tip,obj.(wingname).vectors.spanforpsi(frm,:))*obj.(wingname).vectors.spanforpsi(frm,:));
                        vec2 = secCM - tip;
                        vec4calc = vec2 - vec1;
                        obj.(wingname).vectors.(secname).(LETEname)(frm,:) = (obj.rotmat_EWtoL*(vec4calc/norm(vec4calc))')';
                        
                        % check the chord's direction
                        if dot(signLE * chordEW,strkplnEW) < 0
                            obj.(wingname).vectors.(secname).(LETEname)(frm,:) = -obj.(wingname).vectors.(secname).(LETEname)(frm,:);
                        end
                        
                        if parser.Results.plot == 1
                            CM4plt_1 = dot(secCM,vec1/norm(vec1))*vec1/norm(vec1);
                            CM4plt_2 = (dot(secCM,vec4calc/norm(vec4calc)))*vec4calc/norm(vec4calc);
                            
                            
                            CM4plt = secCM - vec4calc;
                            
                            hull_wing = hull3d.rightwing.hull.hull3d{frm};
                            
                            hull_wing2 = hull3d.leftwing.hull.hull3d{frm};
                            
                            bod = hull3d.body.hull{frm};
                            %                 allBound = [hull3d.(wingname).hull.LE{frm};hull3d.(wingname).hull.TE{frm}];
                            plot3(hull3d.(wingname).hull.LE{frm}(:,1),hull3d.(wingname).hull.LE{frm}(:,2),hull3d.(wingname).hull.LE{frm}(:,3),'k*');hold on;
                            plot3(hull3d.(wingname).hull.TE{frm}(:,1),hull3d.(wingname).hull.TE{frm}(:,2),hull3d.(wingname).hull.TE{frm}(:,3),'b*');hold on;
                            plot3(bod(:,1),bod(:,2),bod(:,3),'.g')
                            %                             plot3(hull_wing(:,1),hull_wing(:,2),hull_wing(:,3),'.c');
                            %                             plot3(hull_wing2(:,1),hull_wing2(:,2),hull_wing2(:,3),'.c');
                            
                            plot3(bnd(:,1),bnd(:,2),bnd(:,3),'*');hold on;
                            plot3(bnd(:,1),bnd(:,2),bnd(:,3),'*');hold on;
                            plot3(coords{k}(:,1),coords{k}(:,2),coords{k}(:,3),'*')
                            %                             quiver3(secCM(1),secCM(2),secCM(3),vec1(1),vec1(2),vec1(3),3,'linewidth',3)
                            %                             quiver3(secCM(1),secCM(2),secCM(3),vec2(1),vec2(2),vec2(3),3)
                            quiver3(CM4plt(1),CM4plt(2),CM4plt(3),vec4calc(1),vec4calc(2),vec4calc(3),3)
                            axis equal
                            grid on ; box on ;
                            grid minor
                        end
                        
                    end
                    obj.(wingname).vectors.spanforpsi(frm,:) = -obj.rotmat_EWtoL*obj.(wingname).vectors.spanforpsi(frm,:)';
                catch
                    for k = 1:1:obj.wingana.numofsec-1
                        secname = sprintf('sec%d',k);
                        obj.(wingname).vectors.(secname).(LETEname)(frm,:) = [nan nan nan];
                        obj.(wingname).vectors.spanforpsi(frm,:) = [nan nan nan];
                        if strcmp(tipprop,'tip') == 0
                            obj.(wingname).vectors.(tipprop)(frm,:) = [nan nan nan];
                        end
                    end
                end
            end
            
            
        end
        
        function [vecCM4plt,vec1,vecCM4plt_wing] =  calcSecChord_frame(obj,LETEname,wingname,hull3d,fr,varargin)
            parser = inputParser;
            addParameter(parser,'plot_sec',0);
            addParameter(parser,'plot',0);
            addParameter(parser,'vecCM4plt',0);
            addParameter(parser,'plotSpn',0);
            
            addParameter(parser,'tipprop','projtip');
            addParameter(parser,'wingsize',30);
            parse(parser, varargin{:})
            vecCM4plt = [];
            signLE = 1;
            tipprop = parser.Results.tipprop;
            if strcmp(LETEname,'TE')
                signLE = -1;
            end
            
            if parser.Results.vecCM4plt ~= 0
                CM4pltvec = parser.Results.vecCM4plt;
            end
            
            frm = (fr == obj.frames);
            try
                bnd =[double(hull3d.(wingname).hull.(LETEname){frm})];
                
                % rotate span, stroke plane and chord to EW axes
                spanEW = obj.rotmat_EWtoL'*obj.(wingname).vectors.span(frm,:)';
                strkplnEW = obj.rotmat_EWtoL'*obj.body.vectors.strkPlan(frm,:)';
                chordEW = obj.rotmat_EWtoL'*obj.(wingname).vectors.chord(frm,:)';
                
                % project LE/TE on span, find min and max value and devide the
                % wing to sections (numberof_sec)
                vmax = dot(obj.(wingname).coords.tip(frm,:)',spanEW');
                vmin = vmax - parser.Results.wingsize;
                spanmat = repmat(spanEW',size(bnd,1),1);
                projOnSpn = dot(bnd',spanmat');
%                                 vmin = min(projOnSpn);
%                                 vmax = max(projOnSpn);

                % start from the second section because the root is noisy
                delroot = (vmax - vmin)/obj.wingana.numofsec;
                vecofsec = linspace(vmin + delroot,vmax,obj.wingana.numofsec);
                
                quart_coor =  (projOnSpn > (vecofsec(1)))';
                LETE4plt = bnd(quart_coor==1,:);
                if strcmp(LETEname,'LE')
                    % use part of the TE and most of the LE as a plane to
                    % calculate the local chord
                    cooTE = cutwingBound(obj,vmax,vmin,double(hull3d.(wingname).hull.TE{frm}),spanEW,obj.wingana.part4TE,1);
                    cooLE = cutwingBound(obj,vmax,vmin,double(hull3d.(wingname).hull.LE{frm}),spanEW,obj.wingana.part4LE,2);
                    coo = [cooTE;cooLE];
                    mean_win = mean([coo;obj.(wingname).coords.tip(frm,:)]);
                    obj.(wingname).vectors.LETEPlane(frm,:) = Functions.affine_fit([coo;obj.(wingname).coords.tip(frm,:)]);
                    projspan = Functions.projection(spanEW',obj.(wingname).vectors.LETEPlane(frm,:))';
                    
                    tpo = obj.(wingname).coords.tip(frm,:);
                    v = tpo-mean_win;
                    d = obj.(wingname).vectors.LETEPlane(frm,:)*v';
                    projectedtip = tpo - d*obj.(wingname).vectors.LETEPlane(frm,:);
                    
                    
                    obj.(wingname).vectors.spanforpsi(frm,:) = -projspan';
                    obj.(wingname).vectors.projtip(frm,:) = projectedtip;
                    
                    if parser.Results.plot == 1
                        allbound = [hull3d.(wingname).hull.LE{frm};hull3d.(wingname).hull.TE{frm}];
                        used4span = [coo;obj.(wingname).coords.tip(frm,:)];
                        figure;plot3(allbound(:,1),allbound(:,2),allbound(:,3),'.b','markersize',5);hold on;
                        plot3(used4span(:,1),used4span(:,2),used4span(:,3),'.r','markersize',5)
                        tp = obj.(wingname).vectors.projtip(frm,:);
                        tpori = obj.(wingname).coords.tip(frm,:);
                        plot3(tpori(:,1),tpori(:,2),tpori(:,3),'.k','markersize',10)
                        plot3(tp(:,1),tp(:,2),tp(:,3),'.m','markersize',10)
                        
                        axis equal
                        grid on ; box on ;
                        grid minor
                        
                    end
                else
                    obj.(wingname).vectors.spanforpsi(frm,:) = -obj.rotmat_EWtoL' *obj.(wingname).vectors.spanforpsi(frm,:)';
                end
                tip = obj.(wingname).vectors.(tipprop)(frm,:);
                
                for k = 1:1:length(vecofsec)-1
                    secname = sprintf('sec%d',k);
                    quart_coor =  (projOnSpn >= (vecofsec(k)) & projOnSpn <= (vecofsec(k+1)))';
                    coords{k} = bnd(quart_coor==1,:);
                    secCM = mean(coords{k});
                    
                    % for each section calculate the direction of vector between the span and the section's CM
                    vec1 = (dot(secCM - tip,obj.(wingname).vectors.spanforpsi(frm,:))*obj.(wingname).vectors.spanforpsi(frm,:));
                    vec2 = secCM - tip;
                    vec4calc = vec2 - vec1;
                    obj.(wingname).vectors.(secname).(LETEname)(frm,:) = (obj.rotmat_EWtoL*(vec4calc/norm(vec4calc))')';
                    
                    % check the chord's direction
                    if dot(signLE * chordEW,strkplnEW) < 0
                        obj.(wingname).vectors.(secname).(LETEname)(frm,:) = -obj.(wingname).vectors.(secname).(LETEname)(frm,:);
                    end
                    
                    if parser.Results.plot_sec == 1
                        ind2mm = 50*10^-6*1000;
                        bod = (ind2mm*obj.rotmat_EWtoL*(double(hull3d.body.hull{frm})'))';
                        meanbod = mean(bod);
                        bod = double(bod) - meanbod;
                        
                        if parser.Results.vecCM4plt == 0
                            CM4plt = secCM - vec4calc;
                        else
                            CM4plt = CM4pltvec(k,:);
                        end
                        
                        
                        secCo = (ind2mm*obj.rotmat_EWtoL*(coords{k})')'- meanbod;
                        
                        %                             Xwing = vec1/norm(vec1);
                        %                             d = dot(repmat(Xwing',1,size(secCo,1)),secCo');
                        %                             projOnX = d'*Xwing;
                        
                        
                        
                        
                        vec4calc = (obj.rotmat_EWtoL*vec4calc')';
                        vec1 = -(obj.rotmat_EWtoL*vec1')';
                        
                        
                        CMpltLab = ind2mm*(obj.rotmat_EWtoL*(CM4plt )')'- meanbod;
                        
                        LE = (ind2mm*obj.rotmat_EWtoL*(double(hull3d.(wingname).hull.LE{frm})'))' - meanbod;
                        TE = (ind2mm*obj.rotmat_EWtoL*(double(hull3d.(wingname).hull.TE{frm})'))' - meanbod;
                        
                        
                        %                 allBound = [hull3d.(wingname).hull.LE{frm};hull3d.(wingname).hull.TE{frm}];
                        plot3(LE(:,1),LE(:,2),LE(:,3),'.r','markersize',8);hold on;
                        plot3(TE(:,1),TE(:,2),TE(:,3),'.b','markersize',8);hold on;
                        plot3(bod(:,1),bod(:,2),bod(:,3),'.g')
                        
                        
                        Xwing = vec1/norm(vec1);
                        Ywing = vec4calc/norm(vec4calc);
                        Zwing = cross(Xwing,Ywing);
                        Zwing = Zwing/norm(Zwing);
                        
                        
                        
                        d = dot(repmat(Xwing',1,size(secCo,1)),secCo');
                        projOnX = d'*Xwing;
                        cmproj = mean(projOnX);
                        
                        bot = projOnX + CMpltLab - cmproj;
                        top = bot + vec4calc*ind2mm;
                        
                        
                        d = dot(repmat(Xwing',1,size(bot,1)),bot');
                        d2 = dot(repmat(Xwing',1,size(top,1)),top');
                        
                        [~,ind] = min(d);
                        [~,ind2] = max(d);
                        
                        
                        pltrec = [bot(ind,:);bot(ind2,:);top(ind2,:);top(ind,:)];
                        clor = colormap('colorcube');
                        c = 43;
                        if strcmp(LETEname,'TE')
                            c = 151;
                        end
                        hold on;s = fill3(pltrec(:,1),pltrec(:,2),pltrec(:,3),clor(k+10 + c ,:))
                        s.FaceAlpha = 0.4;
                        
                        
                        %                             seconpsi = secCo - CM4plt;
                        %                             Y = dot(secCo',repmat(Xwing',1,size(secCo,1)))
                        %                             X = dot(secCo',repmat(Ywing',1,size(secCo,1)))
                        %                             Z = dot(secCo',repmat(Zwing',1,size(secCo,1)))
                        %                             xyz = [X',Y',Z'] ;
                        %                             plot3(projOnX(:,1),projOnX(:,2),projOnX(:,3),'.')
                        %                             plot3(secCo(:,1),secCo(:,2),secCo(:,3),'.')
                        %
                        vec4calc = vec4calc/norm(vec4calc);
                        vec1 = vec1/norm(vec1);
                        quiver3(CMpltLab(1),CMpltLab(2),CMpltLab(3),vec4calc(1),vec4calc(2),vec4calc(3),10*ind2mm,'color','k')
                        if k == 1
                            quiver3(CMpltLab(1),CMpltLab(2),CMpltLab(3),vec1(1),vec1(2),vec1(3),25*ind2mm,'color','k')
                        end
                        axis equal
                        grid on ; box on ;
                        grid minor
                        vecCM4plt(k,1:3) = CM4plt;
                        vecCM4plt_wing(k,1:3) = CMpltLab;
                        
                        
                        
                        
                        
                        
                        
                        
                        %                             X = sqr(:,1);
                        %                             Y = sqr(:,2);
                        %                             Z = -1/pln(3)*(pln(1)*X + pln(2)*Y)
                        %
                        %                             plot3(X,Y,Z,'*')
                        %                             p =fill3(X,Y,Z,'r')
                        %                             p.FaceAlpha = 0.2;
                        %                             quiver3(0,0,0,pln(1),pln(2),pln(3),0.5,'linewidth',1,'color','m');axis equal
                        
                    end
                    
                end
                obj.(wingname).vectors.spanforpsi(frm,:) = -obj.rotmat_EWtoL*obj.(wingname).vectors.spanforpsi(frm,:)';
            catch
                for k = 1:1:obj.wingana.numofsec-1
                    secname = sprintf('sec%d',k);
                    obj.(wingname).vectors.(secname).(LETEname)(frm,:) = [nan nan nan];
                    obj.(wingname).vectors.spanforpsi(frm,:) = [nan nan nan];
                    if strcmp(tipprop,'tip') == 0
                        obj.(wingname).vectors.(tipprop)(frm,:) = [nan nan nan];
                    end
                end
            end
        end
        
        
        
        function calcSecPsi(obj,LETEname,wingname)
            for k = 1:1:obj.wingana.numofsec - 1
                psisec_name = sprintf('sec%d %s',k,LETEname);
                name = split(psisec_name);
                obj.calcPsi(obj.(wingname).vectors.span,obj.(wingname).vectors.(name{1}).(name{2}),obj.body.vectors.strkPlan,wingname,'psiname',psisec_name);
                avpsi_sec_dist = abs(obj.(wingname).angles.(name{1}).(name{2}) - obj.(wingname).angles.psi);
                flipind = abs(avpsi_sec_dist) > abs(avpsi_sec_dist-180);
%                 if isempty(flipind) == 0
%                     obj.(wingname).angles.(name{1}).(name{2})(flipind) = 180-obj.(wingname).angles.(name{1}).(name{2})(flipind);
%                 end
            end
        end
        
        function plotPsi_sec(obj,wingname,varargin)
            parser = inputParser;
            addParameter(parser,'stenfr',[0 0]); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            LETEname = {'LE','TE'};
                        if sum(parser.Results.stenfr) ~= 0

            fr0 = find(obj.frames == parser.Results.stenfr(1));
            fr1 = find(obj.frames == parser.Results.stenfr(2));
                        else
                            fr0 = 1;
                            fr1 = length(obj.frames)
                        end
            for kLETE = 1:1:2
                a(kLETE) = subplot(2,1,kLETE)
                for k = 1:1:obj.wingana.numofsec - 1
                    psisec_name = sprintf('sec%d %s',k,LETEname{kLETE});
                    name = split(psisec_name);
                    hold on ;plot(obj.video.timeframe(fr0:fr1),obj.(wingname).angles.(name{1}).(name{2})(fr0:fr1),'-*');grid on;
                    title(LETEname{kLETE})
                end
            end
            linkaxes([a(1),a(2)])
            sgtitle(wingname);
        end
        
        
        
        function [wingIm,DownIm,UpIm] = Split2Dbound(obj,TwoD,SE)
            % Generate the boundary of the projected hull and split it.
            % make sure there are no holes in the projected image by
            % closeing it (SE)
            % TwoD contains the projection of LE,TE,wing_cut,Tip,Wing
            % CM,body.
            
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
            Up=wingAll(wingAll(:,2)>=yperp_wall,:);
            
            [wingIm] = Functions.ImfromSp(obj.cameras.all.size_image,fliplr(TwoD{3}));
            wingIm=imclose(wingIm,SE);
            
            [DownIm] = Functions.ImfromSp(obj.cameras.all.size_image,fliplr(Down));
            DownIm=imclose(DownIm,SE);
            
            [UpIm] = Functions.ImfromSp(obj.cameras.all.size_image,fliplr(Up));
            UpIm=imclose(UpIm,SE);
        end
    end
    methods (Access = private)
        
        function [idxdwnStrk] = FindUp_downStrk(~,changeSgn,mean_strks,up_down_strk)
            % find where the signs change (compared to 0 on the body axis).
            % in this position the wing are farthest apart.
            
            downstrk=find((changeSgn(1,1:end-1)+changeSgn(2,2:end))==up_down_strk);
            mean_val=[mean_strks(downstrk+1);mean_strks(downstrk)];
            [~,indMin]=min(abs(mean_val),[],1);
            idx_vec=[downstrk+1;downstrk];
            idxdwnStrk=(idx_vec([(indMin==1);(indMin==2)]));
        end
        
        function [idx4StrkPln] = Choose_GrAng_wing1_wing2(~,distSpans,FrbckStrk,angleTH,Search_cell)
            % make sure the chosen point are indeed the greatest. Check +/- 10 cells
            % around the chosen index
            
            idx4StrkPln = FrbckStrk(distSpans(FrbckStrk)>angleTH);
            
            for kind_strk=1:1:length(idx4StrkPln)
                inds2ch=idx4StrkPln(kind_strk)-Search_cell:idx4StrkPln(kind_strk)+Search_cell;
                inds2ch(inds2ch<=0)=[];
                inds2ch(inds2ch>length(distSpans))=[];
                [val_max indMax]=max(distSpans(inds2ch));
                idx4StrkPln(kind_strk)=inds2ch(indMax);
            end
        end
        
        
    end
    
end

