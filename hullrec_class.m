classdef hullrec_class<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cameras
        im4hull
        seed
        Find_vol
        voxelSize
        volLength
        voxelSize4search
        real_coord
        index_seed
        ofst
        parts2run
        framehull
        TseviMethod
        ind_0
        real_3d
        minCameraBody
        openBodySEsize
        voxelSize_body
        real_volume
        real_volume_frame
        ind_0_frame
    end
    
    methods
        function obj = hullrec_class(easyWandData,varargin)
            parser = inputParser;
            addParameter(parser,'ZaxCam',1); % number of camera pointing in Z lab axis
            addParameter(parser,'camvec',[1,2,3,4]); % order of cameras in easywang file
            % x,y,zdeg : rotating the volume to get maximal covverege of
            % voxels. User can change their values.
            addParameter(parser,'xdeg',0);
            addParameter(parser,'ydeg',0);
            addParameter(parser,'zdeg',0);
            % parameters for hull reconstruction
            addParameter(parser,'VxlSize',50e-6); % voxel size
            addParameter(parser,'voxelSize_body',50e-6); % voxel size
            addParameter(parser,'VolLen',5e-3); % volume length
            % size of voxels to perform initial search to find the
            % volume for hull. use if no older volume was used as an
            % input for hull_reconstruction. To reconstruct a dense
            % hull without a pre defined volume set VxlSize4search as 0
            % to use the same value of VxlSize or define VxlSize4search
            % as any number.
            addParameter(parser,'VxlSize4search',40e-5);
            addParameter(parser,'ofst',2*40e-5); % size of voxels to perform initial search to find the volume for hull
            addParameter(parser,'minCameraBody',1);
            addParameter(parser,'TseviMethod',0);
            addParameter(parser,'openBodySEsize',2);
            parse(parser, varargin{:})
            
            obj.openBodySEsize = parser.Results.openBodySEsize;
            obj.minCameraBody = parser.Results.minCameraBody;
            obj.TseviMethod = parser.Results.TseviMethod;
            obj.ofst = parser.Results.ofst;
            obj.voxelSize = parser.Results.VxlSize;
            obj.volLength = parser.Results.VolLen;
            obj.voxelSize4search = parser.Results.VxlSize4search;
            if parser.Results.VxlSize4search == 0
                obj.voxelSize4search = parser.Results.VxlSize;
            end
            obj.cameras.all.camvec = parser.Results.camvec;
            obj.cameras.all.size_image=[easyWandData.imageHeight(1),easyWandData.imageWidth(1)];
            centers_cam=squeeze(easyWandData.DLTtranslationVector);
            obj.cameras.all.ZaxCam = parser.Results.ZaxCam;
            indcam = 1;
            for kcam = sort(obj.cameras.all.camvec)
                % get all camera parameters. input camvec as easywand order (if
                % the first camera is 4 camvec should be [4,1,2,3]. the
                % camera parameters will be orgenized as [1,2,3,4]
                dltloc = find(obj.cameras.all.camvec == kcam);
                obj.cameras.cams_arr(indcam).camera_cnt=easyWandData.DLTtranslationVector(:,:,dltloc);
                obj.cameras.cams_arr(indcam).dlt=easyWandData.coefs(:,dltloc);
                obj.cameras.cams_arr(indcam).reshaped_dlt=reshape([easyWandData.coefs(:,dltloc);1],[4,3])';
                obj.cameras.cams_arr(indcam).invDLT=pinv(obj.cameras.cams_arr(indcam).reshaped_dlt);
                obj.cameras.cams_arr(indcam).camera_dir=-easyWandData.rotationMatrices(3,1:3,dltloc)';
                obj.cameras.all.DLT_coefs(:,indcam)=easyWandData.coefs(:,dltloc);
                obj.cameras.all.centers_cam(:,indcam)=centers_cam(:,dltloc);
                indcam = indcam+ 1;
            end
            
            obj.cameras.all.numOfCams = length(obj.cameras.all.camvec);
            obj.cameras.all.camvec = [1:1:length(obj.cameras.all.camvec)];%sort(obj.cameras.all.camvec); % save new cameravec order
            % find the direction of Z lab in camera axis
            if parser.Results.ZaxCam ~=0
                % if a camera is pointing to Z lab
                % update the index of the camera according to the new order
                % of camvec (if camvec = [4,1,2,3] and Zcam = 4, the new
                % vectors will be: camvec = [1,2,3,4] and Zcam = 1)
                obj.cameras.all.ZaxCam = find(obj.cameras.all.camvec == obj.cameras.all.ZaxCam);
                cameras_up =  obj.cameras.cams_arr(obj.cameras.all.ZaxCam).camera_dir;
            else
                
                cameras_up=zeros(3,1);
                for kcam=sort(obj.cameras.all.camvec)
                    % if the configuration is cameras that have a mean
                    % direction of Z lab and no camera pointing to Z lab
                    cameras_up=cameras_up+obj.cameras.cams_arr(kcam).camera_dir;
                end
            end
            if min(obj.cameras.all.camvec) ~= 1
                obj.cameras.all.camvec = [1:1:length(obj.cameras.all.camvec)];
            end
            obj.cameras.all.Rotation_Matrix=vrrotvec2mat(vrrotvec(cameras_up/norm(cameras_up),[0,0,1])); % calculate the rotation matrix between camera and lab axis
            obj.cameras.Fundamental_matrices=get_fundamental_mats(obj);
            Rx = rotx(parser.Results.xdeg);Ry = roty(parser.Results.ydeg);Rz = rotz(parser.Results.zdeg);
            obj.cameras.all.RotMat_vol = Rz*Ry*Rx;
        end
        
        function [Fs,couples] = get_fundamental_mats(obj)
            % Description:
            % returns array of fundamental matrices with order F12,F13,...,F23,...
            %
            % Required input:
            % obj - all_cameras_class
            %
            % Output:
            % Fs - 3*3*number_of_couples array with fundamental matrix ordered -
            % F12,F13,...,F23,...
            % couples - camera indices for each fundamental matrix
            num_of_cams=length(obj.cameras.cams_arr);
            couples=nchoosek(1:num_of_cams,2);
            
            DLT_mats = zeros(3,4,num_of_cams);
            DLT_invs = zeros(4,3,num_of_cams);
            cam_cents = zeros(3,num_of_cams);
            for cam_ind=sort(obj.cameras.all.camvec)
                DLT_mats(:,:,cam_ind)=obj.cameras.cams_arr(cam_ind).reshaped_dlt;
                DLT_invs(:,:,cam_ind)=obj.cameras.cams_arr(cam_ind).invDLT;
                cam_cents(:,cam_ind)=obj.cameras.cams_arr(cam_ind).camera_cnt;
            end
            
            Fs = zeros(3,3,size(couples,1));
            % calculate fundamental matrices
            for couple_ind=1:size(couples,1)
                A = DLT_mats(:,:,couples(couple_ind,2))*[cam_cents(:,couples(couple_ind,1));1];
                C = [0 -A(3) A(2); A(3) 0 -A(1); -A(2) A(1) 0];% skew-symmetric matrix
                Fs(:,:,couple_ind)=C*DLT_mats(:,:,couples(couple_ind,2))*DLT_invs(:,:,couples(couple_ind,1));
            end
        end
        
        function Seg2hullRIm(obj,seg2D,fr)
            % save current 2D images in hull class
            for kcam = 1:1:obj.cameras.all.numOfCams
                obj.im4hull.sprs.wing1{kcam} = seg2D.wing1{kcam}(fr).indIm;
                obj.im4hull.sprs.wing2{kcam} = seg2D.wing2{kcam}(fr).indIm;
                obj.im4hull.sprs.body{kcam} = seg2D.body{kcam}(fr).indIm;
                body4open = Functions.ImfromSp(obj.cameras.all.size_image,obj.im4hull.sprs.body{kcam});
                
                %                 bodyThik = bwmorph(body4open,'thicken',1);
                %                 [ro co] = find(bodyThik);
                %                 bodyThik = [ro co ones(length(ro),1)];
                
                
                SE = strel('disk',obj.openBodySEsize);
                body4open = imerode(body4open,SE);
                [ro co] = find(body4open);
                body4open = [ro co ones(length(ro),1)];
                obj.im4hull.sprs.body3{kcam} = body4open;
                
                
                
                %                 obj.im4hull.sprs.all{kcam} = [seg2D.wing1{kcam}(fr).indIm;seg2D.wing2{kcam}(fr).indIm;bodyThik];%seg2D.all{kcam}(fr).indIm;
                obj.im4hull.sprs.all{kcam} = seg2D.all{kcam}(fr).indIm;
            end
            
        end
        
        function FindSeed(obj,prop)
            % Description:
            % Constructor
            %
            % Required input:
            % prop - name of the image to hull (string)
            % optional : 'diffIm' - 1 prop is a cell array of the images
            % names. used when the images are different
            
            kcamInd = 1;
            for kcam = sort(obj.cameras.all.camvec)
                % find the blob's center of mass
                
                % find the 2D center of mass
                x = obj.im4hull.sprs.(prop){kcam}(:,2);
                y = obj.im4hull.sprs.(prop){kcam}(:,1);
                
                if kcam == (obj.cameras.all.ZaxCam) && obj.cameras.all.ZaxCam ~=0
                    y = 801-y;
                end
                
                CM=round([mean(x),mean(y)]);
                
                %-------------------------
                
                % calculate the 3d position of the cm
                centerPt=[CM(1);...
                    obj.cameras.all.size_image(1)+1-CM(2);1];
                PB(kcamInd,:)=obj.cameras.cams_arr(kcam).invDLT*centerPt;
                kcamInd = kcamInd+1;
            end
            % generate a line for each 3d point and find the point closest
            % to all lines - this is the seed
            ends=PB(:,1:3)./(PB(:,4));
            if ~isnan(ends)
                obj.seed=Functions.lineIntersect3D(obj.cameras.all.centers_cam(:,obj.cameras.all.camvec)',ends);
            else
                obj.seed=[nan nan nan];
            end
            
        end
        
        function obj=hull_params(obj)
            % Description:
            % Constructor
            %
            % Required input:
            % seed - 3d point to start the reconstruction (asstimated center of mass)
            % voxelSize - size of each voxel in meters
            % volLength  - size of the square sub-vol cube to reconstruct (meters)
            
            
            obj.seed=round(obj.seed/obj.voxelSize) * obj.voxelSize; % allign to grid
            % create a vector of the real coordinates of the grid
            obj.real_coord = cell2mat(arrayfun(@(x) colon(x,obj.voxelSize,x+...
                obj.volLength),obj.seed'-obj.volLength/2,'UniformOutput',false));
            [~ ,I]=min(abs(obj.real_coord-obj.seed'),[],2);
            obj.index_seed=I;
            volNoRot = [obj.seed(1)-obj.volLength/2,obj.seed(2)-obj.volLength/2,obj.seed(3)-obj.volLength/2;...
                obj.seed(1)+obj.volLength/2,obj.seed(2)+obj.volLength/2,obj.seed(3)+obj.volLength/2];
            obj.Find_vol = (volNoRot')';
        end
        function [ind_0,real_volume] =  createVol(obj,vxlSz)
            % create a low density square volume
            x = [-obj.volLength/2 :vxlSz:obj.volLength/2];
            y = [-obj.volLength/2 :vxlSz:obj.volLength/2];
            z = [-obj.volLength/2 :vxlSz:obj.volLength/2];
            
            [X,Y,Z] = meshgrid(x',y',z');
            real_vol = [X(:),Y(:),Z(:)];
            [X_ind,Y_ind,Z_ind] = meshgrid([1:1:length(x)]',[1:1:length(y)]',[1:1:length(z)]');
            ind_0 = uint16([X_ind(:),Y_ind(:),Z_ind(:)])';
            obj.real_3d = {x,y,z};
            
            
            
            
            %             for k = obj.cameras.all.camvec
            %                 rot2cam = obj.cameras.all.Rotation_Matrix*obj.cameras.cams_arr(k).camera_dir;
            %                 r = vrrotvec2mat(vrrotvec([0,0,1],rot2cam));
            %                 rotated_cube{k} = (r*real_vol')';
            %             end
            %             newvol = cell2mat(rotated_cube');
            %
            %             % find the voxels that are inside the low density volume
            %             xyz = [newvol(:,1),newvol(:,2),newvol(:,3)];
            %             [indsx,X1] = discretize(xyz(:,1),[min(xyz(:,1))- vxlSz:vxlSz:max(xyz(:,1))+ vxlSz]);
            %             [indsy,Y1] = discretize(xyz(:,2),[min(xyz(:,2)) - vxlSz:vxlSz:max(xyz(:,2)) + vxlSz]);
            %             [indsz,Z1] = discretize(xyz(:,3),[min(xyz(:,3)) - vxlSz:vxlSz:max(xyz(:,3))+ vxlSz]);
            %             ind_0 = uint16(unique([indsx,indsy,indsz],'rows')');
            %             [res,~] = Functions.countRepRow([indsx,indsy,indsz]);
            % %             ind_0 = uint16(res(res(:,4) ==4 ,1:3)');
            %
            %
            %
            %
            %
            %             obj.real_3d = {X1,Y1,Z1};
            real_volume = ([obj.real_3d{1}(ind_0(1,:))',obj.real_3d{2}(ind_0(2,:))',obj.real_3d{3}(ind_0(3,:))']')';
            
        end
        
        function [hull_recon,ind0_frame] = hull_reconstruction(obj,field_string,createvol,ind0_frame,varargin)
            % Description:
            % generates hull reconstruction from images loaded to
            % obj.im4hull.
            % input: field_string - name of images to generate hull from
            %        createvol - 1 when constructing more than one image,
            %        use a set of images to reconstruct a smaller volume.
            %        the reconstructed volume is saved in the property
            %        real_volume_frame and ind_0_frame
            
            
           
            parser = inputParser;
            addParameter(parser,'plotflg',0);
            addParameter(parser,'plot2Dflg',0);
            parse(parser, varargin{:})
            
             hull_recon = [];
            camvec = obj.cameras.all.camvec;
            % initilize hull and boundary cells
            hull_cell = cell(1,length(camvec)+1);hull_indcell = cell(1,length(camvec)+1);
            bb_real = cell(1,length(camvec)+1);bb_inds = cell(1,length(camvec)+1);
            
            % if createvol = 1, use the current images to create a volume
            % used by another image. when the createvol = 0, the volume
            % generated for the image before is used (when createvol = 1).
            if createvol == 1
                % move the volume to the location of the seed. (CM of all
                % images) and update the real_coords by the seed.
                obj.real_coord = {obj.real_3d{1} + obj.seed(1),obj.real_3d{2} + obj.seed(2),obj.real_3d{3} + obj.seed(3)};
                hull_cell{1}= ([obj.real_coord{1}(ind0_frame(1,:))',obj.real_coord{2}(ind0_frame(2,:))',obj.real_coord{3}(ind0_frame(3,:))']')';
                hull_indcell{1}=ind0_frame'; % assign the first cell the complete volume [indices]
                
            else
                % update ind_0 and real3d
                hull_cell{1}= ([obj.real_coord{1}(ind0_frame(1,:))',obj.real_coord{2}(ind0_frame(2,:))',obj.real_coord{3}(ind0_frame(3,:))']')';
                hull_indcell{1}=ind0_frame';%obj.ind_0_frame'; % assign the first cell the  generated for createvol = 1[indices]
            end
            
            Cmpstr1=strcmp(field_string{1},field_string);
            if length(obj.cameras.all.camvec)>2 && (sum(Cmpstr1==1)==2 || sum(Cmpstr1==0)==2)
                if sum(Cmpstr1==1)==2
                    inicam=find(Cmpstr1==0);
                else
                    inicam=find(Cmpstr1==1);
                end
                delCam = camvec(inicam);
                camvec(inicam)=[];
                camvec=[delCam,camvec];
                
                delfs = field_string(inicam);
                field_string(inicam) = [];
                field_string = [delfs,field_string];
            end
            camvec=[1,camvec];
            kcamInd = 1;
            
            
            for kcam=camvec(2:end)
                % for each camera reproject the volume and find the
                % intersecting 2D coordinates. The coodinates represent the
                % voxels that are defined by the camera
                % first camera: use the hole volume
                %           a. reproject the volume and find the pixels
                %           intersecting with the image.
                %               a.1. define a BB in 2D and delete every voxel
                %                   that coresponds to a pixel that is outside of the BB
                %                   save the voxels that are inside the BB
                %                   to use for the next camera.
                %               b.1. search the remaining pixels for intersection
                %                   and delete the coresponding voxels
                % for the next cameras use the vxels in the BB of the
                % previous cameras
                camvectmp = camvec(2:end);
                if kcamInd == 1
                    % for the first camera, use the full volume to
                    % reproject
                    pixel_cam1 = Functions.dlt_inverse(obj.cameras.cams_arr(kcam).dlt,hull_cell{1} );
                else
                    % for the rest of the cameras use the voxels that are
                    % inside the bou
                    camvectmp(find(camvec(2:end)==kcam));
                    pixel_cam1 = Functions.dlt_inverse(obj.cameras.cams_arr(kcam).dlt,bb_real{(camvectmp(find(camvec(2:end)==oldcam)))+1} );
                end
                pixel_cam1=[pixel_cam1(:,1),obj.cameras.all.size_image(1)+1-pixel_cam1(:,2)];
                
                % get the 2D coordinates from the image
                if isempty(obj.im4hull.sprs.(field_string{kcamInd}){kcam})
                    bb_real{camvectmp(kcamInd)+1}=hull_cell{1};
                    bb_inds{camvectmp(kcamInd)+1}=hull_indcell{1};
                    continue
                end
                im1Coords = [obj.im4hull.sprs.(field_string{kcamInd}){kcam}(:,[2,1])];
                
                % flip the image if its a Z axes camera (because of the
                % mirror)
                if kcam == obj.cameras.all.ZaxCam && obj.cameras.all.ZaxCam ~=0
                    im1Coords = [im1Coords(:,1) 801 - im1Coords(:,2)];
                end
                
                % delete every voxel that have apixel outside the bounding box of the image
                
                min_x=max((min(im1Coords(:,1))) - 1,0);
                min_y=max((min(im1Coords(:,2)))-1,0);
                max_x=min((max(im1Coords(:,1)))+1,obj.cameras.all.size_image(2));
                max_y=min((max(im1Coords(:,2)))+1,obj.cameras.all.size_image(1));
                
                fnd3dHull=uint16(pixel_cam1);
                rows2del=(fnd3dHull(:,1)<min_x)*1+(fnd3dHull(:,1)>max_x)*1+...
                    (fnd3dHull(:,2)>max_y)*1+(fnd3dHull(:,2)<min_y)*1;
                if kcamInd == 1
                    bb_real{camvectmp(kcamInd)+1}=hull_cell{1}(rows2del==0,:);
                    bb_inds{camvectmp(kcamInd)+1}=hull_indcell{1}(rows2del==0,:);
                else
                    bb_real{camvectmp(kcamInd)+1}=bb_real{oldcam+1}(rows2del==0,:);
                    bb_inds{camvectmp(kcamInd)+1}=bb_inds{oldcam+1}(rows2del==0,:);
                end
                
                
                % delete the rows outside of the boundary box
                fnd3dHull=fnd3dHull(rows2del==0,:) ;
                indRows=ismember(fnd3dHull,im1Coords,'rows');
                
                % save the voxels generated by the 2d intersection to
                % hull_cell. the cel contains the hull created by all of the camera
                hull_cell{camvectmp(kcamInd)+1}=(bb_real{camvectmp(kcamInd)+1}(indRows,:)')';
                hull_indcell{camvectmp(kcamInd)+1}=bb_inds{camvectmp(kcamInd)+1}(indRows,:);
                
                
                if parser.Results.plot2Dflg == 1
                    figure(1); p(kcamInd) = plot(pixel_cam1(:,1),801-pixel_cam1(:,2),'o');hold on
                    plot(im1Coords(:,1),im1Coords(:,2),'k*');grid on;xlabel('index');ylabel('index');hold on;
                    plot(fnd3dHull(:,1),fnd3dHull(:,2),'mx');grid on;xlabel('index');ylabel('index');hold on;
                    plot(fnd3dHull(indRows,1),fnd3dHull(indRows,2),'^');hold on
                end
                
                kcamInd = kcamInd+1;
                oldcam = kcam;
            end
            % count voxels that apear in all 4 cameras
            checkifall_cam_ind = cell2mat(hull_indcell(2:end)');
            checkifall_cam_real = cell2mat(hull_cell(2:end)');
            
            [~,maph] = Functions.countRepRow(checkifall_cam_ind);
            
            if createvol == 0
                % save the 3D hull as the voxels that apear in all cameras
                hull_recon=unique(checkifall_cam_ind(maph == 4,:),'rows');
            elseif createvol == 2
                ind0_frame=unique(checkifall_cam_ind(maph == 4,:),'rows')';
            else
                % create a new and denser volume. use the minimum value of the voxels in every axis ( +-
                % ofset ) to generate a new 3d bounding box. the new voxels
                % will be used to generate the rest of the images in this
                % frame
                hull_recon_real=unique(checkifall_cam_real(maph == 4,:),'rows');
                obj.real_volume_frame = hull_recon_real; % assign the first cell the complete volume [lab coordinates]
                
                x = [min(hull_recon_real(:,1)) - obj.ofst*obj.voxelSize4search:obj.voxelSize:max(hull_recon_real(:,1))+ obj.ofst*obj.voxelSize4search];
                y = [min(hull_recon_real(:,2)) - obj.ofst*obj.voxelSize4search:obj.voxelSize:max(hull_recon_real(:,2))+ obj.ofst*obj.voxelSize4search];
                z = [min(hull_recon_real(:,3)) - obj.ofst*obj.voxelSize4search:obj.voxelSize:max(hull_recon_real(:,3))+ obj.ofst*obj.voxelSize4search];
                
                [X_ind,Y_ind,Z_ind] = meshgrid([1:1:length(x)]',[1:1:length(y)]',[1:1:length(z)]');
                ind0_frame= uint16([X_ind(:),Y_ind(:),Z_ind(:)])';
                obj.real_coord = {x,y,z};
            end
            
            
            plt = 0;
            if plt == 1
                if createvol == 0
                    hullreal= ([obj.real_coord{1}(hull_recon(:,1))',obj.real_coord{2}(hull_recon(:,2))',obj.real_coord{3}(hull_recon(:,3))']')';
                elseif createvol == 1
                    hullreal= ([obj.real_coord{1}(ind0_frame(1,:))',obj.real_coord{2}(ind0_frame(2,:))',obj.real_coord{3}(ind0_frame(3,:))']')';
                end
                pixel_cam1 = Functions.dlt_inverse(obj.cameras.cams_arr(1).dlt,hullreal );
                pixel_cam2 = Functions.dlt_inverse(obj.cameras.cams_arr(2).dlt,hullreal );
                pixel_cam3 = Functions.dlt_inverse(obj.cameras.cams_arr(3).dlt,hullreal);
                pixel_cam4 = Functions.dlt_inverse(obj.cameras.cams_arr(4).dlt,hullreal);
                
                
                subplot(2,2,1);scatter(pixel_cam1(:,2),pixel_cam1(:,1),'.');hold on;axis equal
                scatter(obj.im4hull.sprs.all{1}(:,1),obj.im4hull.sprs.all{1}(:,2),'x');axis equal
                subplot(2,2,2);scatter(pixel_cam2(:,1),pixel_cam2(:,2),'.');hold on;axis equal
                scatter(obj.im4hull.sprs.all{2}(:,2),801-obj.im4hull.sprs.all{2}(:,1),'x');axis equal
                
                subplot(2,2,3);scatter(pixel_cam3(:,1),pixel_cam3(:,2),'.');hold on;axis equal
                scatter(obj.im4hull.sprs.all{3}(:,2),801-obj.im4hull.sprs.all{3}(:,1),'x');axis equal
                subplot(2,2,4);scatter(pixel_cam4(:,1),pixel_cam4(:,2),'.');hold on;axis equal
                scatter(obj.im4hull.sprs.all{4}(:,2),801-obj.im4hull.sprs.all{4}(:,1),'x');axis equal
            end
            
            if parser.Results.plot2Dflg == 1
                indx = 1;
                for camleg = camvec(2:end)
                    legc{indx} = sprintf('Camera %d',camleg)
                    indx = indx+1;
                end
                legend(p,legc)
            end
            
            if parser.Results.plotflg == 1
                figure;pcshow(hull_recon_real,'r');
            end
            
        end
        
        function  genNames4rec(obj,part_names)
            
            
            % generate matrix of names of parts to reconstruct.
            % all*4 (used to initilize volume)
            % body*4
            %wing1/wing2/body *1 + all*3 (all combinations)
            lenName = length(obj.cameras.all.camvec);
            obj.parts2run = repmat(part_names(1),1,lenName + 1);
            obj.parts2run(2,:) = repmat({'body3'},1,lenName + 1);
            obj.parts2run(2,lenName + 1) = {'body3'};
            Wop = repmat(part_names(1),lenName,lenName + 1);
            for wrl = 2:length(part_names)
                Wop((eye(lenName)==1)) = part_names(wrl);
                Wop(:,end) = [repmat(part_names(wrl),lenName,1)];
                obj.parts2run = [obj.parts2run;Wop];
            end
            if obj.TseviMethod == 1
                % used to compare to Tsevis method (segment wings only in
                % top view)
                indupcam = obj.cameras.all.ZaxCam;
                obj.parts2run = [obj.parts2run(1:6,:);obj.parts2run(indupcam + 6,:);obj.parts2run(indupcam + 10,:)];
            end
        end
        
        function substrBody(obj,varargin)
            parser = inputParser;
            addParameter(parser,'cutFromBod',1);
            addParameter(parser,'uglywing',0);
            parse(parser, varargin{:})
            
            bodyUgly = cell2mat(obj.framehull.body');
            wingsUgly = cell2mat([obj.framehull.wing1';obj.framehull.wing2']);
            [~,bodyCount] = Functions.countRepRow(bodyUgly);
            [~,wingsCount] = Functions.countRepRow(wingsUgly);
            
            if obj.TseviMethod == 0
                wing4 = wingsUgly(wingsCount > 3,1:3);
            else
                wing3 = wingsUgly;
            end
            body2 = bodyUgly(bodyCount > obj.minCameraBody,:);
            body4plot = unique(bodyUgly(bodyCount > (obj.cameras.all.numOfCams - 1),:),'rows');
            
            
            
            
            % keep as wing the voxels that do not intersect with body amd ,ake
            % sure it contains wing3
            %             wings_all_noBod=1-[ismembertol(double(wingsUgly),double(body2),0.07,'ByRows',true)];
            wings_all_noBod=1-[ismember(double(wingsUgly),double(body2),'rows')];
            
            
            if parser.Results.uglywing == 1
                TW=wingsUgly;
            else
                TW=wingsUgly(wings_all_noBod==1,:);
            end
            TW=unique([TW],'rows');
            obj.framehull.Twowings = {TW};
            obj.framehull.body4plot = {body4plot};
        end
        
        function prerp4save_hullRec(obj,body,wing,realC,frvec,body4plot)
            obj.framehull.body = body;
            obj.framehull.body4plot = body4plot;
            obj.framehull.wing = wing;
            obj.real_coord = realC;
            obj.framehull.frames = frvec;
        end
        
        function plothullLab(obj,fr,varargin)
            parser = inputParser;
            addParameter(parser,'flipcam',1); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            % plot hull in lab axes
            
            % transform indices to lab axis
            frm = (fr == obj.framehull.frames);
            bodyhull = obj.framehull.body{frm};
            winghull = obj.framehull.wing{frm};
            realC = obj.real_coord{frm}';
            [hullRealbody rotmat] = Functions.hullRec2lab(bodyhull,obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            hullRealwing = Functions.hullRec2lab(winghull,obj.cameras.all.Rotation_Matrix,obj.cameras.all.RotMat_vol,realC);
            
            bodycm = mean(hullRealbody);
            
            % plot wing and body
            plot3(hullRealwing(:,1),hullRealwing(:,2),hullRealwing(:,3),'.r');hold on
            plot3(hullRealbody(:,1),hullRealbody(:,2),hullRealbody(:,3),'.g');
            
            for k= 1:1:length(obj.cameras.cams_arr)
                dircam = obj.cameras.all.Rotation_Matrix*(parser.Results.flipcam * obj.cameras.cams_arr(k).camera_dir);
                quiver3(bodycm(1),bodycm(2),bodycm(3),dircam(1),dircam(2),dircam(3),0.005,'color','k','linewidth',2);hold on
            end
            
            xlabel('mm');ylabel('mm');zlabel('mm');
            ttl = sprintf('frame %d',fr);
            title(ttl);
            axis equal
            grid on
            box on
            legend('wing','body','cameras')
        end
    end
end




