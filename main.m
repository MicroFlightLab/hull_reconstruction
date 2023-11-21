

clear
close all
clc
warning off





path='H:\My Drive\data_for_hull_git\'; % location of the movie folder
nameOFeasyFile='wand_data1_19_05_2022_skip5_easyWandData'; % calibration file ( located inside path)

mov = 17
loadSeg = 0 % 1 - load segmentation
loadhull_op =0 % 1 load hull reconstrucion 
loadhull = 0 % 1 - load analyzed hull
waitbar = 1; % show wait bar (for hull reconstrucaion and anlysis)

flythin = 1; % change the amount of pixels to remove around the 2d image ( 0 for synthetic data)
framesSeg = [100 800]%; % if 0 use all sparse frames, [a b] - from frame a to b
frmshullRec = [0]; % if 0 use segmentation frames [a:b]
frms2runAna = [0]; % if 0 use hull reconstruction frames
hullRec4debug = 0; % 1 - debug hull reconstruction(not parallel)
hullAnaDebug = 0; % 1 - debug hull analysis(not parallel)
segdebug = 0;  % 1 - debug segmentation(not parallel)
camvec = [1,2,3,4] % order of camera (according to easywand file)
zcam =  1 % camera pointing to Z world direction


loaders = loaders_class(path,mov,nameOFeasyFile,'hullfile','//hull_op//');
easy = loaders.easywand();

hullRec = loaders.loadhullRec(loadhull_op);
seg = loaders.loadSegfile(loadSeg);
sp  = loaders.loadsparse();
[hull,hull3d] = loaders.loadhull(loadhull);
myCluster = parcluster('local')
delete(myCluster.Jobs)


%% parallel segmentation of wings and body from 2D images

if loadSeg == 0
    
    codepart = 'running segmentation';
    fprintf(codepart)

    % start and end frames. if 0 - run all. user input [start end]
    seg = seg_class(sp,'st_en_fr',framesSeg,'bodyTH',72);

    for in = 1:1:length(seg.camvec)
        if segdebug == 1
            [body_op,wing1_op,wing2_op,image_op,frames] = parallelFun.par_seg(sp,seg.camvec(in),seg);
            body{in} = body_op{1};
            wing1{in} = wing1_op{1};
            wing2{in} = wing2_op{1};
            image{in} = image_op{1};
            %                 st_enfr = frames(1,:);
        else
            f(in) = parfeval(@parallelFun.par_seg, 5,sp,seg.camvec(in),seg);
        end
    end
    if segdebug == 0
        [idx, result] = fetchNext(f);
        [body,wing1,wing2,image,frames] = fetchOutputs(f);
    end
    % save output in class and save class
    seg.prep4save(body,wing1,wing2,image,frames);
    save([loaders.segpath,loaders.segfile],'seg');


end

%% parallel hull reconstruction


if loadhull_op == 0
    if frmshullRec == 0
        frmshullRec = [seg.st_enfr(1):seg.st_enfr(2)];
    end
    codepart = 'running hull reconstruction';
    fprintf(codepart)
    hullRec = hullrec_class(easy,'ofst',2,'VxlSize4search',10e-5,...
        'camvec',camvec,'ZaxCam',zcam,'VxlSize',50e-6,'VolLen',10e-3);

    field_string = {'all','body','wing1','wing2'};
    hullRec.genNames4rec(field_string);
    indpar = 1;
    [ind0_frame,framevolume] = hullRec.createVol(hullRec.voxelSize4search);



    for idx = 1:length(frmshullRec)
        try
            %
            % run parallrl computation for hull reconstruction
            fr = frmshullRec(idx);
            hullRec.Seg2hullRIm(seg,fr);
            %     % build an initial volume from "all" * 4 images
            if sum(cellfun(@isempty,hullRec.im4hull.sprs.body)>0)
                continue
            end
            hullRec.FindSeed('all');
            hullRec.hull_params();

            if hullRec4debug == 0
                parhull(indpar) = parfeval(@parallelFun.par_hullRec,4,hullRec,hullRec.parts2run,ind0_frame);
            else
                [body,Twowings,realC,body4plot] =  parallelFun.par_hullRec(hullRec,hullRec.parts2run,ind0_frame);
            end
            if exist('parhull') && isempty(parhull(indpar).Error) == 0
                continue
            end
            frvec(indpar) = fr;
            indpar = indpar + 1;
        catch
            continue
        end
    end

    parallelFun.fetchNext_waitbar(parhull,frvec,loaders.mov,codepart,'waitbar',1); % run parallel function and show progress bar
    [body,wing,realC,body4plot] = fetchOutputs(parhull); % get outputs from parallel run

    hullRec.prerp4save_hullRec(body,wing,realC,frvec,body4plot); % arange in class

    save([loaders.hullpath,loaders.hullRecfile],'hullRec');
else
    load([loaders.hullpath,loaders.hullRecfile],'hullRec');
end

clear parhull
clear frvec
%% parallel cleaning wing, calculating X body and span
indpar = 1;

if loadhull == 0
    codepart = 'calculating vectors and boundaries';
    fprintf(codepart)
    if frms2runAna == 0
        frms2runAna = [hullRec.framehull.frames(1),hullRec.framehull.frames(length(hullRec.framehull.frames))];
    end

    for fr = frms2runAna(1) : frms2runAna(2)
        try
            % initilize hull analysis class
            hull = hullAna_class(hullRec,'bound_dilate',0,'coneang',15);%,'percWing_distfromBod',[0.35 0.45]);
            % input current frame wing and body from hull reconstruction
            frm = (fr == hullRec.framehull.frames);
            hull.body.hull3d = hullRec.framehull.body{frm};
            bodyfrm{indpar} = hull.body.hull3d;
            winghull = hullRec.framehull.wing{frm};
            realCfrm{indpar} = hullRec.real_coord{frm}';

            if hullAnaDebug == 1
                parallelFun.BoundaryAndVectors(hull,realCfrm{indpar},winghull);
            else
                % run parallel computation
                parclean(indpar) = parfeval(@parallelFun.BoundaryAndVectors,1,hull,realCfrm{indpar},winghull);
            end
            if length(parclean)>1 && ~isempty(parclean(indpar).Error)
                continue
            end
            frvecclean(indpar) = fr;
            indpar = indpar + 1;
        catch
            continue
        end
    end

    parallelFun.fetchNext_waitbar(parclean,frvecclean,loaders.mov,codepart,'waitbar',waitbar); % run parallel function and show progress bar
    [op] = fetchOutputs(parclean,'UniformOutput',false); % get outputs from parallel run
    hull3d = hull.prep4save_spnXTip(op,bodyfrm,realCfrm,frvecclean); % arange in class
    hull3d.body.body4plot = hullRec.framehull.body4plot(1:length(hull.frames));

    hull.body.coords.CM = cell2mat(cellfun(@(x) mean(x,1),hull3d.body.hull' ,'UniformOutput',false));
    save([loaders.hullpath,loaders.hullfile],'hull');
    save([loaders.hullpath,loaders.hull3dname],'hull3d');

end
clear bodyfrm
clear realCfrm
clear parclean
clear frvecclean


%%
%--- continue analysis to calculate Y,Z, define right and
%left wing and calculate the wing vectors from the boundary
hull.video.cineFrame = sp{1}.metaData.startFrame;
hull.video.timeframe = 1/sp{1}.metaData.frameRate * (sp{1}.metaData.startFrame + hull.frames)*1000;
hull.BodyAxesYZV2; % calculate Y and Z axis
hull3d = hull.RightLeftWing(hull3d); % define right and left wing
for fr = hull.frames
    frm = (fr == hull.frames);
    if sum(isnan(hull3d.body.hull{frm})) == 0
        CM = mean(Functions.hullRec2lab(hull3d.body.hull{frm},hull.cameras.all.Rotation_Matrix,hull.cameras.all.RotMat_vol,hull.real_coord{frm}));
    end
    hull3d =  hull.calcNormChord_bound(hull3d,frm); % calculate normal and chord from boundary
end

hull.bodyAngles();
hull.wingAngles('spanProp','span','chordProp','chord_bound');


%% calculate psi angle for TE and LE for multiple sections (defined by numberof_sec)
wingname= {'leftwing','rightwing'};
hull.wingana.numofsec = 5;
LETEname = {'LE','TE'};
for kwing = 1:1:2
    %     figure;
    for kLETE = 1:1:2
        hull.calcSecChord(LETEname{kLETE},wingname{kwing},hull3d,'tipprop','projtip');
        hull.calcSecPsi(LETEname{kLETE},wingname{kwing})
    end
    %             hull.plotPsi_sec(wingname{kwing})
    %             hold on;plot(hull.(wingname{kwing}).angles.psi,'-k','linewidth',3);
end
Shull = struct(hull);
save([loaders.hullpath,loaders.hullfile],'hull');
save([loaders.hullpath,loaders.hull3dname],'hull3d');
save([loaders.hullpath,loaders.struct_hullfile],'Shull');





%%
