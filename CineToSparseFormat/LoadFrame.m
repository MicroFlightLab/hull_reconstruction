function full_im=LoadFrame(mov_path,frame_ind)
    % im=LoadFrame('C:\Users\noaml\Documents\sparse_movies\mov1_cam3.mat',100);
    % figure;imshow(im)
    mf=matfile(mov_path);
    frame=mf.frames(frame_ind,1);
    meta_data=mf.metaData;

    full_im=zeros(size(meta_data.bg),'like',meta_data.bg);
    lin_inds=sub2ind(meta_data.frameSize,frame.indIm(:,1),frame.indIm(:,2));
    full_im(lin_inds)=meta_data.bg(lin_inds)-frame.indIm(:,3);
end