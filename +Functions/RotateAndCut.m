function [new_tail_points,new_tail_points_real] = RotateAndCut(wing,DirVec,PercToLe)
%% align the tail's primary vector with the x-axis
%     body_center=mean(body);
    old_tail_points=wing;
    old_tail_mean=mean(old_tail_points);
    old_tail_points_cm=old_tail_points-old_tail_mean; %center the data
    roti = vrrotvec2mat(vrrotvec(DirVec,[1,0,0]));

    old_tail_cen_main=old_tail_points_cm*roti';
   

   
    %% determine which side is the tip of the tail and cut from the other side
    min_cor = min(old_tail_cen_main(:,1));
    max_cor = max(old_tail_cen_main(:,1));
    length_tail_x=max_cor(1)-min_cor(1);
    
            good_part=old_tail_cen_main(old_tail_cen_main(:,1)>...
            (min_cor(1)+PercToLe*length_tail_x),:);
    % reverse transformation
    parta=(roti'*good_part')';
    new_tail_points=round(parta+old_tail_mean);
    new_tail_points_real = parta+old_tail_mean;
end

