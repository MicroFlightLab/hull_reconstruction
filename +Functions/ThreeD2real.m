function [points,real_part] = ThreeD2real(real_c,ind_part,easyWandData,notLab,ImHight,RotMAt_vol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if notLab == 1
    ind_part = round(ind_part);
    real_part=[real_c{1}(ind_part(:,1));real_c{2}(ind_part(:,2));real_c{3}(ind_part(:,3))]';
else
    real_part = ind_part;
end
real_part = (RotMAt_vol'*(real_part'))';
invDLT = Functions.dlt_inverse(easyWandData, real_part);
points = round(invDLT);


points(:,2)=ImHight-points(:,2);

end

