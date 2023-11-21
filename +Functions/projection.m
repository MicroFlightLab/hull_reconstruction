
function [VecOnPln,notNorm]=projection(spn,strkpln)
% given an plane equation ax+by+cz=d, project points xyz onto the plane
% return the coordinates of the new projected points

projonN = dot(spn',strkpln')'; % project spn on normal to stroke plane
VecOnPln = (spn -projonN.*strkpln);
notNorm = VecOnPln;

VecOnPln = VecOnPln./vecnorm(VecOnPln,2,2);

end
