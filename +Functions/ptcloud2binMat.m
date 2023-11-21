function [mat3d] = ptcloud2binMat(obj_cor)

x=obj_cor(:,1);
y=obj_cor(:,2);
z=obj_cor(:,3);

mat3d=zeros(max(x),max(y),max(z));
linearInd = sub2ind([max(x),max(y),max(z)], x, y, z);
if sum(isnan(linearInd))>0;
linearInd(isnan(linearInd))=[];
end
mat3d(linearInd)=1;
end

