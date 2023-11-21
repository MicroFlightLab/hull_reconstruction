function TwoD = Create2DCoords(realC,coefs,coords_cell,notLab,ImHight,RotMAt_vol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for k=1:1:length(coords_cell)
TwoD{k} = Functions.ThreeD2real(realC,coords_cell{k},coefs,notLab,ImHight,RotMAt_vol);
end


end

