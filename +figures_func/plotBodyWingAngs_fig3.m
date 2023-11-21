function plotBodyWingAngs_fig3(hull,figurepath,figsdir,stenfr,ploters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure;
ploters.BodyAng(hull,'marker','k*')
flname = ['\body angles'];
saveas(gcf,[figurepath,figsdir,flname])

figure;
ploters.wingAng(hull,'rightwing','marker','*-r','time',1,'stenfr',stenfr);hold on
ploters.wingAng(hull,'leftwing','marker','*-b','time',1,'stenfr',stenfr);hold on
flname = ['\wing angles'];
print(gcf,[figurepath,figsdir,'\svg\',flname],'-dsvg')
saveas(gcf,[figurepath,figsdir,'\fig\',flname,'.fig'])

end

