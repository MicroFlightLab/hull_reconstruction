function [vecspn,crd] = plotWings_sec_fig2(hull,hull3d,fr,wingsz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ind2mm = 50e-6*1000;

LETEname = {'LE','TE'};
wingname= {'leftwing','rightwing'};
frm = find(fr == hull.frames);

for kwing = 1:1:2
    %     figure;
    for kLETE = 1:1:2
        if kLETE == 1
            cm4plt = 0;
        end
        
        [cm4plt,vecspn,vecCM4plt_wing] =  hull.calcSecChord_frame(LETEname{kLETE},wingname{kwing},hull3d,fr,...
            'tipprop','projtip','plot_sec',1,'vecCM4plt',0,'wingsize',wingsz);
        if kLETE == 2
            vecspn = vecspn/norm(vecspn);
            pln =  (hull.rotmat_EWtoL*hull.(wingname{kwing}).vectors.LETEPlane(frm,:)')';
            crd = cross(vecspn,pln);

            bod = (ind2mm*hull.rotmat_EWtoL*(double(hull3d.body.hull{frm})'))';
            meanbod = mean(bod);
            
            LE = (ind2mm*hull.rotmat_EWtoL*(double(hull3d.(wingname{kwing}).hull.LE{frm})'))' - meanbod;
            TE = (ind2mm*hull.rotmat_EWtoL*double(hull3d.(wingname{kwing}).hull.TE{frm})')' - meanbod;
            bound = [LE;TE];
            meanB = mean(bound);
            ampspn = 10;
            ampcrd = 12;
            sqr = [vecCM4plt_wing(4,:) + vecspn*ampspn + crd*ampcrd;vecCM4plt_wing(4,:) + vecspn*ampspn - crd*1.1*ampcrd...
                ;vecCM4plt_wing(4,:) - vecspn*ampspn*3 - crd*1.1*ampcrd;vecCM4plt_wing(4,:) - vecspn*ampspn*3 + crd*ampcrd];

            X = sqr(:,1);
            Y = sqr(:,2);
            Z = sqr(:,3);
            
%             plot3(X,Y,Z,'.');
%             p =fill3(X,Y,Z,'r');
            
        end
    end
    p.FaceAlpha = 0.2;
    axis equal
    grid on
    xlabel('mm');ylabel('mm');zlabel('mm');
    box on
    view(3)
end
end

