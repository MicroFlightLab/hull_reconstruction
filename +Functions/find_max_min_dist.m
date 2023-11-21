 function [max_old,ind_max_op,min_old,ind_min_op] = find_max_min_dist(vec3d,comp_points)
            max_old=0;
            min_old=1000;
            for k=1:1:size(vec3d,1)
                dist=sqrt(sum((vec3d(k,:)-comp_points).^2,2));
                [max_dist ~]=max(dist);
                [min_dist ~]=min(dist);
                
                if max_dist>max_old
                    max_old=max_dist;
                    ind_max_op=k;
                end
                if min_dist<min_old
                    min_old=min_dist;
                    ind_min_op=k;
                end
            end
            
        end

