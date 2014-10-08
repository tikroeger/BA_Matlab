function pt2d = func_reproject_BAtest(pt3d, R,t,fc,cc,kc, camcenter)


        if (camcenter==1)  % t is set in world frame 
            if (size(pt3d,2)==1)
                tmp = (R*(pt3d -t));        
            else
                ttmp = bsxfun(@minus, pt3d, t');
                tmp = arrayfun(@(x) R*ttmp(x,:)', 1:size(ttmp,1), 'UniformOutput',0);
                tmp = cat(2,tmp{:});
            end
        else  % t is set in camera reference frame
            if (size(pt3d,2)==1)
                tmp = (R*pt3d + t);        
            else
                ttmp = arrayfun(@(x) R*pt3d(x,:)', 1:size(pt3d,1), 'UniformOutput',0);
                ttmp = cat(2,ttmp{:});
                tmp = bsxfun(@plus, ttmp, t);
            end
        end

        %[R t] * [pt3d(1,:) 1]' 
        %R*pt3d(1,:)'+t
        
        
        %tmp = [(tmp(1) / tmp(3)) (tmp(2) / tmp(3))];   
        tmp = [(tmp(1,:) ./ tmp(3,:)); (tmp(2,:) ./ tmp(3,:))];
         
         
        r2 = sum(tmp.^2,1);        
        if (length(kc)==5) % 2nd, 4th, 6th order sym. rad. dist, + tangential dist.
            dx = [(2*kc(3)*tmp(1)*tmp(2) + kc(4)*(r2 + 2*tmp(1).^2)); ...  
                  (2*kc(4)*tmp(1)*tmp(2) + kc(3)*(r2 + 2*tmp(2).^2))]; % tangential distortion
            rc = 1 + kc(1)*r2 + kc(2)*r2.*r2 + kc(5)*r2.*r2.*r2;
            tmp = bsxfun(@times,tmp,rc) + dx;
        elseif (length(kc)==4) % 2nd, 4th order sym. rad. dist, + tangential dist.
            dx = [(2*kc(3)*tmp(1)*tmp(2) + kc(4)*(r2 + 2*tmp(1).^2)); ...  
                  (2*kc(4)*tmp(1)*tmp(2) + kc(3)*(r2 + 2*tmp(2).^2))]; % tangential distortion
            rc = 1 + kc(1)*r2 + kc(2)*r2.*r2;
            tmp = bsxfun(@times,tmp,rc) + dx;    
        elseif (length(kc)==2) % 2nd, 4th order sym. rad. dist
            rc = 1 + kc(1)*r2 + kc(2)*r2.*r2;            
            tmp = bsxfun(@times,tmp,rc);
        elseif (length(kc)==1) % 2nd order sym. rad. dist
            rc = 1 + kc(1)*r2;            
            tmp = bsxfun(@times,tmp,rc);
        end
        
        
        %r2 = sum(tmp.^2);        
        %rc = 1 + kc(1)*r2 + kc(2)*r2*r2;
        %tmp = rc*tmp;
        
        
        
        pt2d = [(fc(1) * tmp(1,:)); (fc(2) * tmp(2,:))];
        if (size(pt2d,2)==1)
            pt2d(1)  = pt2d(1)  + cc(1);
            pt2d(2)  = pt2d(2)  + cc(2);
        else
            if (size(cc,2)==1)
                pt2d  = bsxfun(@plus, pt2d, cc);
            else
                pt2d  = bsxfun(@plus, pt2d, cc');
            end
        end
% 
% K = eye(3);
% K(1,1) = fc(1); K(2,2) = fc(2); 
% K(1,3) = cc(1); K(2,3) = cc(2);
% 
% P = K * [R t];
% 
% pt2d  = P*[pt3d' 1]';
% 
% pt2d = [pt2d(1)/pt2d(3) pt2d(2)/pt2d(3)];






%         if (camcenter==1)  % t is set in world frame 
%             tmp = R*(pt3d -t);        
%         else
%             tmp = R*pt3d + t;        
%         end
%         
%         tmp = [(tmp(1) / tmp(3)) (tmp(2) / tmp(3))];
%         
%         
%         r2 = sum(tmp.^2);        
%         if (length(kc)==5) % 2nd, 4th, 6th order sym. rad. dist, + tangential dist.
%             dx = [(2*kc(3)*tmp(1)*tmp(2) + kc(4)*(r2 + 2*tmp(1).^2)) ...  
%                   (2*kc(4)*tmp(1)*tmp(2) + kc(3)*(r2 + 2*tmp(2).^2))]; % tangential distortion
%             rc = 1 + kc(1)*r2 + kc(2)*r2*r2 + kc(5)*r2*r2*r2;
%             tmp = rc*tmp + dx;
%         elseif (length(kc)==4) % 2nd, 4th order sym. rad. dist, + tangential dist.
%             dx = [(2*kc(3)*tmp(1)*tmp(2) + kc(4)*(r2 + 2*tmp(1).^2)) ...  
%                   (2*kc(4)*tmp(1)*tmp(2) + kc(3)*(r2 + 2*tmp(2).^2))]; % tangential distortion
%             rc = 1 + kc(1)*r2 + kc(2)*r2*r2;
%             tmp = rc*tmp + dx;              
%         elseif (length(kc)==2) % 2nd, 4th order sym. rad. dist
%             rc = 1 + kc(1)*r2 + kc(2)*r2*r2;            
%             tmp = rc*tmp;
%         elseif (length(kc)==1) % 2nd order sym. rad. dist
%             rc = 1 + kc(1)*r2;            
%             tmp = rc*tmp;            
%         end
%         
%         
%         %r2 = sum(tmp.^2);        
%         %rc = 1 + kc(1)*r2 + kc(2)*r2*r2;
%         %tmp = rc*tmp;
%         
%         
%         
%         pt2d = [(fc(1) * tmp(1)) (fc(2) * tmp(2))];
%         pt2d  = pt2d + cc;
% % 
% % K = eye(3);
% % K(1,1) = fc(1); K(2,2) = fc(2); 
% % K(1,3) = cc(1); K(2,3) = cc(2);
% % 
% % P = K * [R t];
% % 
% % pt2d  = P*[pt3d' 1]';
% % 
% % pt2d = [pt2d(1)/pt2d(3) pt2d(2)/pt2d(3)];
% 
