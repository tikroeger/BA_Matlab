
function [err1_fin err2_fin vplines_proj] = func_reproject_VP_BAtest(ret, camid, viewid)

    if (size(ret.cams(camid).views_orien,2)>1) % if all views of camera share the same orientation
        R = func_rodrigues(ret.cams(camid).views_orien(:,viewid));
    else
        R = func_rodrigues(ret.cams(camid).views_orien(:,1));        
    end

    if (size(ret.cams(camid).views_trans,2)>1)  % if all views of camera share the same translation
        t = ret.cams(camid).views_trans(:,viewid);
    else
        t = ret.cams(camid).views_trans(:,1);        
    end

    err= 0;
    
    
    vp_cell = ret.cams(camid).vplines(viewid).l; 
    
    if (size(vp_cell,1)>0)
        
        err1_fin = cell(1,size(vp_cell,1));
        err2_fin = cell(1,size(vp_cell,1));
        vplines_proj = cell(1,size(vp_cell,1));
        
        for l = 1:size(vp_cell,1)
            vpid = vp_cell{l,1};
            vp = ret.VPconstraint.vp(:,vpid);
            vp = R*vp; % convert into camera centric orientation
            
            vplines = vp_cell{l,2};
            vplines_mid = [(vplines(:,1:2) + vplines(:,3:4))/2];
            nolines = size(vplines,1);
            

% 0. Apply rotation to 3D vanishing direction
% 1/2. bring 2D line endpoints to normalized image coordinates, and undistort points
% [3. compute line midpoint in normalized image coordinates]
% 3. compute original 2D line midpoint, project to normalized image coordinates, undistort
% 4. use line midpoint in normalized image coordinates to compute greater circle (VP plane through origin)
% Var1: 5). compute projection error of line endpoints (in n. i. c) to greater circle, scale proj. error by factor of focal length (TEST 1)

% Var2: 5). project line endpoints (in n. i. c.) to greater circle
% 6) transform projected points on greater circle back to 2D image coordinates, distort!
% 7) Transformed endpoints and original 2D endpoints should overlap. Project transformed endpoints on original 2D image line.
% 8) Projection error vector of both endpoints is residual(TEST 2)


            % 1./2. 
            ptnorm = bsxfun(@rdivide, bsxfun(@minus, vplines(:,1:2), ret.cams(camid).cc), ret.cams(camid).fc);
            ptnorm_undist_a = func_undist_kc(ptnorm',ret.cams(camid).kc)';
            ptnorm_undist_a = [ptnorm_undist_a ones(size(ptnorm_undist_a,1),1)];
            ptnorm_undist_a = bsxfun(@rdivide, ptnorm_undist_a, sqrt(sum(ptnorm_undist_a.^ 2,2)));
            
            ptnorm = bsxfun(@rdivide, bsxfun(@minus, vplines(:,3:4), ret.cams(camid).cc), ret.cams(camid).fc);
            ptnorm_undist_b = func_undist_kc(ptnorm',ret.cams(camid).kc)';
            ptnorm_undist_b = [ptnorm_undist_b ones(size(ptnorm_undist_b,1),1)];
            ptnorm_undist_b = bsxfun(@rdivide, ptnorm_undist_b, sqrt(sum(ptnorm_undist_b .^ 2,2)));
            
            % 3.
            ptnorm = bsxfun(@rdivide, bsxfun(@minus, (vplines(:,1:2) + vplines(:,3:4))/2, ret.cams(camid).cc), ret.cams(camid).fc);
            meanpt3D = func_undist_kc(ptnorm',ret.cams(camid).kc)';
            meanpt3D = [meanpt3D ones(size(meanpt3D,1),1)];
            meanpt3D = bsxfun(@rdivide, meanpt3D, sqrt(sum(meanpt3D.^ 2,2)));
            
            %meanpt3D = (ptnorm_undist_a +ptnorm_undist_b )/2;
            %meanpt3D = bsxfun(@rdivide, meanpt3D, sqrt(sum(meanpt3D .^ 2,2))); % renormalize
            
            % 4. 
            planevec = zeros(3,nolines);
            for k = 1:nolines
                planevec(:,k) = cross(vp, meanpt3D(k,:)');
                planevec(:,k) = planevec(:,k) ./ norm(planevec(:,k));
            end
            
            % 5. 
            err_var1 = zeros(nolines,2); 
            ptnorm_proj_a = zeros(size(ptnorm_undist_a));
            ptnorm_proj_b = zeros(size(ptnorm_undist_b));
            for k = 1:nolines 

%figure
%quiver3(0,0,0,planevec(1,k),planevec(2,k),planevec(3,k), 'linewidth',2); 
%hold on

%quiver3(0,0,0,ptnorm_undist_a(k,1),ptnorm_undist_a(k,2),ptnorm_undist_a(k,3), 'r'); 
%quiver3(0,0,0,ptnorm_undist_b(k,1),ptnorm_undist_b(k,2),ptnorm_undist_b(k,3), 'r'); 
%norm(ptnorm_undist_a(k,:))

                err_a = dot(planevec(:,k), ptnorm_undist_a(k,:)) * 1;   % projection error of pt on interpretation plane
                err_b = dot(planevec(:,k), ptnorm_undist_b(k,:)) * 1;   % projection error of pt on interpretation plane
                

                ptnorm_proj_a(k,:) = ptnorm_undist_a(k,:) - (err_a*planevec(:,k))';
                ptnorm_proj_a(k,:) = ptnorm_proj_a(k,:) ./ norm(ptnorm_proj_a(k,:));
                
                ptnorm_proj_b(k,:) = ptnorm_undist_b(k,:) - (err_b*planevec(:,k))';
                ptnorm_proj_b(k,:) = ptnorm_proj_b(k,:) ./ norm(ptnorm_proj_b(k,:));

%quiver3(0,0,0,ptnorm_proj_a(k,1),ptnorm_proj_a(k,2),ptnorm_proj_a(k,3), 'b'); 
%quiver3(0,0,0,ptnorm_proj_b(k,1),ptnorm_proj_b(k,2),ptnorm_proj_b(k,3), 'b'); 
%norm(ptnorm_proj_a(k,:))
%axis equal                

                % Variant 1: compute projection error in normalized coordinates, scale up by focal length, ignores rad. distortion img. space
                err_var1(k,:) = [sqrt(sum((ptnorm_proj_a(k,:) - ptnorm_undist_a(k,:)).^2,2)) * mean(ret.cams(camid).fc) sqrt(sum((ptnorm_proj_b(k,:) - ptnorm_undist_b(k,:)).^2,2)) * mean(ret.cams(camid).fc)];
            end  
            
            % 6. 
            ptnorm_proj_a = [ptnorm_proj_a(:,1)./ptnorm_proj_a(:,3) ptnorm_proj_a(:,2)./ptnorm_proj_a(:,3)];
            ptnorm_proj_b = [ptnorm_proj_b(:,1)./ptnorm_proj_b(:,3) ptnorm_proj_b(:,2)./ptnorm_proj_b(:,3)];
            
            ptnorm_proj_a = func_dist_kc(ptnorm_proj_a',ret.cams(camid).kc)';
            ptnorm_proj_b = func_dist_kc(ptnorm_proj_b',ret.cams(camid).kc)';

            
            ptnorm_proj_a_img = bsxfun(@plus, bsxfun(@times, ptnorm_proj_a, ret.cams(camid).fc), ret.cams(camid).cc);
            ptnorm_proj_b_img = bsxfun(@plus, bsxfun(@times, ptnorm_proj_b, ret.cams(camid).fc), ret.cams(camid).cc);
            
%[(vplines(:,1:2) + vplines(:,3:4))/2] - [(ptnorm_proj_a_img + ptnorm_proj_b_img)/2]
%meanpt2D = [meanpt3D(:,1)./meanpt3D(:,3) meanpt3D(:,2)./meanpt3D(:,3)];
%meanpt2D = func_dist_kc(meanpt2D',ret.cams(camid).kc)';
%meanpt2D = bsxfun(@plus, bsxfun(@times, meanpt2D, ret.cams(camid).fc), ret.cams(camid).cc);
%[(vplines(:,1:2) + vplines(:,3:4))/2] - meanpt2D
            
            
            % 7.  
            projlines_vec_a = vplines_mid - ptnorm_proj_a_img; 
            projlines_vec_b = vplines_mid - ptnorm_proj_b_img; 
            norm_projlines_vec_a = sqrt(sum(projlines_vec_a.^2,2));
            norm_projlines_vec_b = sqrt(sum(projlines_vec_b.^2,2));
            
            
            vplines_vec_a = vplines_mid - vplines(:,1:2); 
            vplines_vec_a = vplines_vec_a ./ repmat(sqrt(sum(vplines_vec_a.^2,2)),1,2);
            vplines_vec_b = vplines_mid - vplines(:,3:4); 
            vplines_vec_b = vplines_vec_b ./ repmat(sqrt(sum(vplines_vec_b.^2,2)),1,2); 
            
            ang_a = sum(projlines_vec_a .* vplines_vec_a,2) ./ norm_projlines_vec_a;
            %acosd(ang_a) 
            ang_b = sum(projlines_vec_b .* vplines_vec_b,2) ./ norm_projlines_vec_b;
            %acosd(ang_b) 
            
            % 8.
            % projected vector as (cos(angle) * lengthof_projlines_vec_from_intersection) * unit_vector_in_vplines_vec_direction
            errvec_a = repmat(          ang_a .*      norm_projlines_vec_a,1,2)     .* vplines_vec_a;  
            errvec_b = repmat(          ang_b .*      norm_projlines_vec_b,1,2)     .* vplines_vec_b; 

            % convert to endpoint projection error vectors
            errvec_a = projlines_vec_a - errvec_a; 
            errvec_b = projlines_vec_b - errvec_b; 
             
            err_var2 = [errvec_a; errvec_b]; 
 
            err1_fin{l} = err_var1;
            err2_fin{l} = err_var2;
            
            vplines_proj{l} = [[vplines_mid ptnorm_proj_a_img]; [vplines_mid ptnorm_proj_b_img];];
        end
        
        
    else
        err_var1 = 0;
        err_var2 = 0;
        vplines_proj = [];
        
    end