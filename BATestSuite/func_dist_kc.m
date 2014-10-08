%--------------------------------------------------------------------------
% File     : func_dist_kc(Matlab M-Function)
% Author   : Till Kroeger
% Created  : 26.09.2014 (in Matlab 8.3.0.532, R2014a)
% Usage: 
%         ptnorm_dist = func_dist_kc(ptnorm,kc)
%
% Description :
%        Distorts 2D points given in normalized image coordinates.
%        ptnorm: 2D points in normalized image coordinates
%        kc: radial distortion (1x1,2x1, 4x1, 5x1), (Brown, Decentering Distortion of lenses, 1966)
% 
%--------------------------------------------------------------------------


function ptnorm = func_dist_kc(ptnorm,kc)

if (length(kc)>0)
        tmp = ptnorm;
        if (size(tmp,1)==1)  
            tmp = tmp';
        end

        r2 = sum(tmp.^2,1);        

        rc = 1 + kc(1)*r2; 
        
        if (length(kc)>1)
            rc = rc + kc(2)*r2.*r2;
            
            if (length(kc)==5)
                rc = rc + kc(5)*r2.*r2.*r2;
            end
            
            if (length(kc)>=3)
                dx = [(2*kc(3)*tmp(1)*tmp(2) + kc(4)*(r2 + 2*tmp(1).^2)); ...  
                      (2*kc(4)*tmp(1)*tmp(2) + kc(3)*(r2 + 2*tmp(2).^2))]; % tangential distortion
            end
        end
        
        if (length(kc)>=3)
            ptnorm = bsxfun(@times,tmp,rc) + dx;
        else
            ptnorm = bsxfun(@times,tmp,rc);
        end
end
        
        
