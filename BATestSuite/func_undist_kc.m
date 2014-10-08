%--------------------------------------------------------------------------
% File     : func_undist_kc(Matlab M-Function)
% Author   : Till Kroeger
% Created  : 26.09.2014 (in Matlab 8.3.0.532, R2014a)
% Usage: 
%         ptnorm_undist = func_undist_kc(ptnorm , kc)
%
% Description :
%        Undistorts 2D points given in normalized image coordinates.
%        Returns an approximative, iterative solution.
%        ptnorm: 2D points in normalized image coordinates
%        kc: radial distortion (1x1,2x1, 4x1, 5x1), (Brown, Decentering Distortion of lenses, 1966)
% 
%--------------------------------------------------------------------------



function ptnorm = func_undist_kc(ptnorm , kc)

if (size(ptnorm,1)==1)  
    ptnorm = ptnorm';
end

if length(kc) > 0  
    x = ptnorm; 

    for kit=1:20
        delta_x = 0;
        r_2 = sum(x.^2);
        k_radial = 1 + kc(1) * r_2;
        
        if (length(kc)>1)
            k_radial = k_radial + kc(2) * r_2.^2;

            if (length(kc)==5)
                k_radial = k_radial  + kc(5) * r_2.^3;
            end            

            if (length(kc)>=3)
                delta_x = [2*kc(3)*x(1,:).*x(2,:) +      kc(4)*(r_2 + 2*x(1,:).^2); ...
                           kc(3) * (r_2 + 2*x(2,:).^2)+  2*kc(4)*x(1,:).*x(2,:)];
            end

        end

        x = (ptnorm - delta_x)./(ones(2,1)*k_radial); 

    end
    ptnorm = x;

end 