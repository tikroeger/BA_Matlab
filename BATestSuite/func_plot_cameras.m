

%--------------------------------------------------------------------------
% File     : func_plot_cameras(Matlab M-Function)
% Author   : Till Kroeger
% Created  : 26.09.2014 (in Matlab 8.3.0.532, R2014a)
% Usage: 
%        func_plot_cameras(K, R, t, w,h, rgbface, camerascaling, ppoint, lw, img, texturealpha);
%
% Description :
%       3D plot of a camera with internal calibration K, orientation R, world position t, 
%       3x1 RGB color 'rgbface', isotropic scaling 'camerascaling'.
%       ppoint=0: in K principal point is specified as width/2,height/2, ppoint=1, p.point is set to (0,0)
%       lw (optional) : linewidth for drawing, img (optional): draw image as texture in 3D
%       texturealpha (optional): transparency with plotted image
%--------------------------------------------------------------------------


function plothandle = func_plot_cameras(K,R,t, w,h, rgbface, camerascaling, ppoint, lw)

if (nargin < 8)
    ppoint = 1;  % set pp to 0,0 else set to w/2,h/2
    lw = 1;
else
    ppoint = 1;  % set pp to 0,0 else set to w/2,h/2
    lw = 1;
end


 
    % plot new cameras
        if (ppoint==1)
            x = double([[-w/2,-h/2,1]', [-w/2,h/2,1]', [w/2,h/2,1]', [w/2,-h/2,1]']);% 2d homog. point in image plane
        else
            x = double([[0,0,1]', [0,h,1]', [w,h,1]', [w,0,1]']);% 2d homog. point in image plane
        end

                            
        camcent = -inv(R)*t';  % correct camera center in world coordinates, -R't = X, because x = RX+t -> R'x - R't = X and R'x is 0,0,0 (camera center)
        
        X = inv(K) * x;    % inverse projection
        X = X ./ repmat(sqrt(sum(X.^2,1)),3,1) * camerascaling;  % camera scaling
        X = inv(R) * (X - repmat(t',1,4));   % apply rotation and translation -> world coordinates
        
        if ~ishold 
            hold on;
            washold = 0; % remember hold state
        else
            washold = 1;
        end
        
        
        plothandle(1) = scatter3(camcent(1), camcent(2), camcent(3), 25, rgbface, 'filled');
        


        for k = 1:4 
            plothandle(k+1) = plot3([camcent(1) X(1,k)],[camcent(2) X(2,k)],[camcent(3) X(3,k)], 'Color', rgbface, 'linewidth',lw);

        end
            
        plothandle(6) = plot3([X(1,:) X(1,1)],[X(2,:) X(2,1)],[X(3,:) X(3,1)],'Color', rgbface, 'linewidth', lw);
        plothandle(7) = surf(reshape(X(1,[1 2 4 3]),2,2),reshape(X(2,[1 2 4 3]),2,2),reshape(X(3,[1 2 4 3]),2,2),'FaceAlpha',.4,'FaceColor',rgbface,'EdgeColor',rgbface,'EdgeAlpha',.5);

        
        if (washold == 0)
            hold off;
        end
    
end
    
