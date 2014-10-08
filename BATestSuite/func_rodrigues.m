%--------------------------------------------------------------------------
% File     : func_rodrigues (Matlab M-function)
% Author   : Andras Bodis-Szomoru
% Created  : 13.10.2011 (in Matlab 7.6.0.324, R2008a)
% Usage: 
%       R = func_rodrigues(a);
%       R = func_rodrigues(a,NumTol);
%       a = func_rodrigues(R);
%       a = func_rodrigues(R,NumTol);
%
% Description : 
%        Compute the angle-axis representation of a 3D rotation given by
%        its 3x3 rotation matrix R and vice-versa. a is a 3x1 vector, thus
%        angle-axis representation is a minimal parametrisation of the 3
%        degrees-of-freedom rotation. The direction of a is parallel to the
%        rotation axis while its norm gives the rotation angle in radians.
%        NumTol is a numerical tolerance value for internal tests, e.g.
%        rank computation, orthonormality checking for R etc. The default
%        value is 1e-10.

%--------------------------------------------------------------------------

function y = func_rodrigues(x,NumTol)

if nargin==1; NumTol=1e-10; elseif nargin~=1; error('Invalid number of input arguments'); end
[m,n] = size(x);
if m~=3; error('x must be a 3x1 vector or a 3x3 orthonotmal matrix'); end

I = eye(3);
%--------------------------------------------------------------------------
%           Angle-axis (Rodrigues-) vector to rotation matrix
%--------------------------------------------------------------------------
if n==1 
    a = x;
    alpha = norm(a); % the angle in radians
    C = crossMatrix(a);
    D = C*C; % diadic product x'*x
    if alpha==0
        y = eye(3);
    else
        c = sin(alpha)/alpha;
        d = (1-cos(alpha))/alpha^2;
        y = I + c*C + d*D;
    end
%--------------------------------------------------------------------------
%           Rotation matrix to angle-axis (Rodrigues-) vector
%--------------------------------------------------------------------------
elseif n==3
    R = x;
    if norm(R'*R-I,'fro')>NumTol; error('R is not an orthonormal matrix'); end
    if abs(det(R)-1)>NumTol; error('R is not a valid rotation matrix as det(R) is not 1'); end
    
    [U,S,V]=svd(R-I);
    if S(3,3)>NumTol; error('No rotation axis found'); end
    a1 = V(:,end); % the 1st unit-length solution of the equation (R-I)x=0 or Rx=x, with the smallest singular value
    q  = [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
    sinalpha = 0.5*a1'*q;
    cosalpha = 0.5*(R(1,1)+R(2,2)+R(3,3)-1);
    alpha = atan2(sinalpha,cosalpha);    % -pi< alpha <= pi
    % if (alpha,a1) is a solution then (-alpha,-a1) is another one
    % since a rotation around axis a1 with angle alpha is equivalent with a
    % rotation around axis -a1 with angle -alpha, however, this fact does
    % not affect the product a = alpha*a1
    % Note: if alpha<0,then norm(a)=-alpha and (-a1) can  be extracted from a
    a = alpha*a1;
    y = a;
else
    error('x must be a 3x1 vector or a 3x3 orthonormal matrix');
end
    

function S = crossMatrix(v)
vx = v(1);
vy = v(2);
vz = v(3);
S = [0 -vz vy;vz 0 -vx;-vy vx 0];