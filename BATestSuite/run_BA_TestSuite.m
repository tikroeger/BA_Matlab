%--------------------------------------------------------------------------
% File     : run_BA_TestSuite(Matlab M-Script)
% Author   : Till Kroeger
% Created  : 26.09.2014 (in Matlab 8.3.0.532, R2014a)
% Usage: 
%        varargout = run_BA_TestSuite(problem, filename)
%
% Description :
%       Edxample and test code for Bundle Adjustment with BA_Matlab.
%       1) Only points and camera residuals
%       2) Only points, cameras, and points assigned planes (part of the optim.)
%       3) Only points, cameras, and cameras assigned to planes (part of the optim.)
%       4) Only points, cameras and vanishing points
%       5) Only points, cameras, and mutuals camera visibility constraints
%
%--------------------------------------------------------------------------



addpath(['../build'])

% setup randomization options
RandSeed = 22; % make sure all runs use the same 'random' parameter changes
rndopt.randsigma_pt3d = 5; %,  3d uncertainty
rndopt.randsigma_pt2d = [0 0 0]; % 2d measurement uncertainty, (2D observations will not be optimized)
rndopt.randsigma_cam = [6 6 6]; % positional uncertainty
rndopt.randsigma_vpang=6; % Sigma for VP angular uncertainty
rndopt.fc = [100 100 100];
rndopt.cc = [50 50 50]; 
rndopt.kc = [.5 .5 .5]; 
rndopt.maxrandangle_cam = [15 15 15];
nocams = [8 8 6];

%% 1) Only points and camera residuals
% get ground truth
[ret_gt, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, [], [0 0 0 0],nocams);
func_plot_dataset(ret_gt, [1 0 1 0 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, rndopt, [0 0 0 0], nocams);
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 0 1 0 0], wh, ccent);
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

ret_test.BAopt.maxiter = uint64(30);
%func_readwrite_BA_problem(ret_test, '/tmp/BAproblem_test.txt');  % save org. problem to file
%system('../build/BAdjustBin   /tmp/BAproblem_test.txt  /tmp/out.txt');  % call BA in command line, non-mex version
%result = func_readwrite_BA_problem([], '/tmp/out.txt');    % read result
result = BAdjustMex(ret_test);

figure;
func_plot_dataset(result, [1 0 1 0 0], wh, ccent);
set(gcf, 'Position', [1920 1 1280 1000])    

% Test: Render views:
% for k = 1:ret_test.cams(1).noviews
%     figure; func_render_view(result, wh, camptidx, 1, k); title(['cam1 view ' num2str(k)]);
% end
% for k = 1:ret_test.cams(2).noviews
%     figure; func_render_view(result, wh, camptidx, 2, k); title(['cam2 view ' num2str(k)]);
% end
% for k = 1:ret_test.cams(3).noviews
%     figure; func_render_view(result, wh, camptidx, 3, k); title(['cam3 view ' num2str(k)]);
% end

%% 2) Only points, cameras, and point planes
% get ground truth
[ret_gt, wh, ccent, maxaxis] = func_create_dataset(RandSeed, [], [1 0 0 0],nocams);
func_plot_dataset(ret_gt, [1 1 1 0 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])


% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [1 0 0 0],nocams);
ret_test_bak = ret_test;
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 1 1 0 0], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

ret_test.BAopt.maxiter = uint64(30);
result = BAdjustMex(ret_test);
figure;
func_plot_dataset(result, [1 1 1 0 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    

%% 3) Only points, cameras, and camera planes
% get ground truth
[ret_gt, wh, ccent, maxaxis] = func_create_dataset(RandSeed, [], [0 1 0 0],nocams);
func_plot_dataset(ret_gt, [1 0 1 1 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [0 1 0 0],nocams);
ret_test_bak = ret_test;
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 0 1 1 0], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

ret_test.BAopt.maxiter = uint64(30);
ret_test.BAopt.maxiter = uint64(30);
ret_test.cams(1).fixInternals = uint8([1 1 1]);
ret_test.cams(2).fixInternals = uint8([1 1 1]);
ret_test.cams(3).fixInternals = uint8([1 1 1]);
result = BAdjustMex(ret_test);
clf;
figure
func_plot_dataset(result, [1 0 1 1 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    

%% 4) Only points, cameras and vanishing points
% get ground truth
[ret_gt, wh, ccent, maxaxis] = func_create_dataset(RandSeed, [], [0 0 0 1],nocams);
ret_gt.cams(1).fixInternals(:)=1;
ret_gt.cams(2).fixInternals(:)=1;
ret_gt.cams(3).fixInternals(:)=1;

rndopt_t = rndopt;
rndopt_t.fc = [0 0 0];  % was: 200
rndopt_t.cc = [0 0 0];  % was: 50
rndopt_t.kc = [0 0 0];  % was: 1
[ret_test, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, rndopt_t, [0 0 0 1],nocams);
ret_test.cams(1).fixInternals(:)=1;
ret_test.cams(2).fixInternals(:)=1;
ret_test.cams(3).fixInternals(:)=1;


for k = 1:ret_test.cams(1).noviews
    figure; func_render_view(ret_test, wh, camptidx, 1, k); title(['cam1 view ' num2str(k)]);
    set(gcf, 'Position', [2221         178        1049         692])
end


result = BAdjustMex(ret_test);
sum((ret_gt.VPconstraint.vp - ret_test.VPconstraint.vp).^2)
sum((ret_gt.VPconstraint.vp - result.VPconstraint.vp).^2)

for k = 1:result.cams(1).noviews
    figure; func_render_view(result, wh, camptidx, 1, k); title(['cam1 view ' num2str(k)]);
    set(gcf, 'Position', [2221         178        1049         692])
end
figure
func_plot_dataset(result, [1 0 1 0 0], wh, ccent);
set(gcf, 'Position', [1920 1 1280 1000])    

%func_readwrite_BA_problem(ret_test, '/tmp/BAproblem_test.txt');  % save org. problem to file
%system('../build/ceresBAbin   /tmp/BAproblem_test.txt  /tmp/out.txt');  % call BA in command line, non-mex version
%result = func_readwrite_BA_problem([], '/tmp/out.txt');    % read result

%% 5) Only points, cameras, and mutuals camera visibility constraints
% get ground truth
[ret_gt, wh, ccent, maxaxis] = func_create_dataset(RandSeed, [], [0 0 1 0],nocams);
func_plot_dataset(ret_gt, [1 0 1 0 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [0 0 1 0],nocams);
ret_test_bak = ret_test;
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 0 1 0 1], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

ret_test.BAopt.maxiter = uint64(30);
result = BAdjustMex(ret_test);
figure;
func_plot_dataset(result, [1 0 1 0 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    


%% 6) Everything together
% get ground truth
[ret_gt, wh, ccent, maxaxis] = func_create_dataset(RandSeed, [], [1 1 1 1],nocams);
func_plot_dataset(ret_gt, [1 1 1 1 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [1 1 1 1],nocams);
ret_test_bak = ret_test;
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 1 1 1 1], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

ret_test.BAopt.maxiter = uint64(30);
result = BAdjustMex(ret_test);
figure;
func_plot_dataset(result, [1 1 1 1 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    


%% Points, point planes,  camera planes, mutual camera visibility, CAMPATH SMOOTHING + CAUCHY LOSS
ret_test.BAopt.PTLoss = double([3 10]);
ret_test.BAopt.PlaneLoss = double([0 0]);
ret_test.BAopt.DerivLoss = double([0 0]);

ret_test.cams(2).smootherM = double([0.42 -1.12 1.4 -1.12 0.42]');
result = BAdjustMex(ret_test);
clf;
func_plot_dataset(result, [1 1 1 1 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    

%func_readwrite_BA_problem(ret_test, '../src/BAproblem_test.txt');  % save org. problem to file



%% 3D Points (low noise and fixed), point planes, camera planes, mutual camera visibility, strong 2D measurement noise
rndopt.randsigma_pt2d = [10 10 10];
rndopt.randsigma_pt3d = 1;
% get noisy toy data, and solve wt
figure
[ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [1 1 1 0],nocams);
ret_test.points.fixPosition(:) = 1;
axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
func_plot_dataset(ret_test, [1 1 1 1 1], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)
for k = 1:8
     figure; func_render_view(ret_test, wh, camptidx, 1, k); title(['cam1 view ' num2str(k)]);
     figure; func_render_view(ret_test, wh, camptidx, 2, k); title(['cam2 view ' num2str(k)]);
     figure; func_render_view(ret_test, wh, camptidx, 3, k); title(['cam3 view ' num2str(k)]);
end


ret_test.BAopt.maxiter = uint64(30);
result = BAdjustMex(ret_test);
clf;
func_plot_dataset(result, [1 1 1 1 1], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])    

%% 3D Points (low noise and fixed), low 2D measurement noise, video sequence
rndopt.randsigma_pt3d = 6; %was: 5,  3d uncertainty
rndopt.randsigma_pt2d = [0 2 0]; %was: 0,  2d measurement uncertainty, (2D observations will not be optimized)
rndopt.randsigma_cam = [0 15 0]; %was: 10, positional uncertainty
rndopt.fc = [0 0 0];  % was: 200
rndopt.cc = [0 0 0];  % was: 50
rndopt.kc = [0 0 0];  % was: 1
rndopt.maxrandangle_cam = [0 5 0]; % was: 10, orientation uncertainty
nocams = [8 150 8];

% create and display noise-less dataset
[ret_gt, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, [], [0 0 0 0],nocams);
ret_test.points.fixPosition(:) = 1;
func_plot_dataset(ret_gt, [1 0 1 0 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

% compute smoother matrix
hi = 1; 
p = .1; 
n = ret_gt.cams(2).noviews;

delta = zeros(n-2,n); 
for i  = 1:(n-2)
    delta(i,i) =  1 / hi;
    delta(i,i+1) =  - 1 / hi - 1 / hi;
    delta(i,i+2) =  1 / hi;
end
W = zeros(n-2,n-2);  size(W)
for i  = 1:(n-2)
    W(i,i) =  (hi+hi)/3;
    if (i>1)
        W(i-1,i) =  hi/6;
        W(i,i-1) =  hi/6;
    end
end
K = delta' * inv(W) * delta;
%Weig = eye(size(K,1));
%Weig(1,1) = 10; Weig(end,end)=10; %Weig(end,end) = 10;
%LL = inv((Weig +  (1/p)*K - K))*Weig; % Weighted Smoother matrix
LL = inv((eye(n)+  (1/p)*K - K)); % Smoother matrix

% interpolate camera position
cp = ret_gt.cams(2).views_trans;
ct = zeros(3,ret_gt.cams(2).noviews); % camera target
cu = zeros(3,ret_gt.cams(2).noviews); % camera up vector
for k = 1:ret_gt.cams(2).noviews
    R = func_rodrigues(ret_gt.cams(2).views_orien(:,k));
    ct(:,k) = R(3,:);
    cu(:,k) = R(2,:);
end
cp_ip = (LL*cp')';
cu_ip = (LL*cu')';
ct_ip = (LL*ct')';

plot3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), '-ob');
quiver3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), cu_ip(1,:), cu_ip(2,:),cu_ip(3,:),1, 'g');
quiver3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), cu(1,:), cu(2,:),cu(3,:),1, 'r');
%quiver3(tr_ip(1,:),tr_ip(2,:),tr_ip(3,:), ct_ip(1,:), ct_ip(2,:),ct_ip(3,:),1, 'g');
%quiver3(tr_ip(1,:),tr_ip(2,:),tr_ip(3,:), ct(1,:), ct(2,:),ct(3,:),1, 'r');

% get noise data set
[ret_test, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, rndopt, [0 1 0 0],nocams);
ret_test.points.fixPosition(:) = 1; % fix 3D positions
ret_test.cams(1).fixInternals = uint8([1 1 1]); % fix camera 1 
ret_test.cams(1).fixOrientation = uint8(ones(1,ret_test.cams(1).noviews));
ret_test.cams(1).fixTranslation = uint8(ones(1,ret_test.cams(1).noviews));
ret_test.cams(3).fixInternals = uint8([1 1 1]); % fix camera 3
ret_test.cams(3).fixOrientation = uint8(ones(1,ret_test.cams(3).noviews));
ret_test.cams(3).fixTranslation = uint8(ones(1,ret_test.cams(3).noviews));
ret_test.cams(2).fixInternals = uint8([1 1 1]); % fix internals of camera 2

axistmp = abs(maxaxis(:,1)-maxaxis(:,2))/8;
maxaxis = [maxaxis(:,1)-axistmp  maxaxis(:,2)+axistmp];
figure;
func_plot_dataset(ret_test, [1 0 1 1 0], wh, ccent)
axis([maxaxis(1,1) maxaxis(1,2) maxaxis(2,1) maxaxis(2,2) maxaxis(3,1) maxaxis(3,2)])        
set(gcf, 'Position', [1920 1 1280 1000])
drawnow;
pause(.01)

cp = ret_test.cams(2).views_trans;
ct = zeros(3,ret_test.cams(2).noviews); % camera target
cu = zeros(3,ret_test.cams(2).noviews); % camera up vector
for k = 1:ret_test.cams(2).noviews
    R = func_rodrigues(ret_test.cams(2).views_orien(:,k));
    ct(:,k) = R(3,:);
    cu(:,k) = R(2,:);
end
cp_ip = (LL*cp')';
cu_ip = (LL*cu')';
ct_ip = (LL*ct')';

plot3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), '-og', 'linewidth',2);

% optimize without smoothing constraint
ret_test.BAopt.maxiter = uint64(300);
ret_test.BAopt.timeout_sec = uint64(120);
tic
result = BAdjustMex(ret_test);
toc
figure;
func_plot_dataset(result, [1 0 1 1 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])

cp = result.cams(2).views_trans;
ct = zeros(3,result.cams(2).noviews); % camera target
cu = zeros(3,result.cams(2).noviews); % camera up vector
for k = 1:result.cams(2).noviews
    R = func_rodrigues(result.cams(2).views_orien(:,k));
    ct(:,k) = R(3,:);
    cu(:,k) = R(2,:);
end
cp_ip = (LL*cp')';
cu_ip = (LL*cu')';
ct_ip = (LL*ct')';

plot3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), '-ob');
plot3(cp(1,:),cp(2,:),cp(3,:), '-or');

% optimize with smoothing constraint, use spline smoother matrix, elements == 0 encode fully independent cameras
ret_test_sm = ret_test; % refine result with spline
LLmap = diag(ones(1,n));
for l = 1:1 % how many off-diagonal elements should be used
    LLmap = LLmap +diag(ones(1,n-l),-l)+diag(ones(1,n-l),+l);
end
LL_filt = LL.*LLmap;
LL_filt = LL_filt ./ repmat(sum(LL_filt,1),n,1); % filter smoother matrix, and renormalize
%ret_test_sm.cams(2).smootherM = double(LL_filt(25,1:(end-1))'); % select row, but pass as column vector, make odd length
%ret_test_sm.cams(2).smootherM = double([.5 1 .5]'); % select row, but pass as column vector, make odd length
%ret_test_sm.cams(2).smootherM = double(LL_filt'); % pass transposed
%ret_test_sm.cams(2).smootherM = double(eye(size(LL,1))); % pass transposed
%ret_test_sm.cams(2).fixTranslation(1:2:end) = uint8([ones(1,floor(ret_test_sm.cams(2).noviews/2))]);%uint8([0 0 0 0 0]);
%ret_test_sm.cams(2).fixOrientation(1:2:end) = uint8([ones(1,floor(ret_test_sm.cams(2).noviews/2))]);%uint8([0 0 0 0 0]);

K_filt = K .* LLmap;
normfull = sum(K_filt,1);
normabs = sum(abs(K_filt),1);
K_filt = K_filt - (abs(K_filt)./ repmat(normabs,length(normabs),1)) .*  repmat(normfull,length(normfull),1);
ret_test_sm.cams(2).smootherM = double(1000*K_filt');



ret_test_sm.BAopt.maxiter = uint64(300);
ret_test_sm.BAopt.timeout_sec = uint64(120);
tic
result2 = BAdjustMex(ret_test_sm);
toc
%func_readwrite_BA_problem(ret_test_sm, 'BAproblem_test.txt');  % save org. problem to file
%system('../build/ceresBAbin   BAproblem_test.txt  out.txt');  % call BA in command line, non-mex version
%system('valgrind ../build/ceresBAbin   BAproblem_test.txt  out.txt');  % call in command line, non-mex version, use valgrind for memory leaks
figure;
func_plot_dataset(result2, [1 0 1 1 0], wh, ccent)
set(gcf, 'Position', [1920 1 1280 1000])
cp = result2.cams(2).views_trans;
cp_ip = (LL*cp')';
plot3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), '-ob');
plot3(cp(1,:),cp(2,:),cp(3,:), '-or');

sum(sum((K * result2.cams(2).views_trans').^2,2))
sum(sum((K * ret_test_sm.cams(2).views_trans').^2,2))
sum(sum((K * ret_test.cams(2).views_trans').^2,2))


%     cp = result.cams(2).views_trans;
%     ct = zeros(3,result.cams(2).noviews); % camera target
%     cu = zeros(3,result.cams(2).noviews); % camera up vector
%     for k = 1:result.cams(2).noviews
%         R = func_rodrigues(result.cams(2).views_orien(:,k));
%         ct(:,k) = R(3,:);
%         cu(:,k) = R(2,:);
%     end
%     cp_ip = (LL_filt*cp')';
%     cu_ip = (LL_filt*cu')';
%     ct_ip = (LL_filt*ct')';
%     RR = cell(1,result.cams(2).noviews);
%     for k = 1:result.cams(2).noviews
%         camline = ct_ip(:,k);
%         vecb = cu_ip(:,k);
%         veca = cross(vecb, camline)';
%         vecb = cross(camline,veca);
%         veca = veca ./ norm(veca);
%         vecb = vecb ./ norm(vecb);
%         camline = camline' ./ norm(camline);
%         result.cams(2).views_orien(:,k) = double(func_rodrigues([veca ;  vecb ; camline]));
%     end
%     result.cams(2).views_trans =  double(cp_ip);
% 
% figure;
% func_plot_dataset(result, [1 0 1 1 0], wh, ccent)
% set(gcf, 'Position', [1920 1 1280 1000])
% plot3(cp_ip(1,:),cp_ip(2,:),cp_ip(3,:), '-ob');


pointno = 1;
reprojno = 8;
pt3d = ret_test_sm.points.pt3d(:,pointno);
pt2d = ret_test_sm.points.reproj_pos{pointno}(:,reprojno);
camid = camind(1,ret_test_sm.points.reproj_view{pointno}(:,reprojno));
viewid = camind(2,ret_test_sm.points.reproj_view{pointno}(:,reprojno));
%ret_test_sm.cams(camid).views_trans(:,viewid:(viewid+1))
%val=double([0 1 2 3 4 5 4 3 2 1 0])';%LL_filt(1,1:8);
val=double(LL_filt(1,1:8));
val = val ./ sum(val);
rr_ct = [ 0 0 0];
rr_cu = [ 0 0 0];
t = [ 0 0 0];
for l = 1:length(val)
    R = ret_test_sm.cams(camid).views_orien(:,viewid+l-1);
    R = func_rodrigues(R);
    rr_ct = rr_ct + val(l)*R(3,:);
    rr_cu = rr_cu + val(l)*R(2,:);
    t = t + val(l)*ret_test_sm.cams(camid).views_trans(:,viewid+l-1)';
end
rr_cc =  cross(rr_cu, rr_ct);
rr_cu =  cross(rr_ct, rr_cc);
rr_cc = rr_cc ./ norm(rr_cc);
rr_ct = rr_ct ./ norm(rr_ct);
rr_cu = rr_cu ./ norm(rr_cu);


R = [rr_cc; rr_cu; rr_ct];
t
R


func_reproject_BAtest(pt3d, R,t',ret_test_sm.cams(camid).fc,ret_test_sm.cams(camid).cc,ret_test_sm.cams(camid).kc, 1) - pt2d




%%  Convergence test
numitertotal=10;

% only pt
ptonly=zeros(2,numitertotal); % residual and time
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandSeed, rndopt, [0 0 0 0], nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, ptonly(1,numiter)]= BAdjustMex(ret_test);
    ptonly(2,numiter) = toc;
    ptonly(1,:) = max(ptonly(1,:),1e-1);
end


% pt and pt-plane
ptplane=zeros(2,numitertotal); 
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [1 0 0 0],nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    ret_test.planeconstraint.planew(:)=1e5
    [~, ptplane(1,numiter)]= BAdjustMex(ret_test);
    ptplane(2,numiter) = toc;
    ptplane(1,:) = max(ptplane(1,:),1e-1);
end

% pt and cam-plane
ptcamplane=zeros(2,numitertotal);
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [0 1 0 0],nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, ptcamplane(1,numiter)]= BAdjustMex(ret_test);
    ptcamplane(2,numiter) = toc;
    ptcamplane(1,:) = max(ptcamplane(1,:),1e-1);
end


% pt and mutual visibility
ptvis=zeros(2,numitertotal); 
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [0 0 1 0],nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, ptvis(1,numiter)]= BAdjustMex(ret_test);
    ptvis(2,numiter) = toc;
    ptvis(1,:) = max(ptvis(1,:),1e-1);
end

% pt and vp - fixed internals
ptvpfix=zeros(2,numitertotal); 
for numiter = 1:numitertotal
    tic;
    rndopt_t = rndopt; rndopt_t.fc = [0 0 0]; rndopt_t.cc = [0 0 0]; rndopt_t.kc = [0 0 0];  % was: 1
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt_t, [0 0 0 1],nocams);
    ret_test.cams(1).fixInternals(:)=1; ret_test.cams(2).fixInternals(:)=1; ret_test.cams(3).fixInternals(:)=1;
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, ptvpfix(1,numiter)]= BAdjustMex(ret_test);
    ptvpfix(2,numiter) = toc;
    ptvpfix(1,:) = max(ptvpfix(1,:),1e-1);
end

% pt and vp - free internals
ptvpfree=zeros(2,numitertotal); 
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [0 0 0 1],nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, ptvpfree(1,numiter)]= BAdjustMex(ret_test);
    ptvpfree(2,numiter) = toc;
    ptvpfree(1,:) = max(ptvpfree(1,:),1e-1);
end


% everything
scall=zeros(2,numitertotal); 
for numiter = 1:numitertotal
    tic;
    [ret_test, wh, ccent, maxaxis] = func_create_dataset(RandSeed, rndopt, [1 1 1 1],nocams);
    ret_test.BAopt.maxiter = uint64(numiter); ret_test.BAopt.timeout_sec  = uint64(1000);
    [~, scall(1,numiter)]= BAdjustMex(ret_test);
    scall(2,numiter) = toc;
    scall(1,:) = max(scall(1,:),1e-1);
end



semilogy(1:numitertotal, ptonly(1,:), '-b', 'linewidth',2); hold on
semilogy(1:numitertotal, ptplane(1,:), '-k', 'linewidth',2); 
semilogy(1:numitertotal, ptcamplane(1,:), '-c', 'linewidth',2); 
semilogy(1:numitertotal, ptvis(1,:), '-g', 'linewidth',2); 
semilogy(1:numitertotal, ptvpfix(1,:), '-m', 'linewidth',2);
semilogy(1:numitertotal, ptvpfree(1,:), '-y', 'linewidth',2);
semilogy(1:numitertotal, scall(1,:), '-r', 'linewidth',2); 
grid on
set(gca,'xlim',[1 numitertotal]);
set(gcf,'color','w');
title('\bf BA residuals over time for multiple BA variants')
xlabel('\bf Iter.number');
ylabel('\bf log(residual)');
legend('pts only','pts + pt.plane','pts + cam.plane','pts + mutual.vis', 'pts + vp-fix', 'pts + vp-free', 'all');
set(gcf, 'Position', [1 1 1024 600])


%% TEST for matlab writing script: struct_in -> text file, text file -> struct_out, assert struct_in == struct_out
% func_readwrite_BA_problem(ret_test, '/tmp/BAproblem_test.txt');  % save org. problem to file
% stcttxt = func_readwrite_BA_problem([], '/tmp/BAproblem_test.txt');    % read result
% [match, er1, er2, erc]  = comp_struct(stcttxt, ret_test,2);   % compare both structs, should be exactly the same
% 
% func_readwrite_BA_problem(ret_test, '/tmp/BAproblem_test.txt');  % save org. problem to file
% system('../build/ceresBAbin   /tmp/BAproblem_test.txt  /tmp/out.txt');  % call BA in command line, non-mex version
% stcttxt = func_readwrite_BA_problem([], '/tmp/out.txt');    % read result
% [match, er1, er2, erc]  = comp_struct(result, stcttxt,2);   % compare both structs, should be exactly the same

