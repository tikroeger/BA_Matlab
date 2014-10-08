function [ret, wh, ccent, maxaxis, camind, camptidx] = func_create_dataset(RandStreamSeed, rndopt, opt, nocams, opt_nopts)

if (nargin < 5)
    opt_nopts = 10;
end
%opt(1): use planes for points
%opt(2): use camera planes
%opt(3): use camera mutual visibility
%opt(4): use vanishing point constraints


camplaneweight = 10000; % model scale dependent, if too small, ineffective, too large: overtakes reprojection residuals
vpweight = 10; % VP image-base end-point error (in px), scaling
cammutualvisibilityweight = 10;


s = RandStream('mrg32k3a', 'Seed', RandStreamSeed);
RandStream.setGlobalStream(s);

% number of views in camera blocks
no_c1=nocams(1);
no_c2=nocams(2);
no_c3=nocams(3); 

[randrot ~] = qr(randi(10,3)); % random rotation for cube

% set ground truth camera parameters to these values
fc = [1200 1200];
fcfact1  = 1; %factor for ring1
fcfact2  = 1.4; %factor for ring2
fcfact3  = 3; %factor for static cams
cc = [20 20];
kc = [-.04 .02]; % 2nd and 4th order radial distortion
%kc = [0 0]; % 2nd and 4th order radial distortion
wh = [2000 1200];

% Construct cube
ccent = [124 345 567];
cubesize = 20;
nopoints = opt_nopts; %15;

[a b] = meshgrid(linspace(1, cubesize,nopoints), linspace(1, cubesize,nopoints));
a = a(:); b = b(:);

face{1} = [a-cubesize/2,b-cubesize/2, repmat(cubesize/2,length(a),1)];
face{2} = [a-cubesize/2,b-cubesize/2, -repmat(cubesize/2,length(a),1)+1];
face{3} = [a-cubesize/2, repmat(cubesize/2,length(a),1), b-cubesize/2];
face{4} = [a-cubesize/2,-repmat(cubesize/2,length(a),1)+1, b-cubesize/2];
face{5} = [repmat(cubesize/2,length(a),1),a-cubesize/2, b-cubesize/2];
face{6} = [-repmat(cubesize/2,length(a),1)+1,a-cubesize/2, b-cubesize/2];

for i = 1:6
    face{i} = bsxfun(@plus, (randrot * face{i}')', ccent); % random rotation, recentering
end

allpt = cat(1,face{1},face{2},face{3}, face{4}, face{5},face{6});
allpt = unique(allpt, 'rows');

ptface{1} = find(ismember(allpt, face{1}, 'rows'));
ptface{2} = find(ismember(allpt, face{2}, 'rows'));
ptface{3} = find(ismember(allpt, face{3}, 'rows'));
ptface{4} = find(ismember(allpt, face{4}, 'rows'));
ptface{5} = find(ismember(allpt, face{5}, 'rows'));
ptface{6} = find(ismember(allpt, face{6}, 'rows'));

pt_edge13 = intersect(ptface{1},ptface{3});
pt_edge14 = intersect(ptface{1},ptface{4});
pt_edge15 = intersect(ptface{1},ptface{5});
pt_edge16 = intersect(ptface{1},ptface{6});

pt_corner135 = intersect(pt_edge13,pt_edge15);
pt_corner136 = intersect(pt_edge13,pt_edge16);
pt_corner145 = intersect(pt_edge14,pt_edge15);
pt_corner146 = intersect(pt_edge14,pt_edge16);

pt_edge23 = intersect(ptface{2},ptface{3});
pt_edge24 = intersect(ptface{2},ptface{4});
pt_edge25 = intersect(ptface{2},ptface{5});
pt_edge26 = intersect(ptface{2},ptface{6});

pt_corner235 = intersect(pt_edge23,pt_edge25);
pt_corner236 = intersect(pt_edge23,pt_edge26);
pt_corner245 = intersect(pt_edge24,pt_edge25);
pt_corner246 = intersect(pt_edge24,pt_edge26);

pt_edge35 = intersect(ptface{3},ptface{5});
pt_edge36 = intersect(ptface{3},ptface{6});
pt_edge45 = intersect(ptface{4},ptface{5});
pt_edge46 = intersect(ptface{4},ptface{6});

allcorners = [pt_corner135; pt_corner136; pt_corner145; pt_corner146; pt_corner235; pt_corner236; pt_corner245; pt_corner246];

for i = 1:6
    facecorners{i} = intersect(allcorners, ptface{i});
end


% compute Vanishing points as 3D directions for all cube edges
totalvps = 2;
VP_cornerids = cell(1,totalvps);
VP_GT_dir = cell(1,totalvps);
VP_cornerids{1} = allcorners([[1 2];[3 4]; [5 6]; [7 8]]);
VP_cornerids{2} = allcorners([[1 3]; [2 4]; [5 7]; [6 8];]);
VP_cornerids{3} = allcorners([[1 5]; [2 6]; [3 7]; [4 8];]);
for vpid = 1:totalvps
    tmp = allpt(VP_cornerids{vpid}(1,1),:)-allpt(VP_cornerids{vpid}(1,2),:);
    VP_GT_dir{vpid} = tmp ./ norm(tmp);
    if (~isempty(rndopt)) % add noise if randomized toy dataset is required 
        ang_rad = 2*pi*(rndopt.randsigma_vpang/360);
        rat = rand([1 3]);
        rat = rat ./ sum(rat);
        Rrand = func_get_rot_matrix(rat(1)*ang_rad,rat(2)*ang_rad,rat(3)*ang_rad);
        npnew = (Rrand*VP_GT_dir{vpid}')';
        %acosd(dot(npnew,VP_GT_dir{vpid}))
        VP_GT_dir{vpid} = npnew;
    end
    
end
% % plot
% scatter3(allpt(allcorners,1), allpt(allcorners,2),allpt(allcorners,3), 'b','filled'); hold on
% lines = [allpt(VP_cornerids{1}(:,1),:) allpt(VP_cornerids{1}(:,2),:)];
% plot3([lines(:,1) lines(:,4)]', [lines(:,2) lines(:,5)]', [lines(:,3) lines(:,6)]', '-r')
% lines = [allpt(VP_cornerids{2}(:,1),:) allpt(VP_cornerids{2}(:,2),:)];
% plot3([lines(:,1) lines(:,4)]', [lines(:,2) lines(:,5)]', [lines(:,3) lines(:,6)]', '-g')
% lines = [allpt(VP_cornerids{3}(:,1),:) allpt(VP_cornerids{3}(:,2),:)];
% plot3([lines(:,1) lines(:,4)]', [lines(:,2) lines(:,5)]', [lines(:,3) lines(:,6)]', '-b')
% quiver3(ccent(1),ccent(2),ccent(3), VP_GT_dir{1}(1),VP_GT_dir{1}(2),VP_GT_dir{1}(3), 5, 'r')
% quiver3(ccent(1),ccent(2),ccent(3), VP_GT_dir{2}(1),VP_GT_dir{2}(2),VP_GT_dir{2}(3), 5, 'g')
% quiver3(ccent(1),ccent(2),ccent(3), VP_GT_dir{3}(1),VP_GT_dir{3}(2),VP_GT_dir{3}(3), 5, 'b')
% axis equal


% Compute planes
for i = 1:6
    [coeff,score,latent] = pca(face{i});
    [~, mnidx] = min(latent);
    planevec{i} = coeff(:,mnidx);
    planevec{i} = planevec{i} * (mean(face{i}) * planevec{i});

    mf = mean(face{i});
    dmfp = mf*planevec{i}./norm(planevec{i}) - norm(planevec{i}); % distance of mean of pts on face from plane
    mf = mf - (dmfp*planevec{i}./norm(planevec{i}))';  % mean projected on plane

    diff = sum(bsxfun(@times, allpt(facecorners{i},:), (planevec{i}./norm(planevec{i}))'),2) - norm(planevec{i}); % distance of corner points on plane
    cornersprojected = allpt(facecorners{i},:) - bsxfun(@times, diff, (planevec{i}./norm(planevec{i}))') ; % project corners onto plane
    cornersprojected = cornersprojected + bsxfun(@minus, cornersprojected, mf) / 2; % make boundary .25 * cubesize larger for visualization purpose
end


% Generate cameras
ringvec1 = rand([1 3])-.5;
ringvec1= ringvec1./norm(ringvec1);
ringvec2 = cross(ringvec1, ringvec1+rand([1 3])-.5);
ringvec2 = ringvec2./norm(ringvec2);
ringvec1 = ringvec1;
ringvec2 = ringvec2;

ring1campos = zeros(3,no_c1);
ring1camrot = zeros(3,no_c1);
vec = cross(ringvec1, (rand([1 3])-.5));
for k = 1:no_c1
    ang = 360/no_c1*(k-1);
    resvec = vec*cosd(ang)+ cross(vec,ringvec1)*sind(ang) + ringvec1*(sum(ringvec1.*vec))*(1-cosd(ang));% rodrigues rotation
    ring1campos(:,k) = resvec ./ norm(resvec);
end
% ring1campos(:,1) = vec ./ norm(vec);
% ring1campos(:,3) = cross(vec,ringvec1./norm(ringvec1));
% ring1campos(:,3) = ring1campos(:,3)./norm(ring1campos(:,3));
% ring1campos(:,5) = -ring1campos(:,1);
% ring1campos(:,7) = -ring1campos(:,3);
% ring1campos(:,2) = cosd(45)*ring1campos(:,1) + sind(45)*ring1campos(:,3);
% ring1campos(:,4) = cosd(45)*ring1campos(:,3) + sind(45)*ring1campos(:,5);
% ring1campos(:,6) = cosd(45)*ring1campos(:,5) + sind(45)*ring1campos(:,7);
% ring1campos(:,8) = cosd(45)*ring1campos(:,7) + sind(45)*ring1campos(:,1);

ring2campos = zeros(3,no_c2);
ring2camrot = zeros(3,no_c2);
vec = cross(ringvec2, (rand([1 3])-.5));
for k = 1:no_c2
    ang = 360/no_c2*(k-1);
    resvec = vec*cosd(ang)+ cross(vec,ringvec2)*sind(ang) + ringvec2*(sum(ringvec2.*vec))*(1-cosd(ang));% rodrigues rotation
    ring2campos(:,k) = resvec ./ norm(resvec);
end

% ring2campos(:,1) = vec ./ norm(vec);
% ring2campos(:,3) = cross(vec,ringvec2./norm(ringvec2));
% ring2campos(:,3) = ring2campos(:,3)./norm(ring2campos(:,3));
% ring2campos(:,5) = -ring2campos(:,1);
% ring2campos(:,7) = -ring2campos(:,3);
% ring2campos(:,2) = cosd(45)*ring2campos(:,1) + sind(45)*ring2campos(:,3);
% ring2campos(:,4) = cosd(45)*ring2campos(:,3) + sind(45)*ring2campos(:,5);
% ring2campos(:,6) = cosd(45)*ring2campos(:,5) + sind(45)*ring2campos(:,7);
% ring2campos(:,8) = cosd(45)*ring2campos(:,7) + sind(45)*ring2campos(:,1);

ring1campos = bsxfun(@plus, ring1campos * cubesize * 2.1, ccent');  % distances of cameras from cube center
ring2campos = bsxfun(@plus, ring2campos * cubesize * 3.2, ccent');
ringvec1 = ringvec1 * (ccent * ringvec1');  % scale plane normal of ring
ringvec2 = ringvec2 * (ccent * ringvec2');

for i = 1:no_c1
    vecnorm = (ccent'-ring1campos(:,i));
    vecnorm = vecnorm./norm(vecnorm);
    %vecnorm_orth = cross(vecnorm, vecnorm + rand([3 1])-.5);
    %vecnorm_orth = vecnorm_orth./norm(vecnorm_orth);
    %vecnorm_orth2 = cross(vecnorm, vecnorm_orth);
    %vecnorm_orth = -cross(vecnorm, vecnorm_orth2);
    
    ang = (4*pi/no_c1)*(i-1); % do 2 full rolls
    vecnorm_orth = ringvec1';
    vecnorm_orth = vecnorm_orth*cos(ang)+ cross(vecnorm_orth,vecnorm)*sin(ang) + vecnorm*(sum(vecnorm.*vecnorm_orth))*(1-cos(ang)); % rodrigues rotation

    vecnorm_orth = vecnorm_orth./norm(vecnorm_orth);
    vecnorm_orth2 = -cross(vecnorm, vecnorm_orth);
    %vecnorm_orth = -cross(vecnorm, vecnorm_orth2);
    R = [vecnorm_orth2 vecnorm_orth vecnorm]';
    Rrand = func_get_rot_matrix((rand-.5) * pi/180*20,(rand-.5) * pi/180*20,(rand-.5) * pi/180*20);
    R = (R*Rrand);
    ring1camrot(:,i) = func_rodrigues(R);
end

for i = 1:no_c2
    vecnorm = (ccent'-ring2campos(:,i));
    vecnorm = vecnorm./norm(vecnorm);
    %vecnorm_orth = cross(vecnorm, vecnorm + rand([3 1])-.5);
    %vecnorm_orth = vecnorm_orth./norm(vecnorm_orth);
    %vecnorm_orth2 = cross(vecnorm, vecnorm_orth);
    %vecnorm_orth = -cross(vecnorm, vecnorm_orth2);
    
    ang = (4*pi/no_c2)*(i-1); % do 2 full rolls
    vecnorm_orth = ringvec2';
    vecnorm_orth = vecnorm_orth*cos(ang)+ cross(vecnorm_orth,vecnorm)*sin(ang) + vecnorm*(sum(vecnorm.*vecnorm_orth))*(1-cos(ang)); % rodrigues rotation
    vecnorm_orth = vecnorm_orth./norm(vecnorm_orth);
    vecnorm_orth2 = -cross(vecnorm, vecnorm_orth);    
    R = [vecnorm_orth2 vecnorm_orth vecnorm]';
    Rrand = func_get_rot_matrix((rand-.5) * pi/180*20,(rand-.5) * pi/180*20,(rand-.5) * pi/180*20);
    R =(R*Rrand);    
    ring2camrot(:,i) = func_rodrigues(R);
end

% Generate static camera with no_c3 views
staticcampos = ccent + norm([3 1])* cubesize/2 * 2;
staticcamrot = zeros(3,no_c3);
for i = 1:no_c3
    vecnorm = (ccent'-staticcampos');
    vecnorm = vecnorm./norm(vecnorm);
    vecnorm_orth = cross(vecnorm, vecnorm + rand([3 1])-.5);
    vecnorm_orth = vecnorm_orth./norm(vecnorm_orth);
    vecnorm_orth2 = cross(vecnorm, vecnorm_orth);
    vecnorm_orth = -cross(vecnorm, vecnorm_orth2);
    R = [vecnorm_orth vecnorm_orth2 vecnorm]';
    Rrand = func_get_rot_matrix((rand-.5) * pi/180*20,(rand-.5) * pi/180*20,(rand-.5) * pi/180*20);
    R = (R*Rrand);        
    staticcamrot(:,i) = func_rodrigues(R);
end



 
%sum(bsxfun(@times, ring1campos, (ringvec1./norm(ringvec1))'),1) - norm(ringvec1)

% Collect visible points in each camera, collect mutual visibility of cameras on ring1
ptvisiblering1 = cell(1,size(ring1camrot,2)) ;
linesvisiblering1 = cell(totalvps,size(ring1camrot,2));
camvisiblering1 = cell(1,size(ring1camrot,2)); % cameras visible from ring1, ring2, static
for i = 1:no_c1
    fct = fc*fcfact1;    
    R = func_rodrigues(ring1camrot(:,i));
    t = ring1campos(:,i);
    pt2d_repro = zeros(3,size(allpt,1));

    %pt2d_repro = arrayfun(@(k) func_reproject_BAtest(allpt(k,:)', R,t,fct,cc,kc,1), 1:size(allpt,1), 'UniformOutput',0);
    %pt2d_repro = arrayfun(@(k) func_reproject_BAtest(allpt(k,:)', R,-R*t,fct,cc,kc,0), 1:size(allpt,1), 'UniformOutput',0); % 
    %pt2d_repro = cat(1,pt2d_repro{:});
    
    pt2d_repro = func_reproject_BAtest(allpt, R,-R*t,fct,cc,kc,0)';

    pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);
    idx = find(pt2d_repro(:,1) <= 0 | pt2d_repro(:,2) <= 0 | pt2d_repro(:,1) > wh(1) | pt2d_repro(:,2) > wh(2));
    pt2d_repro = bsxfun(@minus, pt2d_repro, [wh(1)/2 wh(2)/2]);
    %pt2d_repro(idx,:) = [];
    
    % compute VP lines    
    for vpid = 1:totalvps;
        vpids_tmp = VP_cornerids{vpid}(sum(ismember(VP_cornerids{vpid},idx),2)==0,:); % remove lines for which pt lies outside of image frustum
        lines_vp = [pt2d_repro(vpids_tmp(:,1),:) pt2d_repro(vpids_tmp(:,2),:)]; % 2D lines for
        %direc_vp = R*VP_GT_dir{vpid}';   % 3D direction for VP, in camera reference frame
        linesvisiblering1{vpid,i} = lines_vp;
    end
%     lines = [pt2d_repro(vpids_tmp(:,1),:) pt2d_repro(vpids_tmp(:,2),:)]; % 2D lines for
%     lines_re = bsxfun(@minus, lines, [cc cc]);
%     lines_re = bsxfun(@rdivide, lines_re, [fc fc]);
%     lineweight = (lines_re(:,1)-lines_re(:,3)).^2+(lines_re(:,2)-lines_re(:,4)).^2; % squared line length
%     lines_re = [lines_re(:,1:2) ones(size(lines_re,1),1) lines_re(:,3:4) ones(size(lines_re,1),1)];
%     a=sqrt(sum(lines_re(:,1:3).^2,2));
%     b=sqrt(sum(lines_re(:,4:6).^2,2));
%     lines_re = bsxfun(@rdivide, lines_re, [a a a b b b]);
%     % interpretation planes
%     planevec = cross(lines_re(:,1:3),lines_re(:,4:6));
%     planevec = bsxfun(@rdivide, planevec, sqrt(sum(planevec.^2,2)));
%     func_SVD_compute_VP(planevec, lineweight)
%     scatter(pt2d_repro(:,1),pt2d_repro(:,2)); hold on
%     plot([lines(:,1) lines(:,3)]', [lines(:,2) lines(:,4)]', '-r', 'linewidth',2)
%     axis equal
    
    
    %find planes of cube oriented towards camera, add points from plane
	ptvisible = zeros(size(allpt,1),1);
    camvec = R(3,:)./norm(R(3,:));
    for k = 1:6 
        planevectmp{k} = ccent - mean(allpt(ptface{k},:));
        ang = sum(camvec   .*   (planevectmp{k}./norm(planevectmp{k})));
        acosd(ang);
        if (abs(acosd(ang))<90)
            ptvisible(ptface{k})=1;
        end
    end
%                 scatter3(allpt(find(~ptvisible),1), allpt(find(~ptvisible),2), allpt(find(~ptvisible),3), 15, 'filled', 'b'); hold on
%                 scatter3(allpt(find(ptvisible),1), allpt(find(ptvisible),2), allpt(find(ptvisible),3), 15, 'filled', 'r');
%                 K = eye(3);
%                 K(1,1) = fcfact2*fc(1); K(2,2) = fcfact2*fc(2); 
%                 K(1,3) = cc(1); K(2,3) = cc(2);
%                 func_plot_cameras(K,func_rodrigues(ring2camrot(:,i)), (-func_rodrigues(ring2camrot(:,i))*ring2campos(:,i))', wh(1), wh(2), [1 0 0], 15, 1, 1);
%                 axis equal
    
    ptvisiblering1{i} = [ptvisible pt2d_repro allpt (1:size(allpt,1))'];
    ptvisiblering1{i}(idx,:) = [];   % delete those, that lie outside image frustum
    ptvisiblering1{i} = ptvisiblering1{i}(find(ptvisiblering1{i}(:,1)),:);
            
% %  Plot reprojected points
%             figure
%             sc = scatter(pt2d_repro(:,1), pt2d_repro(:,2), 25, 'filled', 'b'); hold on
%             sc = scatter(pt2d_repro(ptvisible,1), pt2d_repro(ptvisible,2), 25, 'filled', 'r');
%             xlim([1 wh(1)])
%             ylim([1 wh(2)])
% %             

    K = eye(3);
    K(1,1) = fct(1); K(2,2) = fct(2); 
    K(1,3) = cc(1); K(2,3) = cc(2);  
    camvisiblering1{i} = [];
    for k = 1:no_c1
        if (k~=i)
            pt2d_repro = func_reproject_BAtest(ring1campos(:,k), R,-R*t,fct,cc,kc,0)';
            pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);            
            if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring1campos(:,k)' 1]',K*[R -R*t]) > 0)
              camvisiblering1{i} = cat(1, camvisiblering1{i},[1 k pt2d_repro-wh/2]);
            end
        end
    end
    for k = 1:no_c2
        pt2d_repro = func_reproject_BAtest(ring2campos(:,k), R,-R*t,fct,cc,kc,0)';
        pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
        if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring2campos(:,k)' 1]',K*[R -R*t]) > 0)
          camvisiblering1{i} = cat(1, camvisiblering1{i},[2 k pt2d_repro-wh/2]);
        end
    end    
    for k = 1:no_c3
        pt2d_repro = func_reproject_BAtest(staticcampos', R,-R*t,fct,cc,kc,0)';
        pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
        if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([staticcampos 1]',K*[R -R*t]) > 0)
          camvisiblering1{i} = cat(1, camvisiblering1{i},[3 k pt2d_repro-wh/2]);
        end
    end        
end


% Collect visible points in each camera, collect mutual visibility of cameras on ring2
ptvisiblering2 = cell(1,size(ring2camrot,2));
linesvisiblering2 = cell(totalvps,size(ring2camrot,2));
camvisiblering2 = cell(1,size(ring2camrot,2)); % cameras visible from ring1, ring2, static
for i = 1:no_c2
    fct = fc*fcfact2;    
    R = func_rodrigues(ring2camrot(:,i));
    t = ring2campos(:,i);
    pt2d_repro = zeros(3,size(allpt,1));
    %pt2d_repro = arrayfun(@(k) func_reproject_BAtest(allpt(k,:)', R,-R*t,fct,cc,kc,0), 1:size(allpt,1), 'UniformOutput',0); % 
    %pt2d_repro = cat(1,pt2d_repro{:});
    pt2d_repro = func_reproject_BAtest(allpt, R,-R*t,fct,cc,kc,0)';
    func_reproject_BAtest(allpt(1,:), R,-R*t,fct,cc,kc,0);
    pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);
    idx = find(pt2d_repro(:,1) <= 0 | pt2d_repro(:,2) <= 0 | pt2d_repro(:,1) > wh(1) | pt2d_repro(:,2) > wh(2));
    pt2d_repro = bsxfun(@minus, pt2d_repro, [wh(1)/2 wh(2)/2]);    
    %pt2d_repro(idx,:) = [];
    
    % compute VP lines    
    for vpid = 1:totalvps;
        vpids_tmp = VP_cornerids{vpid}(sum(ismember(VP_cornerids{vpid},idx),2)==0,:); % remove lines for which pt lies outside of image frustum
        lines_vp = [pt2d_repro(vpids_tmp(:,1),:) pt2d_repro(vpids_tmp(:,2),:)]; % 2D lines for
        %direc_vp = R*VP_GT_dir{vpid}';   % 3D direction for VP, in camera reference frame
        linesvisiblering2{vpid,i} = lines_vp;
    end

    %find planes of cube oriented towards camera, add points from plane
	ptvisible = zeros(size(allpt,1),1);
    camvec = R(3,:)./norm(R(3,:));
    for k = 1:6 
        planevectmp{k} = ccent - mean(allpt(ptface{k},:));
        ang = sum(camvec   .*   (planevectmp{k}./norm(planevectmp{k})));
        if (abs(acosd(ang))<90)
            ptvisible(ptface{k})=1;
        end
    end
    ptvisiblering2{i} = [ptvisible pt2d_repro allpt (1:size(allpt,1))'];
    ptvisiblering2{i}(idx,:) = [];   % delete those, that lie outside image frustum
    ptvisiblering2{i} = ptvisiblering2{i}(find(ptvisiblering2{i}(:,1)),:);

    K = eye(3);
    K(1,1) = fct(1); K(2,2) = fct(2); 
    K(1,3) = cc(1); K(2,3) = cc(2);  
    camvisiblering2{i} = [];
    for k = 1:no_c1
        pt2d_repro = func_reproject_BAtest(ring1campos(:,k), R,-R*t,fct,cc,kc,0)';
        pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);            
        if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring1campos(:,k)' 1]',K*[R -R*t]) > 0)
          camvisiblering2{i} = cat(1, camvisiblering2{i},[1 k pt2d_repro-wh/2]);
        end
    end
    for k = 1:no_c2
        if (k~=i)        
            pt2d_repro = func_reproject_BAtest(ring2campos(:,k), R,-R*t,fct,cc,kc,0)';
            pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
            if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring2campos(:,k)' 1]',K*[R -R*t]) > 0)
              camvisiblering2{i} = cat(1, camvisiblering2{i},[2 k pt2d_repro-wh/2]);
            end
        end
    end    
    for k = 1:no_c3
        pt2d_repro = func_reproject_BAtest(staticcampos', R,-R*t,fct,cc,kc,0)';
        pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
        if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([staticcampos 1]',K*[R -R*t]) > 0)
          camvisiblering2{i} = cat(1, camvisiblering2{i},[3 k pt2d_repro-wh/2]);
        end
    end        
end


% Collect visible points in each camera, collect mutual visibility of cameras on static cameras
ptvisiblestatic = cell(1,size(staticcamrot,2));
linesvisiblestatic = cell(totalvps,size(staticcamrot,2));
camvisiblestatic = cell(1,size(staticcamrot,2)); % cameras visible from ring1, ring2, static
for i = 1:no_c3
    fct = fc*fcfact3;    
    R = func_rodrigues(staticcamrot(:,i));
    t = staticcampos';
    pt2d_repro = zeros(3,size(allpt,1));
    %pt2d_repro = arrayfun(@(k) func_reproject_BAtest(allpt(k,:)', R,-R*t,fct,cc,kc,0), 1:size(allpt,1), 'UniformOutput',0); % 
    %pt2d_repro = cat(1,pt2d_repro{:});
    pt2d_repro = func_reproject_BAtest(allpt, R,-R*t,fct,cc,kc,0)';
    pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);
    idx = find(pt2d_repro(:,1) <= 0 | pt2d_repro(:,2) <= 0 | pt2d_repro(:,1) > wh(1) | pt2d_repro(:,2) > wh(2));
    pt2d_repro = bsxfun(@minus, pt2d_repro, [wh(1)/2 wh(2)/2]);    
    %pt2d_repro(idx,:) = [];
    
    % compute VP lines    
    for vpid = 1:totalvps;
        vpids_tmp = VP_cornerids{vpid}(sum(ismember(VP_cornerids{vpid},idx),2)==0,:); % remove lines for which pt lies outside of image frustum
        lines_vp = [pt2d_repro(vpids_tmp(:,1),:) pt2d_repro(vpids_tmp(:,2),:)]; % 2D lines for
        %direc_vp = R*VP_GT_dir{vpid}';   % 3D direction for VP, in camera reference frame
        linesvisiblestatic{vpid,i} = lines_vp;
    end
    
    %find planes of cube oriented towards camera, add points from plane
	ptvisible = zeros(size(allpt,1),1);
    camvec = R(3,:)./norm(R(3,:));
    for k = 1:6 
        planevectmp{k} = ccent - mean(allpt(ptface{k},:));
        ang = sum(camvec   .*   (planevectmp{k}./norm(planevectmp{k})));
        if (abs(acosd(ang))<90)
            ptvisible(ptface{k})=1;
        end
    end
    ptvisiblestatic{i} = [ptvisible pt2d_repro allpt (1:size(allpt,1))'];
    ptvisiblestatic{i}(idx,:) = [];   % delete those, that lie outside image frustum
    ptvisiblestatic{i} = ptvisiblestatic{i}(find(ptvisiblestatic{i}(:,1)),:);

    K = eye(3);
    K(1,1) = fct(1); K(2,2) = fct(2); 
    K(1,3) = cc(1); K(2,3) = cc(2);  
    camvisiblestatic{i} = [];
    for k = 1:no_c1
        pt2d_repro = func_reproject_BAtest(ring1campos(:,k), R,-R*t,fct,cc,kc,0)';
        pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);            
        if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring1campos(:,k)' 1]',K*[R -R*t]) > 0)
          camvisiblestatic{i} = cat(1, camvisiblestatic{i},[1 k pt2d_repro-wh/2]);
        end
    end
    for k = 1:no_c2
            pt2d_repro = func_reproject_BAtest(ring2campos(:,k), R,-R*t,fct,cc,kc,0)';
            pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
            if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([ring2campos(:,k)' 1]',K*[R -R*t]) > 0)
              camvisiblestatic{i} = cat(1, camvisiblestatic{i},[2 k pt2d_repro-wh/2]);
            end
    end    
    for k = 1:no_c3
        if (k~=i)                
            pt2d_repro = func_reproject_BAtest(staticcampos', R,-R*t,fct,cc,kc,0)';
            pt2d_repro = bsxfun(@plus, pt2d_repro, [wh(1)/2 wh(2)/2]);                    
            if (pt2d_repro(:,1) > 0 && pt2d_repro(:,2) > 0 && pt2d_repro(:,1) < wh(1) && pt2d_repro(:,2) < wh(2) && func_ptdepth([staticcampos 1]',K*[R -R*t]) > 0)
              camvisiblestatic{i} = cat(1, camvisiblestatic{i},[3 k pt2d_repro-wh/2]);
            end
        end
    end        
end


% fill struct
clear ret;
[status, out] = system('grep -c processor /proc/cpuinfo');
ret.formatversion = double(1.0);
ret.BAopt.cores = uint64(str2num(out));
ret.BAopt.timeout_sec = uint64(20);
ret.BAopt.maxiter = uint64(20);
ret.BAopt.PTLoss = double([0 0]); % use trivial loss
ret.BAopt.PlaneLoss = double([0 0]); % use trivial loss
ret.BAopt.DerivLoss = double([0 0]); % use trivial loss
ret.BAopt.VPLoss = double([0 0]); % use trivial loss

ret.planeconstraint.plane = double([]);
ret.planeconstraint.planew = double([]);
ret.planeconstraint.fixPlane = uint8([]);
if (opt(1)==1)
    ret.planeconstraint.plane = double(cat(2,planevec{:}));
end


if (opt(4)==1)
    ret.VPconstraint.vp = double(cat(1,VP_GT_dir{:})');
    ret.VPconstraint.vpw = double(vpweight*ones(1,size(ret.VPconstraint.vp,2))); % 1xN: VP weights, if empty assumed to be one for all VPs
    ret.VPconstraint.fixVP = uint8(zeros(1,size(ret.VPconstraint.vp,2))); % 1xN: Fix VP or leave it as free variable
else
    ret.VPconstraint.vp = double([]);
    ret.VPconstraint.vpw = double([]);
    ret.VPconstraint.fixVP = uint8([]);
end


 
    ret.cams(1).noviews = uint64(no_c1);
    ret.cams(1).fc = double(fc*fcfact1);
    ret.cams(1).cc = double(cc);
    ret.cams(1).kc = double(kc);
    ret.cams(1).viewids = uint64([1:no_c1]);
    ret.cams(1).fixInternals = uint8([0 0 0]);
    ret.cams(1).camweight = double(ones(1,ret.cams(1).noviews));
    for k = 1:no_c1
        ret.cams(1).views_orien(:,k) = double(ring1camrot(:,k));
        ret.cams(1).views_trans(:,k) = double(ring1campos(:,k));
    end
    if (opt(2)==1)
        if (opt(1)==1)
            ret.cams(1).OnPlane = uint64(7*ones(1,1)); % link to 7th plane (ring1vec) because 6 cube-planes are included as well
        else
            ret.cams(1).OnPlane = uint64(1*ones(1,1)); % link to 1th plane (ring1vec)            
        end
    else
        ret.cams(1).OnPlane = uint64(0*ones(1,1)); % not on plane
    end
    ret.cams(1).fixOrientation = uint8(zeros(1,size(ret.cams(1).views_orien,2))); % do not fix orientation of camera
    ret.cams(1).fixTranslation = uint8(zeros(1,size(ret.cams(1).views_trans,2))); % do not fix translation of camera
    
    ret.cams(1).smootherM = double([]);
    
    ret.cams(1).vplines = double([]);
    if (opt(4)==1)
        for l = 1:ret.cams(1).noviews
            ret.cams(1).vplines(l).l = cell(0, 2);
            cnt=0;
            for k = 1:totalvps 
                if (l~=4)  % leave one frame without VP lines for testing purposes                                
                    if (~isempty(linesvisiblering1{k,l}))
                        cnt = cnt + 1;
                        ret.cams(1).vplines(l).l{cnt,1} = uint64(k);
                        ret.cams(1).vplines(l).l{cnt,2} = linesvisiblering1{k,l};
                    end
                end
            end
        end
    end
    

    if (opt(2)==1)
        ret.planeconstraint.plane = double(cat(2,ret.planeconstraint.plane, ringvec1'));
        %ret.planeconstraint.planew = double([ret.planeconstraint.planew camplaneweight]);
        %ret.planeconstraint.fixPlane = uint8([ret.planeconstraint.fixPlane 0]);
    end

    
    
    ret.cams(2).noviews = uint64(no_c2);
    ret.cams(2).fc = double(fc*fcfact2);
    ret.cams(2).cc = double(cc);
    ret.cams(2).kc = double(kc);
    ret.cams(2).viewids = uint64([(no_c1+1):(no_c1+no_c2)]);
    ret.cams(2).fixInternals = uint8([0 0 0]);
    ret.cams(2).camweight = double(ones(1,ret.cams(2).noviews));
    for k = 1:no_c2
        ret.cams(2).views_orien(:,k) = double(ring2camrot(:,k));
        ret.cams(2).views_trans(:,k) = double(ring2campos(:,k));
    end
    if (opt(2)==1)
        if (opt(1)==1)
            ret.cams(2).OnPlane = uint64(8*ones(1,1)); % link to 8th plane (ring2vec) because 6 cube-planes are included as well
        else
            ret.cams(2).OnPlane = uint64(2*ones(1,1)); % link to 2th plane (ring2vec)            
        end
    else
        ret.cams(2).OnPlane = uint64(0*ones(1,1)); % not on plane
    end    
    ret.cams(2).fixOrientation = uint8(zeros(1,size(ret.cams(2).views_orien,2))); % do not fix orientation of camera
    ret.cams(2).fixTranslation = uint8(zeros(1,size(ret.cams(2).views_trans,2))); % do not fix translation of camera
    
    ret.cams(2).smootherM = double([]);
    
    ret.cams(2).vplines = double([]);
    if (opt(4)==1)
        for l = 1:ret.cams(2).noviews
            ret.cams(2).vplines(l).l = cell(0, 2);
            cnt=0;
            for k = 1:totalvps 
                if (l~=4)  % leave one frame without VP lines for testing purposes                
                    if (~isempty(linesvisiblering2{k,l}))
                        cnt = cnt + 1;
                        ret.cams(2).vplines(l).l{cnt,1} = uint64(k);
                        ret.cams(2).vplines(l).l{cnt,2} = linesvisiblering2{k,l};
                    end
                end
            end
        end
    end
    
    if (opt(2)==1)
        ret.planeconstraint.plane = double(cat(2,ret.planeconstraint.plane, ringvec2'));
        %ret.planeconstraint.planew = double([ret.planeconstraint.planew camplaneweight]);
        %ret.planeconstraint.fixPlane = uint8([ret.planeconstraint.fixPlane 0]);
    end

    
    
    ret.cams(3).noviews = uint64(no_c3);
    ret.cams(3).fc = double(fc*fcfact3);
    ret.cams(3).cc = double(cc);
    ret.cams(3).kc = double(kc);
    ret.cams(3).viewids = uint64([(no_c1+no_c2+1):(no_c1+no_c2+no_c3)]);
    ret.cams(3).fixInternals = uint8([0 0 0]);
    ret.cams(3).camweight = double(ones(1,ret.cams(3).noviews));
    ret.cams(3).views_trans = double(staticcampos');  % shared translation
    for k = 1:no_c3
        ret.cams(3).views_orien(:,k) = double(staticcamrot(:,k));
    end
    ret.cams(3).OnPlane = uint64(zeros(1,1)); % do not link to plane
    ret.cams(3).fixOrientation = uint8(zeros(1,size(ret.cams(3).views_orien,2))); % do not fix orientation of camera
    ret.cams(3).fixTranslation = uint8(zeros(1,size(ret.cams(3).views_trans,2))); % do not fix translation of camera

    ret.cams(3).smootherM = double([]);

    ret.cams(3).vplines = double([]);
    if (opt(4)==1)
        for l = 1:ret.cams(3).noviews
            ret.cams(3).vplines(l).l = cell(0, 2);
            cnt=0;
            for k = 1:totalvps 
                if (l~=4)  % leave one frame without VP lines for testing purposes
                    if (~isempty(linesvisiblestatic{k,l}))
                        cnt = cnt + 1;
                        ret.cams(3).vplines(l).l{cnt,1} = uint64(k);
                        ret.cams(3).vplines(l).l{cnt,2} = linesvisiblestatic{k,l};
                    end
                end
            end
        end
    end
    
    
no3Dpoints = size(allpt,1);
ret.points.pt3d = allpt';
ret.points.fixPosition = uint8(zeros(1, no3Dpoints)); % is point is fixed, or a free variable?, if empty set uniformly to zero
ret.points.pointw = double(ones(1, no3Dpoints)); % points weights, if empty set uniformly to one
ret.points.OnPlane = cell(1, no3Dpoints);
if (opt(1)==1)
    for k = 1:6
        for i = ptface{k}'
            ret.points.OnPlane{i} = uint64([ret.points.OnPlane{i} k]);
        end
    end
end

ret.points.reproj_view = cell(1, no3Dpoints);
ret.points.reproj_pos = cell(1, no3Dpoints);
for i = 1:no_c1
    for k = 1:size(ptvisiblering1{i},1)
        pt3did = ptvisiblering1{i}(k,end);
        pt2d = ptvisiblering1{i}(k,2:3);
        ret.points.reproj_view{pt3did} = uint64([ret.points.reproj_view{pt3did} i+0]);
        ret.points.reproj_pos {pt3did} = double(cat(2,ret.points.reproj_pos {pt3did}, pt2d'));
    end
end
for i = 1:no_c2
    for k = 1:size(ptvisiblering2{i},1)
        pt3did = ptvisiblering2{i}(k,end);
        pt2d = ptvisiblering2{i}(k,2:3);
        ret.points.reproj_view{pt3did} = uint64([ret.points.reproj_view{pt3did} i+no_c1]);
        ret.points.reproj_pos {pt3did} = double(cat(2,ret.points.reproj_pos {pt3did}, pt2d'));
    end
end
for i = 1:no_c3
    for k = 1:size(ptvisiblestatic{i},1)
        pt3did = ptvisiblestatic{i}(k,end);
        pt2d = ptvisiblestatic{i}(k,2:3);
        ret.points.reproj_view{pt3did} = uint64([ret.points.reproj_view{pt3did} i+(no_c1+no_c2)]);
        ret.points.reproj_pos {pt3did} = double(cat(2,ret.points.reproj_pos {pt3did}, pt2d'));
    end
end

if (size(ret.planeconstraint.plane,2)>0)
    ret.planeconstraint.planew = double(camplaneweight*ones(1,size(ret.planeconstraint.plane,2))); % 1xN: Plane weights, if empty assumed to be one for all planes
    ret.planeconstraint.fixPlane = uint8(zeros(1,size(ret.planeconstraint.plane,2))); % 1xN: Fix plane or leave it as free variable
end

ret.cam_reproj.view = uint64([]);
ret.cam_reproj.pos = double([]);
ret.cam_reproj.weight = double([]);
if (opt(3) == 1)
    for i = 1:no_c1
        for k = 1:size(camvisiblering1{i},1)
            if (camvisiblering1{i}(k,1)== 1)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering1{i}(k,2)  i]));
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering1{i}(k,3:4)));
            elseif (camvisiblering1{i}(k,1)== 2)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering1{i}(k,2)+no_c1  i]));       
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering1{i}(k,3:4)));
            else
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering1{i}(k,2)+(no_c1+no_c2)  i]));  
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering1{i}(k,3:4)));
            end
        end
    end
    for i = 1:no_c2
        for k = 1:size(camvisiblering2{i},1)
            if (camvisiblering2{i}(k,1)== 1)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering2{i}(k,2)  i+no_c1]));
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering2{i}(k,3:4)));
            elseif (camvisiblering2{i}(k,1)== 2)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering2{i}(k,2)+no_c1  i+no_c1]));       
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering2{i}(k,3:4)));
            else
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblering2{i}(k,2)+(no_c1+no_c2)  i+no_c1]));  
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblering2{i}(k,3:4)));
            end
        end
    end
    for i = 1:no_c3
        for k = 1:size(camvisiblestatic{i},1)
            if (camvisiblestatic{i}(k,1)== 1)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblestatic{i}(k,2)  i+(no_c1+no_c2)]));
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblestatic{i}(k,3:4)));
            elseif (camvisiblestatic{i}(k,1)== 2)
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblestatic{i}(k,2)+no_c1  i+(no_c1+no_c2)]));       
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblestatic{i}(k,3:4)));
            else
                ret.cam_reproj.view = uint64(cat(1,ret.cam_reproj.view, [camvisiblestatic{i}(k,2)+(no_c1+no_c2)  i+(no_c1+no_c2)]));  
                ret.cam_reproj.pos = double(cat(1,ret.cam_reproj.pos, camvisiblestatic{i}(k,3:4)));
            end
        end
    end
    ret.cam_reproj.weight = double(cammutualvisibilityweight*ones(1,size(ret.cam_reproj.pos,1)));
end

% map from continuous index to cams,view 
camind = zeros(2,0); 
camind = cat(2, camind, [ones(1,size(ring1campos,2)); 1:size(ring1campos,2)]);
camind = cat(2, camind, [2*ones(1,size(ring2campos,2)); 1:size(ring2campos,2)]);
camind = cat(2, camind, [3*ones(1,size(staticcamrot,2)); 1:size(staticcamrot,2)]);

% map from cams,view to all visible points (i.e. index ret.points.pt3d{k}, ret.points.reproj_pos{k}(:,l)
camptidx = cell(1,length(ret.cams));
for i = 1:length(ret.cams)
    camptidx{i} = cell(1, ret.cams(i).noviews);
end
for k = 1:length(ret.points.reproj_view)
    for l = 1:length(ret.points.reproj_view{k})
        camid = camind(1,ret.points.reproj_view{k}(l));
        viewid = camind(2,ret.points.reproj_view{k}(l));
        camptidx{camid}{viewid} = cat(1,camptidx{camid}{viewid}, [k l]);
    end
end

%ret.BAopt.maxiter = uint64(1); 
%result = ceresBAmex(ret); 
 
% add noise to points, planes,  camera positions and rotations
if (~isempty(rndopt))
    % add noise to 2D observations, make noise dependent on 3D pt distance to camera
    if (sum(rndopt.randsigma_pt2d)>0)
        for k = 1:length(ret.points.reproj_view)
            for l = 1:length(ret.points.reproj_view{k})
                camid = camind(1,ret.points.reproj_view{k}(l));
                viewid = camind(2,ret.points.reproj_view{k}(l));
                if (size(ret.cams(camid).views_orien,2)>1)
                    R = func_rodrigues(ret.cams(camid).views_orien(:,viewid));
                else
                    R = func_rodrigues(ret.cams(camid).views_orien(:,1));        
                end
                if (size(ret.cams(camid).views_trans,2)>1)
                    t = ret.cams(camid).views_trans(:,viewid);
                else
                    t = ret.cams(camid).views_trans(:,1);        
                end

                pt3d = ret.points.pt3d(:,k);

                tmp = R*(pt3d-t);
                dis = norm(tmp) / mean(ret.cams(camid).fc); % distance from camera

                %func_render_view(ret, wh, camptidx, camid, viewid)
                %func_reproject_BAtest(pt3d, R,t,ret.cams(camid).fc,ret.cams(camid).cc,ret.cams(camid).kc,1)
                rndadd = mvnrnd([0,0], diag([rndopt.randsigma_pt2d(camid) rndopt.randsigma_pt2d(camid)]./dis))';
                %if (norm(rndadd) > tmptmp)
                %    tmptmp = norm(rndadd)
                %end
                %rndadd = mvnrnd([0,0], diag([rndopt.randsigma_pt2d(camid) rndopt.randsigma_pt2d(camid)]))';
                ret.points.reproj_pos{k}(:,l) = ret.points.reproj_pos{k}(:,l) + rndadd;
            end 
        end 
    end
    
    % add noise to 3D points
    ret.points.pt3d = ret.points.pt3d + normrnd(0,rndopt.randsigma_pt3d,size(ret.points.pt3d));
    
    ret.cams(1).views_trans = ret.cams(1).views_trans + normrnd(0,rndopt.randsigma_cam(1),size(ret.cams(1).views_trans));
    ret.cams(2).views_trans = ret.cams(2).views_trans + normrnd(0,rndopt.randsigma_cam(2),size(ret.cams(2).views_trans));
    ret.cams(3).views_trans = ret.cams(3).views_trans + normrnd(0,rndopt.randsigma_cam(3),size(ret.cams(3).views_trans));

    ret.cams(1).fc = ret.cams(1).fc + normrnd(0, rndopt.fc(1), [1 2]);
    ret.cams(2).fc = ret.cams(2).fc + normrnd(0, rndopt.fc(2), [1 2]);
    ret.cams(3).fc = ret.cams(3).fc + normrnd(0, rndopt.fc(3), [1 2]);

    ret.cams(1).cc = ret.cams(1).cc + normrnd(0, rndopt.cc(1), [1 2]);
    ret.cams(2).cc = ret.cams(2).cc + normrnd(0, rndopt.cc(2), [1 2]);
    ret.cams(3).cc = ret.cams(3).cc + normrnd(0, rndopt.cc(3), [1 2]);

    ret.cams(1).kc = ret.cams(1).kc + normrnd(0, abs(ret.cams(1).kc*rndopt.kc(1)), [1 2]);
    ret.cams(2).kc = ret.cams(2).kc + normrnd(0, abs(ret.cams(2).kc*rndopt.kc(2)), [1 2]);
    ret.cams(3).kc = ret.cams(3).kc + normrnd(0, abs(ret.cams(3).kc*rndopt.kc(3)), [1 2]);

    for k = 1:no_c1
        Rrand = func_get_rot_matrix((rand-.5) * pi/180*rndopt.maxrandangle_cam(1),(rand-.5) * pi/180*rndopt.maxrandangle_cam(1),(rand-.5) * pi/180*rndopt.maxrandangle_cam(1));
        R = ret.cams(1).views_orien(:,k);
        ret.cams(1).views_orien(:,k) = func_rodrigues(Rrand*func_rodrigues(R));
    end
    for k = 1:no_c2
        Rrand = func_get_rot_matrix((rand-.5) * pi/180*rndopt.maxrandangle_cam(2),(rand-.5) * pi/180*rndopt.maxrandangle_cam(2),(rand-.5) * pi/180*rndopt.maxrandangle_cam(2));
        R = ret.cams(2).views_orien(:,k);
        ret.cams(2).views_orien(:,k) = func_rodrigues(Rrand*func_rodrigues(R));
    end
    for k = 1:no_c3
        Rrand = func_get_rot_matrix((rand-.5) * pi/180*rndopt.maxrandangle_cam(3),(rand-.5) * pi/180*rndopt.maxrandangle_cam(3),(rand-.5) * pi/180*rndopt.maxrandangle_cam(3));
        R = ret.cams(3).views_orien(:,k);
        ret.cams(3).views_orien(:,k) = func_rodrigues(Rrand*func_rodrigues(R));
    end
    
    % Re-Compute planes for points
    if (opt(1)==1)
        for i = 1:6
            [coeff,score,latent] = pca(face{i});
            [~, mnidx] = min(latent);
            planevec{i} = coeff(:,mnidx); 
            planevec{i} = planevec{i} * (mean(face{i}) * planevec{i});
        end
        ret.planeconstraint.plane(:,1:6) = double(cat(2,planevec{:}));
    end
    
    % Re-Compute planes for cameras
    if (opt(2)==1)
        [coeff,score,latent] = pca(ret.cams(1).views_trans');
        [~, mnidx] = min(latent);
        ringvec1 = coeff(:,mnidx);
        ringvec1 = ringvec1 * (mean(ret.cams(1).views_trans') * ringvec1);

        [coeff,score,latent] = pca(ret.cams(2).views_trans');
        [~, mnidx] = min(latent);
        ringvec2 = coeff(:,mnidx);
        ringvec2 = ringvec2 * (mean(ret.cams(1).views_trans') * ringvec2);
        if (opt(1)==1)
            ret.planeconstraint.plane(:,7) = ringvec1;
            ret.planeconstraint.plane(:,8) = ringvec2;
        else
            ret.planeconstraint.plane(:,1) = ringvec1;
            ret.planeconstraint.plane(:,2) = ringvec2;
        end
    end    
end


maxaxis = [ min([min(ret.points.pt3d,[],2) min(ret.cams(1).views_trans,[],2) min(ret.cams(2).views_trans,[],2) min(ret.cams(3).views_trans,[],2)],[],2) ...
max([max(ret.cams(3).views_trans,[],2) max(ret.points.pt3d,[],2) max(ret.cams(1).views_trans,[],2) max(ret.cams(2).views_trans,[],2)],[],2)  ];






