function func_plot_dataset(ret, opt, wh, ccent)

%opt(1): plot points
%opt(2): plot point planes
%opt(3): plot cameras
%opt(4): plot camera planes
%opt(5): plot camera mutual visibility

% image sizes


allpt = ret.points.pt3d';
if (opt(1)==1)
    scatter3(allpt(:,1), allpt(:,2), allpt(:,3), 15, 'filled', 'b');  hold on;
end
axis equal

% Compute and plot planes
if (opt(2)==1)

    for i = 1:6
        ptface{i} = find(arrayfun(@(x) ismember(i,ret.points.OnPlane{x}), 1:length(ret.points.OnPlane)));
    end
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

    for i = 1:6
        planevec{i} = ret.planeconstraint.plane(:,i);
        
        mf = mean(allpt(ptface{i},:));
        dmfp = mf*planevec{i}./norm(planevec{i}) - norm(planevec{i}); % distance of mean of pts on face from plane
        mf = mf - (dmfp*planevec{i}./norm(planevec{i}))';

        v = planevec{i}
        
        
        sum(  mf'/norm(mf) .* (v/norm(v)))
        
        
        
        
        diffe = sum(bsxfun(@times, allpt(facecorners{i},:), (planevec{i}./norm(planevec{i}))'),2) - norm(planevec{i}); % distance of corner points on plane
        cornersprojected = allpt(facecorners{i},:) - bsxfun(@times, diffe, (planevec{i}./norm(planevec{i}))') ; % project corners onto plane
        cornersprojected = cornersprojected + bsxfun(@minus, cornersprojected, mf) / 2; % make boundary .25 * cubesize larger for visualization purpose

        pt1 = cornersprojected(1,:);
        pt2 = cornersprojected(2,:);
        pt3 = cornersprojected(3,:);
        pt4 = cornersprojected(4,:);

        clear fv 
        fv.vertices =  [pt1; pt2; pt3; pt4];
        fv.faces(1,:) = [1 2 4 3]';
        pa = patch(fv,'FaceColor','red');
        alpha(pa,.1);
        axis equal;  
    
    end
end 


if (opt(3)==1)
    % plot cameras
    val = mean([sqrt(sum(diff(ret.cams(1).views_trans')'.^2)) sqrt(sum(diff(ret.cams(2).views_trans')'.^2))]);
    camerascaling = val/3;

    for i = 1:ret.cams(1).noviews
        K = eye(3);
        K(1,1) = ret.cams(1).fc(1); K(2,2) = ret.cams(1).fc(2); 
        K(1,3) = ret.cams(1).cc(1); K(2,3) = ret.cams(1).cc(2);
        R = func_rodrigues(ret.cams(1).views_orien(:,i));
        t = ret.cams(1).views_trans(:,i);
        func_plot_cameras(K,R, (-R*t)', wh(1), wh(2), [1 0 0], camerascaling, 1, 1);
    end
    for i = 1:ret.cams(2).noviews
        K = eye(3);
        K(1,1) = ret.cams(2).fc(1); K(2,2) = ret.cams(2).fc(2); 
        K(1,3) = ret.cams(2).cc(1); K(2,3) = ret.cams(2).cc(2);
        R = func_rodrigues(ret.cams(2).views_orien(:,i));
        t = ret.cams(2).views_trans(:,i);
        func_plot_cameras(K,R, (-R*t)', wh(1), wh(2), [0 0 1], camerascaling, 1, 1);
    end    
    for i = 1:ret.cams(3).noviews
        K = eye(3);
        K(1,1) = ret.cams(3).fc(1); K(2,2) = ret.cams(3).fc(2); 
        K(1,3) = ret.cams(3).cc(1); K(2,3) = ret.cams(3).cc(2);
        R = func_rodrigues(ret.cams(3).views_orien(:,i));
        t = ret.cams(3).views_trans;
        func_plot_cameras(K,R, (-R*t)', wh(1), wh(2), [0 0 0], camerascaling, 1, 1);
    end    
    axis equal
end


% plot mutual visibility
if (opt(5)==1)
    for i = 1:size(ret.cam_reproj.view,1)
        k1 = ret.cam_reproj.view(i,1);
        k2 = ret.cam_reproj.view(i,2);
        
        if (k1>=1 && k1<=ret.cams(1).noviews)
            t1 = ret.cams(1).views_trans(:,k1);
        elseif (k1>ret.cams(1).noviews && k1<=(ret.cams(1).noviews+ret.cams(2).noviews))
            t1 = ret.cams(2).views_trans(:,k1-ret.cams(1).noviews);
        elseif (k1>(ret.cams(1).noviews+ret.cams(2).noviews) && k1<=(ret.cams(1).noviews+ret.cams(2).noviews+ret.cams(3).noviews))
            t1 = ret.cams(3).views_trans;
        end
         
        if (k2>=1 && k2<=ret.cams(1).noviews)
            t2 = ret.cams(1).views_trans(:,k2);
        elseif (k2>ret.cams(1).noviews && k2<=(ret.cams(1).noviews+ret.cams(2).noviews))
            t2 = ret.cams(2).views_trans(:,k2-ret.cams(1).noviews);
        elseif (k2>(ret.cams(1).noviews+ret.cams(2).noviews) && k2<=(ret.cams(1).noviews+ret.cams(2).noviews+ret.cams(3).noviews))
            t2 = ret.cams(3).views_trans;
        end
        plot3([t1(1) t2(1)], [t1(2) t2(2)], [t1(3) t2(3)], '-k');
    end
end


% plot camera planes
if (opt(4)==1) 
    ring1campos = ret.cams(1).views_trans;
    if (opt(2)==1)  % are point planes used?
        ringvec1 = ret.planeconstraint.plane(:,7)';
    else
        ringvec1 = ret.planeconstraint.plane(:,1)';
    end
    ring2campos = ret.cams(2).views_trans;
    if (opt(2)==1)  % are point planes used?
        ringvec2 = ret.planeconstraint.plane(:,8)';
    else
        ringvec2 = ret.planeconstraint.plane(:,2)' ;
    end
    
    diffe = sum(bsxfun(@times, ring1campos, (ringvec1./norm(ringvec1))'),1) - norm(ringvec1); % distance of all camera points on plane
    correction = -bsxfun(@times, diffe, (ringvec1./norm(ringvec1))') ;
    planepts = ring1campos + correction; % project corners onto plane
    planepts = planepts + bsxfun(@minus, planepts, ccent') / 3;
    clear fv
    viewno = double(ret.cams(1).noviews);
    planepts = planepts(:,[floor(linspace(1,floor(viewno/2),4)) floor(linspace(floor(viewno/2)+1,viewno,4))]);
    fv.vertices =  planepts';
    %fv.faces(1,:) = floor(linspace(1,floor(viewno/2),4))';
    %fv.faces(2,:) = [floor(viewno/2) floor(viewno/2)+1 viewno 1 ]';
    %fv.faces(3,:) = floor(linspace(floor(viewno/2)+1,viewno,4))';
    fv.faces(1,:) = [1 2 3 4]';
    fv.faces(2,:) = [4 5 8 1]';
    fv.faces(3,:) = [5 6 7 8]';
    pa = patch(fv,'FaceColor','red', 'EdgeColor', 'none');
    alpha(pa,.1);
    plot3([ring1campos(1,:); ring1campos(1,:)+correction(1,:)], [ring1campos(2,:); ring1campos(2,:)+correction(2,:)], [ring1campos(3,:); ring1campos(3,:)+correction(3,:)], '-r');
    scatter3(ring1campos(1,:)+correction(1,:), ring1campos(2,:)+correction(2,:), ring1campos(3,:)+correction(3,:), 20, 'r');

    planepts = cat(2,planepts, planepts(:,1));
    plot3(planepts(1,:), planepts(2,:), planepts(3,:), '-k'); 
    
    diffe = sum(bsxfun(@times, ring2campos, (ringvec2./norm(ringvec2))'),1) - norm(ringvec2); % distance of all camera points on plane
    correction = -bsxfun(@times, diffe, (ringvec2./norm(ringvec2))') ;
    planepts = ring2campos + correction; % project corners onto plane
    planepts = planepts + bsxfun(@minus, planepts, ccent') / 6;
    clear fv
    viewno = double(ret.cams(2).noviews);
    planepts = planepts(:,[floor(linspace(1,floor(viewno/2),4)) floor(linspace(floor(viewno/2)+1,viewno,4))]);
    fv.vertices =  planepts';
    fv.faces(1,:) = [1 2 3 4]';
    fv.faces(2,:) = [4 5 8 1]';
    fv.faces(3,:) = [5 6 7 8]';
    pa = patch(fv,'FaceColor','blue', 'EdgeColor', 'none');
    alpha(pa,.1);
    plot3([ring2campos(1,:); ring2campos(1,:)+correction(1,:)], [ring2campos(2,:); ring2campos(2,:)+correction(2,:)], [ring2campos(3,:); ring2campos(3,:)+correction(3,:)], '-b');
    scatter3(ring2campos(1,:)+correction(1,:), ring2campos(2,:)+correction(2,:), ring2campos(3,:)+correction(3,:), 20, 'b');

    planepts = cat(2,planepts, planepts(:,1));
    plot3(planepts(1,:), planepts(2,:), planepts(3,:), '-k'); 

    axis equal
end

