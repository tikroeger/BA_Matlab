function func_render_view(ret, wh, camptidx, camid, viewid)


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


    % get 3d points
    pt3d = ret.points.pt3d(:,camptidx{camid}{viewid}(:,1))';
            
    
    % get 2d points
    pt2d_reproj = func_reproject_BAtest(pt3d, R,t,ret.cams(camid).fc,ret.cams(camid).cc,ret.cams(camid).kc,1);
    pt2d_reproj = bsxfun(@plus, pt2d_reproj, wh'/2); % center points on [width height]/2 
    
    % reproject 3d -> 2D
    pt2d_measured = arrayfun(@(x)  ret.points.reproj_pos{camptidx{camid}{viewid}(x,1)}(:,camptidx{camid}{viewid}(x,2))  ,1:size(camptidx{camid}{viewid},1), 'UniformOutput',0);
    pt2d_measured = [pt2d_measured{:}];
    pt2d_measured= bsxfun(@plus, pt2d_measured, wh'/2);
    
     
    if (~isempty(ret.cams(camid).vplines))
        tt=ret.cams(camid).vplines(viewid).l(:,2);
        lines = cat(1,tt{:});
        if (~isempty(lines)) 
            lines = bsxfun(@plus, lines, [wh/2 wh/2]); % center points on [width height]/2 
            
            [err1, err2, vplines_proj] = func_reproject_VP_BAtest(ret, camid, viewid);        
            vpdisp = cat(1,vplines_proj{:});
            vpdisp = bsxfun(@plus, vpdisp, [wh/2 wh/2]);
            errvecdisp = cat(1,err2{:});
            err1{1} ;
            err2{1};
        end
    end
   
     
    clear clf;  
    h1 = plot(pt2d_measured(1,:),pt2d_measured(2,:), '.b'); hold on
    h2 = plot(pt2d_reproj(1,:),pt2d_reproj(2,:), 'or'); 
    if (~isempty(lines))
        h3 = plot([lines(:,1) lines(:,3)]', [lines(:,2) lines(:,4)]', '-*g','linewidth',2);
        h4 = plot([vpdisp(:,1) vpdisp(:,3)]', [vpdisp(:,2) vpdisp(:,4)]', '-ok','linewidth',2);
        legend([h1 h2 h3(1) h4(1)],{'measured pt', 'reprojected pt', 'GT VP line', 'reprojected VP lines'});
        
        quiver(vpdisp(:,3),vpdisp(:,4), errvecdisp(:,1),errvecdisp(:,2),'-k', 'AutoScale', 'off')
    else
        legend([h1 h2], {'measured pt', 'reprojected pt'})
    end
    axis image
    axis([1 wh(1) 2 wh(2)])

    