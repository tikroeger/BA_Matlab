%--------------------------------------------------------------------------
% File     : func_readwrite_BA_problem(Matlab M-Function)
% Author   : Till Kroeger
% Created  : 26.09.2014 (in Matlab 8.3.0.532, R2014a)
% Usage: 
%        varargout = func_readwrite_BA_problem(problem, filename)
%
% Description :
%       Write a BA_Matlab problem struct to file 'filename' in the input format to BAdjustBin.
%       If 'problem' is empty a BA_Matlab problem struct is read from 'filename'.
%
%--------------------------------------------------------------------------


function varargout = func_readwrite_BA_problem(problem, filename)


if (~isempty(problem))  % write problem to filename
    fprintf('Write BA problem to file %s\n',filename);
    % write to file
    no3Dpoints = size(problem.points.pt3d,2);  
    dlmwrite(filename, [problem.formatversion],  'delimiter', ' ', 'precision', 20);
    dlmwrite(filename, [problem.BAopt.cores problem.BAopt.timeout_sec problem.BAopt.maxiter], '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, [problem.BAopt.PTLoss],    '-append', 'delimiter', ' ', 'precision', 20);
    dlmwrite(filename, [problem.BAopt.PlaneLoss], '-append', 'delimiter', ' ', 'precision', 20);
    dlmwrite(filename, [problem.BAopt.DerivLoss], '-append', 'delimiter', ' ', 'precision', 20);
    dlmwrite(filename, [problem.BAopt.VPLoss], '-append', 'delimiter', ' ', 'precision', 20);
    
    dlmwrite(filename, [no3Dpoints cellfun(@(x) size(x,2), problem.points.reproj_view)], '-append', 'delimiter',' ', 'precision', 20); % no 3d points, and no reprojections for each point
    dlmwrite(filename, [size(problem.cams,2) [problem.cams.noviews]], '-append', 'delimiter',' ', 'precision', 20);  % no cameras, and no views for each camera
    dlmwrite(filename, size(problem.cam_reproj.view,1), '-append', 'delimiter',' ', 'precision', 20);  % number of camera reprojections
    dlmwrite(filename, size(problem.planeconstraint.plane,2), '-append', 'delimiter',' ', 'precision', 20);  % number of camera reprojections
    dlmwrite(filename, problem.planeconstraint.plane(:)', '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, problem.planeconstraint.planew(:)', '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, problem.planeconstraint.fixPlane(:)', '-append', 'delimiter',' ', 'precision', 20);

    dlmwrite(filename, size(problem.VPconstraint.vp,2), '-append', 'delimiter',' ', 'precision', 20);  % number of camera reprojections
    dlmwrite(filename, problem.VPconstraint.vp(:)', '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, problem.VPconstraint.vpw(:)', '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, problem.VPconstraint.fixVP(:)', '-append', 'delimiter',' ', 'precision', 20);
    
    for k = 1:length(problem.cams)
        dlmwrite(filename, problem.cams(k).fc(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).cc(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).kc(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).fixInternals , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).views_orien(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).fixOrientation , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).views_trans(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).fixTranslation , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).viewids(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).OnPlane(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).camweight(:)' , '-append', 'delimiter',' ', 'precision', 20);
        dlmwrite(filename, problem.cams(k).smootherM(:)' , '-append', 'delimiter',' ', 'precision', 20);
        
        if (isempty(problem.cams(k).vplines)) 
            dlmwrite(filename, [' '] , '-append', 'delimiter',' ', 'precision', 20);
        else
            tmp_ou = cell(1, problem.cams(k).noviews);
            for fr = 1:problem.cams(k).noviews % write all lines for all vps in one line per frame
                novpthis = size(problem.cams(k).vplines(fr).l,1);
                
                tmp = cell(1,novpthis);
                for kk = 1:novpthis 
                   nolinesthis = double(size(problem.cams(k).vplines(fr).l{kk,2},1));
                   vpid = double(problem.cams(k).vplines(fr).l{kk,1});
                   tmp{kk} = [vpid nolinesthis problem.cams(k).vplines(fr).l{kk,2}(:)']; 
                end 
                tmp_ou{fr} = [double(novpthis) [tmp{:}]];
            end
            dlmwrite(filename, [tmp_ou{:}] , '-append', 'delimiter',' ', 'precision', 20); % save VP count, VP id, no lines, linearized line endpoints
        end
    end
    tmp = problem.cam_reproj.view'; tmp = tmp(:);
    dlmwrite(filename, tmp', '-append', 'delimiter',' ', 'precision', 20);  
    tmp = problem.cam_reproj.pos'; tmp = tmp(:);
    dlmwrite(filename, tmp', '-append', 'delimiter',' ', 'precision', 20);  
    dlmwrite(filename, problem.cam_reproj.weight(:)', '-append', 'delimiter',' ', 'precision', 20);  
    dlmwrite(filename, problem.points.pt3d(:)', '-append', 'delimiter',' ', 'precision', 20);
    dlmwrite(filename, problem.points.fixPosition, '-append', 'delimiter',' ', 'precision', 20);  
    dlmwrite(filename, problem.points.pointw, '-append', 'delimiter',' ', 'precision', 20);  
    dlmwrite(filename, sum(~cellfun(@isempty, problem.points.OnPlane)) , '-append', 'delimiter',' ', 'precision', 20);  
    for k = 1:length(problem.points.OnPlane)
        if (~isempty(problem.points.OnPlane{k}))
            dlmwrite(filename, [k problem.points.OnPlane{k}], '-append', 'delimiter',' ', 'precision', 20);  
        end
    end
    for k = 1:no3Dpoints
        dlmwrite(filename, problem.points.reproj_view{k}(:)', '-append', 'delimiter',' ', 'precision', 20);  
        dlmwrite(filename, problem.points.reproj_pos{k}(:)', '-append', 'delimiter',' ', 'precision', 20);          
    end
    
else
    fprintf('Read BA problem from file %s\n',filename);
    % read from file
    fid = fopen(filename);
    allData = textscan(fid,'%s','Delimiter','\n');
    fclose('all'); 

    clear problem;
    it = 1;
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.formatversion = tmp(1);
    tmp = uint64(str2num(allData{1}{it})); it = it + 1;
    problem.BAopt.cores = tmp(1);
    problem.BAopt.timeout_sec = tmp(2);
    problem.BAopt.maxiter = tmp(3);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.BAopt.PTLoss = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.BAopt.PlaneLoss = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.BAopt.DerivLoss = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.BAopt.VPLoss = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    no3Dpoints = tmp(1);
    noreproject = tmp(2:end);
    tmp = str2num(allData{1}{it}); it = it + 1;
    nocams = tmp(1);
    noviews = tmp(2:end); 
    tmp = str2num(allData{1}{it}); it = it + 1;
    nocam_reproj = tmp(1);    
    
    tmp = str2num(allData{1}{it}); it = it + 1;
    noplanes = tmp(1);    
    tmp = str2num(allData{1}{it}); it = it + 1;
    if (length(tmp)>0)
        problem.planeconstraint.plane = reshape(double(tmp),3,length(tmp)/3);
    else
        problem.planeconstraint.plane = double([]);
    end
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.planeconstraint.planew = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.planeconstraint.fixPlane = uint8(tmp);

    tmp = str2num(allData{1}{it}); it = it + 1;
    novps = tmp(1);    
    tmp = str2num(allData{1}{it}); it = it + 1;
    if (length(tmp)>0)
        problem.VPconstraint.vp = reshape(double(tmp),3,length(tmp)/3);
    else
        problem.VPconstraint.vp = double([]);
    end
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.VPconstraint.vpw = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.VPconstraint.fixVP = uint8(tmp);

    view_cnt = 1;
    for k = 1:nocams
        problem.cams(k).noviews = uint64(noviews(k));
        %problem.cams(k).viewids = uint64((1:problem.cams(k).noviews)+view_cnt-1);
        view_cnt = view_cnt+problem.cams(k).noviews;
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).fc = double(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).cc = double(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).kc = double(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).fixInternals = uint8(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).views_orien = reshape(double(tmp),3,length(tmp)/3);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).fixOrientation = uint8(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).views_trans = reshape(double(tmp),3,length(tmp)/3);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).fixTranslation = uint8(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).viewids= uint64(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).OnPlane = uint64(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.cams(k).camweight = double(tmp);        
        tmp = str2num(allData{1}{it}); it = it + 1;
        if (length(tmp)==problem.cams(k).noviews^2) % full smoother matrix
            problem.cams(k).smootherM = reshape(double(tmp),problem.cams(k).noviews,problem.cams(k).noviews);
        elseif (length(tmp)==problem.cams(k).noviews) % instead of full matrix, same smoothing vector for all cameras 
            problem.cams(k).smootherM = double(tmp)'; % make column vector
        else
            problem.cams(k).smootherM = double([]);
        end
        
        tmp = str2num(allData{1}{it}); it = it + 1;
        if (length(tmp)>problem.cams(k).noviews) % string with all vp ids and lines in all views
            cnt=1;
            for fr = 1:problem.cams(k).noviews % write all lines for all vps in one line per frame
                novpsthis = tmp(cnt);
                cnt = cnt+1;
                problem.cams(k).vplines(fr).l = cell(novpsthis,2);
                for vpcnt = 1:novpsthis
                    problem.cams(k).vplines(fr).l{vpcnt,1} = uint64(tmp(cnt));
                    cnt = cnt+1;
                    nolines = tmp(cnt);
                    cnt = cnt+1;
                    problem.cams(k).vplines(fr).l{vpcnt,2} = reshape(tmp(cnt:(cnt + 4*nolines - 1)), nolines, 4);
                    cnt = cnt + 4*nolines;
                end
            end
        else
            problem.cams(k).vplines = double([]);
        end
        
    end
    tmp = str2num(allData{1}{it}); it = it + 1;
    if (length(tmp)>0)
        problem.cam_reproj.view = reshape(uint64(tmp), 2,length(tmp)/2)';
    else
        problem.cam_reproj.view = uint64([]);
    end 
    tmp = str2num(allData{1}{it}); it = it + 1;
    if (length(tmp)>0)
        problem.cam_reproj.pos = reshape(double(tmp), 2,length(tmp)/2)';
    else
        problem.cam_reproj.pos = double([]);
    end
    tmp = str2num(allData{1}{it}); it = it + 1;    
    problem.cam_reproj.weight = double(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.points.pt3d = reshape(double(tmp), 3,length(tmp)/3);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.points.fixPosition = uint8(tmp);
    tmp = str2num(allData{1}{it}); it = it + 1;
    problem.points.pointw = double(tmp);

    problem.points.OnPlane = cell(1,size(problem.points.pt3d,2));
    tmp = str2num(allData{1}{it}); it = it + 1;
    nopointsonplanes = tmp(1);
    for k = 1:nopointsonplanes
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.points.OnPlane{tmp(1)} = uint64(tmp(2:end));    
    end
    %tmp = str2num(allData{1}{it}); it = it + 1;
    %problem.points.OnPlane = uint64(tmp);    
    
    problem.points.reproj_view = cell(1,no3Dpoints);
    problem.points.reproj_pos = cell(1,no3Dpoints);
    for k = 1:no3Dpoints
        tmp = str2num(allData{1}{it}); it = it + 1;
        problem.points.reproj_view{k} = uint64(tmp);
        tmp = str2num(allData{1}{it}); it = it + 1;        
        problem.points.reproj_pos{k} = reshape(double(tmp), 2,length(tmp)/2);
    end
    
    varargout{1} = problem;
    
    if (it <= size(allData{1},1)) % check if residual is written, if yes, return it.
        tmp = str2num(allData{1}{it}); it = it + 1;
        varargout{2} = double(tmp);
    end
    
end

%match test
%[match, er1, er2, erc]  = comp_struct(problem, problem,2);
