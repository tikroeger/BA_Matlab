
function depth = func_ptdepth(pt3d,P)

pt2d = P*pt3d;
depth = sign(det(  P(1:3,1:3)  ))/norm( P(3,1:3)  ,2)* pt2d(3,:) ./   pt3d(4,:) ;
depth(find(isnan(depth))) = 0;          