
function depth = func_ptdepth(M,P)


m = P*M;
H = M(4,:); 
h = m(3,:); 

depth = sign(det(  P(1:3,1:3)  ))/norm( P(3,1:3)  ,2)*h./H;

depth(find(isnan(depth))) = 0;          