function score = func_eval_residual(ret_test, ret_gt)




pterr = sum((sum((ret_test.points.pt3d - ret_gt.points.pt3d).^2)));

c1err = sum((sum((ret_test.cams(1).views_trans - ret_gt.cams(1).views_trans).^2)));
c2err = sum((sum((ret_test.cams(2).views_trans - ret_gt.cams(2).views_trans).^2)));
c3err = 8*sum((sum((ret_test.cams(3).views_trans - ret_gt.cams(3).views_trans).^2)));

score = pterr + c1err + c2err + c3err;


