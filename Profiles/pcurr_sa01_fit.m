function pcurr_out = pcurr_sa01_fit(ac_vmec, xx)
pcurr_out = pcurr(xx, [0 1 ac_vmec(1:end)], [], [], 'sum_atan');
