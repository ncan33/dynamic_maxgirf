function [first_est,para] = get_initial_estimation_NUFFT(sens_map,kSpace_all,theta_all,phase_mod_all,RayPosition,kSpace,phase_mod,NUFFT,para)

img_FE = NUFFT'* bsxfun(@times, kSpace, conj(phase_mod));
img_FE = img_FE/mean(abs(img_FE(:)))*0.0215;
first_est = single(sum(bsxfun(@times, img_FE, conj(sens_map)),4));