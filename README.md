# Dynamic MaxGIRF
A repository in progress, which demonstrates the feasibility of [MaxGIRF](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29232) for dynamic MRI.

## STCR routine
Call `dual_te_STCR_parameter_sweep()` to create a parameter sweep. Within this function, `dual_te_STCR_wrapper()` is called to perform dual-echo STCR. This function splits kspace data into two echoes, then performs `STCR()` on each echo separately.

The `dual_te_STCR_wrapper()` function is a modified version of [Ye Tian](https://scholar.google.com/citations?user=ia7_fu8AAAAJ&hl=en&oi=ao)'s `dual_te_recon.m`, and `STCR()` is a modified version of [Yongwan Lim](https://scholar.google.com/citations?user=z_uvzWIAAAAJ&hl=en&oi=ao)'s `demo_recon_MATLAB.m`. Original versions of these functions can be found in the `./ye_and_yongwan_functions` directory.

## Note to self
Add details about the RTHawk sequence. Figures demonstrating how dual-echo is encoded and explaining the variable density trajectory will be helpful.
