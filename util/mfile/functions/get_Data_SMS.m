function [Data, para] = get_Data_SMS(Data, para)
%--------------------------------------------------------------------------
%   [Data, para] = get_Data_SMS(Data, para)
%--------------------------------------------------------------------------
%   Create some required fields in Data and para structures
%--------------------------------------------------------------------------
%   Inputs:
%       - Data                      [structure]
%           Data.kSpace             [sx, sy, nof, nc, nset, 1, nSMS]
%       - para                      [structure]
%           para.Recon              [structure]
%               Recon.nSMS          [scalar]
%
%           'sx'    number of readout per ray
%           'sy'    same as 'sx' for radial case
%           'nof'   number of time frames
%           'nc'    number of coils
%           'nset'  number of SMS sets
%           'nSMS'  number of simultaneous multi-slice
%
%       - Data                      see 'help STCR_conjugate_gradient.m'
%           Data.kSpace             pre-interpolated SMS k-space data
%       - para                      see 'help STCR_conjugate_gradient.m'
%           para.Recon.nSMS         number of simultaneous multi-slice
%--------------------------------------------------------------------------
%   Output:
%       - Data                      [structure]
%           Data.kSpace             [sx, sy, nof, nc, nset, 1,    nSMS]
%           Data.mask               [sx, sy, nof, 1,  nset, 1,    nSMS]
%           Data.SMS                [1,  1,  1,   1,  1,    nSMS, nSMS]
%           Data.first_est          [sx, sy, nof, 1,  nset, nSMS]
%           Data.filter             [sx, sy]
%
%       - para                      [structure]
%           para.Recon              [structure]
%               Recon.kSpace_size   [1, 2]
%               Recon.image_size    [1, 2]
%
%       - Data                      see 'help STCR_conjugate_gradient.m'
%           Data.kSpace             processed after fftshift
%           Data.mask               undersample mask
%           Data.SMS                SMS forward/backword operator [1]
%           Data.first_est          initial estimation "A^H d"
%           Data.filter             density compensation function for 
%                                   pre-interpolation
%
%       - para                      see 'help STCR_conjugate_gradient.m'              
%           Recon.kSpace_size   	first two dimensions of k-space 
%           Recon.image_size        first two dimensions of image
%--------------------------------------------------------------------------
%   Reference:
%       [1] Feasibility of multiple-view myocardial perfusion MRI using 
%           radial simultaneous multi-slice acquisitions. PLoS One, 2019, 
%           14(2):e0211738.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

Data.kSpace = fftshift2(Data.kSpace);
Data.mask   = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));
Data.kSpace = ifft2(Data.kSpace);
Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = fft2(Data.kSpace);
Data.kSpace = Data.kSpace.*Data.mask;
nSMS = para.Recon.nSMS;
Data.SMS = exp(1i*(0:nSMS-1)*2*pi/nSMS.*(0:nSMS-1).');
Data.SMS = permute(Data.SMS,[3,4,5,6,1,7,2]);
kSpace_sms = sum(Data.kSpace.*conj(Data.SMS),7);
Data.first_est = ifft2(kSpace_sms);
if ~isfield(Data,'sens_map')
    Data.sens_map = get_sens_map(Data.first_est,'SMS');
end
para.Recon.kSpace_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.image_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.sx = size(Data.kSpace,1);
Data.filter = ramp_filter_for_pre_interp(para);
Data.first_est = ifft2(kSpace_sms.*Data.filter);
Data.first_est = Data.first_est.*Data.sens_map;
Data.first_est = sum(Data.first_est,4);
