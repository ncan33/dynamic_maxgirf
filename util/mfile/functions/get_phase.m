function phase_mod = get_phase(phase_mod)

nSMS = max(phase_mod(:))+1;
N = 0:nSMS-1;
N = permute(N,[3,4,5,1,2]);
phase_mod= exp(1i*phase_mod*2*pi/nSMS.*N);

end
