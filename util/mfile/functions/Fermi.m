function r = Fermi(t,F,t0,k)

r = F.*(1./(exp((t-t0)*k)+1));
r = r.*(exp(-t0*k)+1);