%
% 2008-09-22 Martin Uecker <muecker@gwdg.de>
% Biomedizinische NMR Forschungs GmbH am
% Max-Planck-Institut fuer biophysikalische Chemie
%


function W = NLINV_weights(x,y)

    [X,Y] = meshgrid(1:x,1:y);
    d = ((X-1)/x-0.5).^2+((Y-1)/y-0.5).^2;
    W = 1./(1+220*d).^16;
	
end

