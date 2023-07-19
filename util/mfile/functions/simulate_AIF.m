function AIF = simulate_AIF(time,M)

if nargin == 1
    M = 'REST';
end

switch M
    case 'REST'
        tao = 9.3;%rest
    case 'STRESS'
        tao = 4.6;%stress
end

AIF = time.^2.*exp(-time/tao);