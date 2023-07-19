function u_ = TGV_primal_dual(Data, para)

ifplot  = para.setting.ifplot;
sens    = Data.sens_map;
kSpace  = Data.kSpace;
N       = Data.N;
maxits  = para.Recon.noi;

% regularization weight
alpha0 = para.Recon.weight_tTV * 2;
alpha1 = para.Recon.weight_tTV;

% step size
sigma = 1/8 * para.Recon.step_size;
tau = 1/16 * para.Recon.step_size;

[sx, sy, nof] = size(Data.first_est);

p = zeros(sx, sy, nof, 'single');
u_ = Data.first_est;
v_ = zeros(sx, sy, nof, 'single');
q = zeros(sx, sy, nof, 'single');
r = zeros(size(kSpace), 'single');
u = u_;
v = v_;

clear Data para

for i = 1:maxits
    % p
    p = p - sigma * (-dtp(u_) + v_);
    
    denom = max(1, abs(p) / alpha1);
    p = p ./ denom;
    
    % q
    q = q - sigma * -dtm(v_);
    
    denom = max(1, abs(q) / alpha0);
    q = q ./ denom;
    
    % r
    kSpace_update = NUFFT.NUFFT(u_ .* sens, N);
    kSpace_update =  kSpace_update - kSpace;
    r = (r + sigma * kSpace_update) / (1 + sigma);
    
    % u
    u_old = u;
    ww = sum(NUFFT.NUFFT_adj(r, N) .* conj(sens), 4);
    u = u - tau * (-dtm(p) + ww);
    u_ = 2 * u - u_old;
    
    % v
    v_old = v;
    v = v - tau * (-p - dtp(q));
    v_ = 2 * v - v_old;
    
    if ifplot
        fidelity_norm(i+1) = sum(abs(vec(kSpace_update)).^2) / 126^2;
        TGV1 = alpha1 * sum(abs(vec(dtp(u_) - v)));
        TGV2 = alpha0 * sum(abs(vec(dtm(v))));
        TGV(i+1) = TGV1 + TGV2;
        
        figure(100)
        clf;
        subplot(2,2,[1 3])
        fprintf('TGV2-L2-2D-PD: it = %4d\n', i);
        imagesc(abs(u_(:,:,90)));
        axis image
        axis off
        colormap gray
        drawnow
        
        subplot(2,2,[2 4])
        loglog(fidelity_norm, 'LineWidth', 2)
        hold on
        plot(TGV, 'LineWidth', 2)
        plot(fidelity_norm + TGV, 'LineWidth', 2)
        drawnow
    end
end