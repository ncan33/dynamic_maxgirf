function [tTV_sweep, sTV_sweep] = generate_anchored_sweep_vectors(n_tTV_steps, ...
    n_sTV_steps, tTV_step_factor, sTV_step_factor, tTV_anchor, max_sTV)
    % Creates two sweep vectors: sTV sweep and tTV sweep, which has the
    % lambda values for each constraint.
    
    tTV_sweep = zeros(1,n_tTV_steps);
    
    for i = 1:length(tTV_sweep)
        k = -length(tTV_sweep) + 1 + i;
        tTV_sweep(i) = tTV_anchor * (tTV_step_factor^(k));
    end
    
    sTV_sweep = zeros(1,n_sTV_steps);
    for i = 1:n_sTV_steps
        sTV_sweep(i) = max_sTV*sTV_step_factor^(-i+1);
    end
    sTV_sweep = fliplr(sTV_sweep);
    
    tTV_sweep = [0, tTV_sweep]; % add zero column
    sTV_sweep = [0, sTV_sweep]; % add zero column
end