function kSpace = GEReadP(path)

%% Use 'GERecon' to read k-space data from P-file
pfile = GERecon('Pfile.Load',path);
for p = 1:pfile.phases
    for s = 1:pfile.slices
        for e = 1:pfile.echoes
            for c = 1:pfile.channels
                kSpace(:,:,s,e,c) = GERecon('Pfile.KSpace', s, e, c);
            end
        end
    end
end