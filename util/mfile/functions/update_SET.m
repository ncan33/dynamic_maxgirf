function SET = update_SET
load('/Users/ytian/Documents/MATLAB/mfile/functions/SET.mat','SET')
SET = SET+1;
save('/Users/ytian/Documents/MATLAB/mfile/functions/SET.mat','SET')
