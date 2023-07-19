function out = sfth(in, param)

out = in;
d = in-param;
idx = abs(in)>param;
out(idx) = d(idx);
out(~idx) = 0;

end