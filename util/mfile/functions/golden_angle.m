function GA = golden_angle(varargin)

GA = ((sqrt(5)-1)/2)*pi;
switch length(varargin)
    case 0
        return
    case 1
        GA = mod(GA*(varargin{1}-1),pi);
        return
    case 2
        GA = mod(GA*(varargin{1}-1),varargin{2});
        return
    otherwise
        error('too many inputs')
        return
end
