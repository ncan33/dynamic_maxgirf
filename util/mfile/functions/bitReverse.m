function brOrder = bitReverse(interleaves)

z = 0:interleaves - 1;
min2power = 2^ceil(log2(length(z)));
brOrder = bitrevorder(1:min2power);
brOrder = brOrder(brOrder<=length(z))-1;