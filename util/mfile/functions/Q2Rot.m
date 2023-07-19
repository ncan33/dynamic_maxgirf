function rot = Q2Rot(w, x, y, z)

rot(1, 1) = 1 - 2 * y * y - 2 * z * z;
rot(2, 1) =     2 * x * y + 2 * w * z;
rot(3, 1) =     2 * x * z - 2 * w * y;

rot(1, 2) =     2 * x * y - 2 * w * z;
rot(2, 2) = 1 - 2 * x * x - 2 * z * z;
rot(3, 2) =     2 * y * z + 2 * w * x;

rot(1, 3) =     2 * x * z + 2 * w * y;
rot(2, 3) =     2 * y * z - 2 * w * x;
rot(3, 3) = 1 - 2 * x * x - 2 * y * y;