function SAMS_prescription(vec1, vec2)

z0 = [0; 0; 1];

% YT: vector between 2 vectors
vecBetween12_1 = [0; 0; 0];
vecBetween12_1(1) = vec1(1) + vec2(1);
vecBetween12_1(2) = vec1(2) + vec2(2);
vecBetween12_1(3) = vec1(3) + vec2(3);
vecBetween12_1_Norm = norm(vecBetween12_1);

if vecBetween12_1_Norm == 0
    vecBetween12_1 = vec1;
else
    vecBetween12_1(1) = vecBetween12_1(1) / vecBetween12_1_Norm;
    vecBetween12_1(2) = vecBetween12_1(2) / vecBetween12_1_Norm;
    vecBetween12_1(3) = vecBetween12_1(3) / vecBetween12_1_Norm;
end

% YT: vector between 2 vectors, with one flipped
vecBetween12_2 = [0; 0; 0];
vecBetween12_2(1) = vec1(1) - vec2(1);
vecBetween12_2(2) = vec1(2) - vec2(2);
vecBetween12_2(3) = vec1(3) - vec2(3);
vecBetween12_2_Norm = norm(vecBetween12_2);

if vecBetween12_2_Norm == 0
    vecBetween12_2 = vec1;
else
    vecBetween12_2(1) = vecBetween12_2(1) / vecBetween12_2_Norm;
    vecBetween12_2(2) = vecBetween12_2(2) / vecBetween12_2_Norm;
    vecBetween12_2(3) = vecBetween12_2(3) / vecBetween12_2_Norm;
end

% YT: find whether vecBetween12_1 or vecBetween12_2 gives the smallest angle
dotProd1 = sum(vecBetween12_1 .* vec1);
dotProd2 = sum(vecBetween12_2 .* vec1);

if dotProd1 > dotProd2
    prescribedVector = vecBetween12_2;
else
    prescribedVector = vecBetween12_1;
end

% YT: rotation matrix that rotates [0, 0, 1] to prescribedVector
rot_ = zeros(3);
rot_(1, 1) = 1 - prescribedVector(1) * prescribedVector(1) / (1 + prescribedVector(3));
rot_(1, 2) =   - prescribedVector(1) * prescribedVector(2) / (1 + prescribedVector(3));
rot_(1, 3) =     prescribedVector(1);
rot_(2, 1) =   - prescribedVector(1) * prescribedVector(2) / (1 + prescribedVector(3));
rot_(2, 2) = 1 - prescribedVector(2) * prescribedVector(2) / (1 + prescribedVector(3));
rot_(2, 3) =     prescribedVector(2);
rot_(3, 1) =   - prescribedVector(1);
rot_(3, 2) =   - prescribedVector(2);
rot_(3, 3) =     prescribedVector(3);

% YT: convert rotation matrix to Pitch, Roll, Yaw
prescribedAngle(1) =   atan2( rot_(2, 1), rot_(1, 1)) * 180 / pi;
prescribedAngle(2) =   atan2(-rot_(3, 1), sqrt(rot_(3, 2) * rot_(3, 2) + rot_(3, 3) * rot_(3, 3))) * 180 / pi;
prescribedAngle(3) = - atan2( rot_(3, 2), rot_(3, 3))* 180 / pi;

% YT: inverse rotation matrix
rot_inv = rot_';

% YT: rot vector 1 back
vec1_ = rot_inv * vec1;

% YT: rot vector 2 back
vec2_ = rot_inv * vec2;

prescribedAngle
vec1
vec2
vecBetween12_1
vecBetween12_2
prescribedVector
vec1_
vec2_

figure
plot3([0, prescribedVector(1)], [0, prescribedVector(2)], [0, prescribedVector(3)], 'black')
hold on
plot3([0, vec1(1)], [0, vec1(2)], [0, vec1(3)], 'Color', colors(1))
plot3([0, vec2(1)], [0, vec2(2)], [0, vec2(3)], 'Color', colors(1))

plot3([0, 0], [0, 0], [0, 1], '--', 'Color', 'black')
plot3([0, vec1_(1)], [0, vec1_(2)], [0, vec1_(3)], '--', 'Color', colors(1))
plot3([0, vec2_(1)], [0, vec2_(2)], [0, vec2_(3)], '--', 'Color', colors(1))

xlabel 'x'
ylabel 'y'
zlabel 'z'

