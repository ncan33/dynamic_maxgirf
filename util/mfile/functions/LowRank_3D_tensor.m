function tensor_LowRank = LowRank_3D_tensor(tensor0)

[s1,s2,s3] = size(tensor0);

tensor1 = reshape(tensor0,[s1,s2*s3]);
[U1,S1,V1] = svd(tensor1,0);

tensor2 = reshape(permute(tensor0,[2,1,3]),[s2,s1*s3]);
[U2,S2,V2] = svd(tensor2,0);

tensor3 = reshape(permute(tensor0,[3,1,2]),[s3,s1*s2]);
[U3,S3,V3] = svd(tensor3,0);

% R1 = S1*V1';
% R2 = S2*V2';
% R3 = S3*V3';

S = U1'*reshape(tensor0,[s1,s2*s3]);
S = reshape(S,[s1,s2,s3]);

S = permute(S,[2,1,3]);
S = reshape(S,s2,s1*s3);
S = U2'*S;
S = reshape(S,s2,s1,s3);
S = permute(S,[2,1,3]);

S = permute(S,[3,1,2]);
S = reshape(S,s3,s1*s2);
S = U3'*S;
S = reshape(S,s3,s1,s2);
S = permute(S,[2,3,1]);

S_abs = abs(S);
S(S_abs<0.002*S_abs(1)) = 0;


tensor_LowRank = U1*reshape(S,[s1,s2*s3]);
tensor_LowRank = reshape(tensor_LowRank,[s1,s2,s3]);

tensor_LowRank = permute(tensor_LowRank,[2,1,3]);
tensor_LowRank = reshape(tensor_LowRank,s2,s1*s3);
tensor_LowRank = U2*tensor_LowRank;
tensor_LowRank = reshape(tensor_LowRank,s2,s1,s3);
tensor_LowRank = permute(tensor_LowRank,[2,1,3]);

tensor_LowRank = permute(tensor_LowRank,[3,1,2]);
tensor_LowRank = reshape(tensor_LowRank,s3,s1*s2);
tensor_LowRank = U3*tensor_LowRank;
tensor_LowRank = reshape(tensor_LowRank,s3,s1,s2);
tensor_LowRank = permute(tensor_LowRank,[2,3,1]);
