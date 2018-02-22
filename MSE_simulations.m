%% Replicate the MSE simulations from the paper

%% Setup general
p=100;

S_i_j_diffs = zeros(p);
for i=1:p
    S_i_j_diffs(i, i:p) = abs(i - (i:p));
    S_i_j_diffs(i:p, i) = S_i_j_diffs(i, i:p);
end
figure
imagesc(S_i_j_diffs)
title('S_i_j_diffs')


%% Setup AR1 covariance matrices

S_AR1_r_0_1 = 0.1.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_1)
title('S_AR1_r_0_1')
S_AR1_r_0_5 = 0.5.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_5)
title('S_AR1_r_0_5')
S_AR1_r_0_9 = 0.9.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_9)
title('S_AR1_r_0_9')

%% setup FBM covariance matrices

FBM_matrix = @(h, S_i_j_diffs) 0.5*((S_i_j_diffs + 1).^(2*h) - 2*(S_i_j_diffs.^(2*h)) + (S_i_j_diffs - 1).^(2*h));

S_FBM_h_0_6 = FBM_matrix(0.6, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_6)
title('S_FBM_h_0_6')

S_FBM_h_0_7 = FBM_matrix(0.7, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_7)
title('S_FBM_h_0_7')

S_FBM_h_0_8 = FBM_matrix(0.8, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_8)
title('S_FBM_h_0_8')


% delete(findall(0,'Type','figure'))