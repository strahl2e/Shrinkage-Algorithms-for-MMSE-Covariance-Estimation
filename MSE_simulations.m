%% Replicate the MSE simulations from the paper

%% Setup AR1 data matrices
p=100;

S_AR1_powers = zeros(p);
for i=1:p
    S_AR1_powers(i, i:p) = abs(i - (i:p));
    S_AR1_powers(i:p, i) = S_AR1_powers(i, i:p);
end
figure
imagesc(S_AR1_powers)
title('S_AR1_powers')

S_AR1_r_0_1 = 0.1.^S_AR1_powers
figure
imagesc(S_AR1_r_0_1)
title('S_AR1_r_0_1')
S_AR1_r_0_5 = 0.5.^S_AR1_powers
figure
imagesc(S_AR1_r_0_5)
title('S_AR1_r_0_5')
S_AR1_r_0_9 = 0.9.^S_AR1_powers
figure
imagesc(S_AR1_r_0_9)
title('S_AR1_r_0_9')

