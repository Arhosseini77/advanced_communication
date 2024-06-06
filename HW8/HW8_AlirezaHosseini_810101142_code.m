clc; clear; close all;

% Define Parameters
SNR_db = 0 : 20;                       % SNR in dB
SNR_lin = 10.^(SNR_db/10);             % Linear SNR
N = 1e5;                               % N
x = [0.0281 0.3437 0.75 0.3437 0.0281]'; % x(t)
xx = conv(x, x);                        % Auto-correlation of x(t)

% Theoretical BER for BPSK in AWGN
theoretical_ber = qfunc(sqrt(2 * SNR_lin));

% Preallocate BER and SINR arrays
ber_zf = zeros(3, length(SNR_lin));
sinr_zf = zeros(3, length(SNR_lin));
ber_mmse = zeros(3, length(SNR_lin));
sinr_mmse = zeros(3, length(SNR_lin));
ber_dfe = zeros(3, length(SNR_lin));
sinr_dfe = zeros(3, length(SNR_lin));

% Define Equalizer Lengths
K_values = [2, 4, 6];
K2_values = [3, 5, 7];

% Calculate Viterbi BER
error_probabilities_bpsk = simulate_viterbi_bpsk(SNR_db, 1e6, 6);

% Run Simulations for ZF, MMSE, and DFE Equalizers
for snr_idx = 1:length(SNR_lin)
    for k_idx = 1:length(K_values)
        % ZF Equalizer
        zf_eq = zf_equalizer(K_values(k_idx), x);
        [ber_zf(k_idx, snr_idx), sinr_zf(k_idx, snr_idx)] = simulate_ber(SNR_lin(snr_idx), N, x, zf_eq);

        % MMSE Equalizer
        mmse_eq = mmse_equalizer(xx + [zeros(floor(length(x)/2), 1); x; zeros(floor(length(x)/2), 1)]/SNR_lin(snr_idx), x, K_values(k_idx));
        [ber_mmse(k_idx, snr_idx), sinr_mmse(k_idx, snr_idx)] = simulate_ber(SNR_lin(snr_idx), N, x, mmse_eq);

        % DFE Equalizer
        dfe_eq = dfe_equalizer(xx + [zeros(floor(length(x)/2), 1); x; zeros(floor(length(x)/2), 1)]/SNR_lin(snr_idx), x, K2_values(k_idx));
        [ber_dfe(k_idx, snr_idx), sinr_dfe(k_idx, snr_idx)] = simulate_ber_feedback(SNR_lin(snr_idx), N, x, dfe_eq);
    end
end

% Plot Results for ZF Equalizer
figure;
subplot(1, 2, 1);
hold on; grid on;
semilogy(SNR_db, ber_zf(1,:), '-o', 'DisplayName', 'ZF 5');
semilogy(SNR_db, ber_zf(2,:), '-x', 'DisplayName', 'ZF 9');
semilogy(SNR_db, ber_zf(3,:), '-s', 'DisplayName', 'ZF 13');
semilogy(SNR_db, theoretical_ber, '-.', 'DisplayName', 'Theoretical BER');
semilogy(SNR_db, error_probabilities_bpsk, 'o-', 'DisplayName', 'Viterbi BPSK');
xlabel('SNR (dB)');
ylabel('BER');
legend('show');
title('BER vs. SNR for ZF Equalizer');

subplot(1, 2, 2);
hold on; grid on;
semilogy(SNR_db, -sinr_zf(1,:), '-o', 'DisplayName', 'ZF 5');
semilogy(SNR_db, -sinr_zf(2,:), '-x', 'DisplayName', 'ZF 9');
semilogy(SNR_db, -sinr_zf(3,:), '-s', 'DisplayName', 'ZF 13');
xlabel('SNR (dB)');
ylabel('SINR');
legend('show');
title('SINR vs. SNR for ZF Equalizer');

% Plot Results for MMSE Equalizer
figure;
subplot(1, 2, 1);
hold on; grid on;
semilogy(SNR_db, ber_mmse(1,:), '--o', 'DisplayName', 'MMSE 5');
semilogy(SNR_db, ber_mmse(2,:), '--x', 'DisplayName', 'MMSE 9');
semilogy(SNR_db, ber_mmse(3,:), '--s', 'DisplayName', 'MMSE 13');
semilogy(SNR_db, theoretical_ber, '-.', 'DisplayName', 'Theoretical BER');
semilogy(SNR_db, error_probabilities_bpsk, 'o-', 'DisplayName', 'Viterbi BPSK');
xlabel('SNR (dB)');
ylabel('BER');
legend('show');
title('BER vs. SNR for MMSE Equalizer');

subplot(1, 2, 2);
hold on; grid on;
semilogy(SNR_db, sinr_mmse(1,:), '--o', 'DisplayName', 'MMSE 5');
semilogy(SNR_db, sinr_mmse(2,:), '--x', 'DisplayName', 'MMSE 9');
semilogy(SNR_db, sinr_mmse(3,:), '--s', 'DisplayName', 'MMSE 13');
xlabel('SNR (dB)');
ylabel('SINR');
legend('show');
title('SINR vs. SNR for MMSE Equalizer');

% Plot Results for DFE Equalizer
figure;
subplot(1, 2, 1);
hold on; grid on;
semilogy(SNR_db, ber_dfe(1,:), ':o', 'DisplayName', 'DFE 5');
semilogy(SNR_db, ber_dfe(2,:), ':x', 'DisplayName', 'DFE 9');
semilogy(SNR_db, ber_dfe(3,:), ':s', 'DisplayName', 'DFE 13');
semilogy(SNR_db, theoretical_ber, '-.', 'DisplayName', 'Theoretical BER');
semilogy(SNR_db, error_probabilities_bpsk, 'o-', 'DisplayName', 'Viterbi BPSK');
xlabel('SNR (dB)');
ylabel('BER');
legend('show');
title('BER vs. SNR for DFE Equalizer');

subplot(1, 2, 2);
hold on; grid on;
semilogy(SNR_db, sinr_dfe(1,:), ':o', 'DisplayName', 'DFE 5');
semilogy(SNR_db, sinr_dfe(2,:), ':x', 'DisplayName', 'DFE 9');
semilogy(SNR_db, sinr_dfe(3,:), ':s', 'DisplayName', 'DFE 13');
xlabel('SNR (dB)');
ylabel('SINR');
legend('show');
title('SINR vs. SNR for DFE Equalizer');

% Combined BER Plot for All Equalizers
figure;
hold on; grid on;
semilogy(SNR_db, ber_zf(2,:), '-x', 'DisplayName', 'ZF 9');
semilogy(SNR_db, ber_mmse(2,:), '--x', 'DisplayName', 'MMSE 9');
semilogy(SNR_db, ber_dfe(2,:), ':x', 'DisplayName', 'DFE 9');
semilogy(SNR_db, theoretical_ber, '-.', 'DisplayName', 'Theoretical BER');
semilogy(SNR_db, error_probabilities_bpsk, 'o-', 'DisplayName', 'Viterbi BPSK');
xlabel('SNR (dB)');
ylabel('BER');
legend('show');
title('BER vs. SNR for All Equalizers');

% Combined SINR Plot for All Equalizers
% figure;
% hold on; grid on;
% semilogy(SNR_db, sinr_zf(2,:), '-x', 'DisplayName', 'ZF 9');
% semilogy(SNR_db, sinr_mmse(2,:), '--x', 'DisplayName', 'MMSE 9');
% semilogy(SNR_db, sinr_dfe(2,:), ':x', 'DisplayName', 'DFE 9');
% xlabel('SNR (dB)');
% ylabel('SINR');
% legend('show');
% title('SINR vs. SNR for All Equalizers');

% Functions

% Simulate BER and SINR for given equalizer
function [ber, sinr] = simulate_ber(snr, pkt_size, x, d)
    % Generate random BPSK symbols
    I = randi([0 1], pkt_size, 1); 
    I(I == 0) = -1;
    
    % Transmit symbols through channel
    u = conv(I, x, 'same');
    noise = sqrt(0.5/snr) * (randn(size(u)) + 1j*randn(size(u)));
    y = u + noise;
    
    % Equalize the received signal
    I_hat = conv(y, d, 'same');
    
    % Decision device
    I_1 = abs(I_hat - 1);
    I_0 = abs(I_hat + 1);
    [~, idxs] = min([I_0, I_1], [], 2);
    I_hat(idxs == 1) = -1;
    I_hat(idxs == 2) = 1;
    
    % Calculate BER
    ber = sum(I_hat ~= I) / pkt_size;
    
    % Calculate SINR
    q = conv(x, d, 'same'); 
    q0 = q(floor(length(q)/2)+1);
    sinr = q0 / (1 - q0);
end

% Simulate BER and SINR for given DFE equalizer with feedback
function [ber, sinr] = simulate_ber_feedback(snr, pkt_size, x, d)
    Len = length(d);
    
    % Generate random BPSK symbols
    I = randi([0 1], pkt_size, 1); 
    I(I == 0) = -1;
    
    % Transmit symbols through channel
    u = conv(I, x);
    noise = sqrt(0.5/snr) * (randn(numel(u), 1) + 1j*randn(numel(u), 1));
    y = u + noise;
    
    % Equalize the received signal
    I_hat = conv(y, [d(1:floor(Len/2)); zeros(floor(Len/2), 1)]) + ...
            conv(conv(I, [0 0 1 0 0]), [zeros(floor(Len/2), 1); d(floor(Len/2)+1:end)]);
    I_hat = I_hat(floor((length(I_hat)-pkt_size)/2)+1:end-floor((length(I_hat)-pkt_size)/2)-1);
    
    % Decision device
    I_1 = abs(I_hat - 1);
    I_0 = abs(I_hat + 1);
    dists = [I_0, I_1];
    [~, idxs] = min(dists, [], 2);
    I_hat(idxs == 1) = -1;
    I_hat(idxs == 2) = 1;
    
    % Calculate BER
    ber = sum(I_hat ~= I) / pkt_size;
    
    % Calculate SINR
    q = conv(x, [d(1:floor(Len/2)); zeros(floor(Len/2), 1)]); 
    q0 = q(floor(length(q)/2));
    sinr = q0 / (1 - q0);
end

% ZF Equalizer design
function d = zf_equalizer(K, x)
    A = zeros(2*K+1);
    for i = 1:2*K+1
        A(i,i) = x(3);
        if i+1 <= 2*K+1, A(i,i+1) = x(4); end
        if i+2 <= 2*K+1, A(i,i+2) = x(5); end
        if i-1 >= 1, A(i,i-1) = x(2); end
        if i-2 >= 1, A(i,i-2) = x(1); end
    end
    b = zeros(2*K+1, 1);
    b(K+1) = 1;
    d = pinv(A) * b;
end

% MMSE Equalizer design
function d = mmse_equalizer(R_y, R_IY, K)
    A = zeros(2*K+1);
    for i = 1:2*K+1
        A(i,i) = R_y(5);
        if i+1 <= 2*K+1, A(i,i+1) = R_y(4); end
        if i+2 <= 2*K+1, A(i,i+2) = R_y(3); end
        if i+3 <= 2*K+1, A(i,i+3) = R_y(2); end
        if i+4 <= 2*K+1, A(i,i+4) = R_y(1); end
        if i-1 >= 1, A(i,i-1) = R_y(6); end
        if i-2 >= 1, A(i,i-2) = R_y(7); end
        if i-3 >= 1, A(i,i-3) = R_y(8); end
        if i-4 >= 1, A(i,i-4) = R_y(9); end
    end
    b = zeros(2*K+1, 1);
    b(K+1-floor(length(R_IY)/2):K+1+floor(length(R_IY)/2)) = R_IY;
    d = pinv(A) * b;
end

% DFE Equalizer design
function d = dfe_equalizer(R_y, R_IY, K2)
    A1 = zeros(K2, K2);
    A2 = zeros(K2, K2);
    A3 = diag(ones(K2, 1));
    for i = 1:K2
        A1(i,i) = R_y(5);
        if i+1 <= K2, A1(i,i+1) = R_y(4); end
        if i+2 <= K2, A1(i,i+2) = R_y(3); end
        if i+3 <= K2, A1(i,i+3) = R_y(2); end
        if i+4 <= K2, A1(i,i+4) = R_y(1); end
        if i-1 >= 1, A1(i,i-1) = R_y(6); end
        if i-2 >= 1, A1(i,i-2) = R_y(7); end
        if i-3 >= 1, A1(i,i-3) = R_y(8); end
        if i-4 >= 1, A1(i,i-4) = R_y(9); end
    end
    A2(K2, 1) = R_IY(2);
    A2(K2-1, 1) = R_IY(1);
    A2(K2, 2) = R_IY(1);
    b = zeros(2*K2, 1);
    b(K2) = R_IY(3); b(K2-1) = R_IY(2); b(K2-2) = R_IY(1);
    d = pinv([A1 A2; A2' A3]) * b;
end

% Viterbi Functions
function error_probabilities = simulate_viterbi_bpsk(SNR_dB_range, N, L)
    error_probabilities = zeros(size(SNR_dB_range));
    for i = 1:length(SNR_dB_range)
        SNR_dB = SNR_dB_range(i);
        total_errors = 0;
        total_bits = 0;
        
        while total_bits < N
            I_n = generate_bpsk_sequence(N);
            y_n = add_noise_bpsk(I_n, SNR_dB);
            decoded = viterbi_algorithm_bpsk(y_n, L);
            errors = sum(decoded ~= I_n);
            total_errors = total_errors + errors;
            total_bits = total_bits + length(I_n);
        end
        
        ber = total_errors / total_bits;
        error_probabilities(i) = ber;
%         fprintf('BPSK SNR: %d dB, BER: %f\n', SNR_dB, ber);
    end
end

function I_n = generate_bpsk_sequence(N)
    % Generate random BPSK symbols
    I_n = randsrc(1, N, [1, -1]);
end

function y_n = add_noise_bpsk(I_n, SNR_dB)
    SNR_linear = 10^(SNR_dB / 10);
    sigma_n = sqrt(1 / SNR_linear);
    noise = sigma_n * randn(size(I_n));
    y_n = I_n + noise;
end

function decoded = viterbi_algorithm_bpsk(y_n, L)
    N = length(y_n);
    trellis = inf(2, N+1);
    trellis(:, 1) = 0; % Starting state with zero path metric
    path = zeros(2, N);
    
    states = [1, -1];
    for i = 2:N+1
        for curr_state = 1:2
            for prev_state = 1:2
                path_metric = trellis(prev_state, i-1) + (y_n(i-1) - states(curr_state))^2;
                if path_metric < trellis(curr_state, i)
                    trellis(curr_state, i) = path_metric;
                    path(curr_state, i-1) = prev_state;
                end
            end
        end
    end
    
    decoded = zeros(1, N);
    [~, state] = min(trellis(:, N+1));
    for i = N:-1:1
        decoded(i) = states(state);
        state = path(state, i);
    end
end

