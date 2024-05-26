% Initial sequence generation and printing
N_init = 1000;  % Length of the initial sequence
SNR_target_init = 10;  % Desired SNR (example value)

% Generate i.i.d sequence I_n with equal probability from {-1, 1}
I_n_init = randsrc(1, N_init, [1, -1]);

% Calculate noise variance
sigma_n2_init = 1 / SNR_target_init;

% Generate Gaussian noise with variance sigma_n2
w_n_init = sqrt(sigma_n2_init) * randn(1, N_init);

% Generate the sequence y_n
y_n_init = I_n_init + w_n_init;

% Print the first few values as a sample
fprintf('I_n (first 10 values):\n');
disp(I_n_init(1:10));
fprintf('w_n (first 10 values):\n');
disp(w_n_init(1:10));
fprintf('y_n (first 10 values):\n');
disp(y_n_init(1:10));

% Main script to run the simulations
% Parameters
N = 1e6; % Length of the sequence, increased to better estimate low BER
L = 6; % Depth of Viterbi algorithm
SNR_dB_range = 0:20; % SNR range from 0dB to 20dB

% Simulate and plot the BER curve for BPSK
error_probabilities_bpsk = simulate_viterbi_bpsk(SNR_dB_range, N, L);
figure;
semilogy(SNR_dB_range, error_probabilities_bpsk, 'o-', 'LineWidth', 2);
hold on;

% Simulate and plot the BER curve for QPSK
error_probabilities_qpsk = simulate_viterbi_qpsk(SNR_dB_range, N, L);
semilogy(SNR_dB_range, error_probabilities_qpsk, 's-', 'LineWidth', 2);

% Plot settings
title('Bit Error Rate vs SNR for BPSK and QPSK Viterbi Algorithm (6L depth)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('BPSK', 'QPSK');
grid on;
hold off;

% BPSK Functions
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
        fprintf('BPSK SNR: %d dB, BER: %f\n', SNR_dB, ber);
    end
end

% QPSK Functions
function I_n = generate_qpsk_sequence(N)
    % Generate random QPSK symbols
    real_part = randsrc(1, N, [1/sqrt(2), -1/sqrt(2)]);
    imag_part = randsrc(1, N, [1/sqrt(2)*1i, -1/sqrt(2)*1i]);
    I_n = real_part + imag_part;
end

function y_n = add_noise_qpsk(I_n, SNR_dB)
    SNR_linear = 10^(SNR_dB / 10);
    sigma_n = sqrt(1 / SNR_linear);
    noise = sigma_n * (randn(size(I_n)) + 1i * randn(size(I_n)));
    y_n = I_n + noise;
end

function decoded = viterbi_algorithm_qpsk(y_n, L)
    N = length(y_n);
    trellis = inf(4, N+1);
    trellis(:, 1) = 0; % Starting state with zero path metric
    path = zeros(4, N);
    
    states = [1/sqrt(2) + 1i/sqrt(2), 1/sqrt(2) - 1i/sqrt(2), ...
              -1/sqrt(2) + 1i/sqrt(2), -1/sqrt(2) - 1i/sqrt(2)];
          
    for i = 2:N+1
        for curr_state = 1:4
            for prev_state = 1:4
                path_metric = trellis(prev_state, i-1) + abs(y_n(i-1) - states(curr_state))^2;
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

function error_probabilities = simulate_viterbi_qpsk(SNR_dB_range, N, L)
    error_probabilities = zeros(size(SNR_dB_range));
    for i = 1:length(SNR_dB_range)
        SNR_dB = SNR_dB_range(i);
        total_errors = 0;
        total_bits = 0;
        
        while total_bits < N
            I_n = generate_qpsk_sequence(N);
            y_n = add_noise_qpsk(I_n, SNR_dB);
            decoded = viterbi_algorithm_qpsk(y_n, L);
            errors = sum(decoded ~= I_n);
            total_errors = total_errors + errors;
            total_bits = total_bits + length(I_n);
        end
        
        ber = total_errors / total_bits;
        error_probabilities(i) = ber;
        fprintf('QPSK SNR: %d dB, BER: %f\n', SNR_dB, ber);
    end
end
