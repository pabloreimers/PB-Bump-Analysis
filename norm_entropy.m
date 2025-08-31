function norm_entropy = norm_entropy(data, edges)
    % NORMALIZED_SHANNON_ENTROPY computes the normalized Shannon entropy of a 1D signal.
    %   norm_entropy = normalized_shannon_entropy(data, num_bins)
    %
    %   Input:
    %     data:      A 1D numerical array (e.g., a signal).
    %     num_bins:  The number of bins to use for the histogram.
    %
    %   Output:
    %     norm_entropy: The normalized Shannon entropy value (0 to 1).

    % 1. Estimate the probability distribution using a histogram
    % 'Normalization', 'probability' ensures the sum of histogram values is 1.
    [N, ~] = histcounts(data, edges, 'Normalization', 'probability');
    num_bins = length(N);

    % 2. Calculate the standard Shannon entropy
    % Formula: H(X) = -sum(p .* log2(p))
    % Avoid log2(0) by replacing zeros with a negligible value.
    p = N;
    p(p == 0) = []; % Remove zero probabilities
    shannon_entropy = -sum(p .* log2(p));

    % 3. Calculate the maximum possible entropy
    % H_max = log2(number of states)
    max_entropy = log2(num_bins);

    % 4. Normalize the entropy
    if max_entropy > 0
        norm_entropy = shannon_entropy / max_entropy;
    else
        % Handle cases with a single bin, where max_entropy is zero.
        norm_entropy = 0;
    end
end