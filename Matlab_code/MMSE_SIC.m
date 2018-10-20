function estMMSE=MMSE_SIC(r,H,sigma2,sigmas2,N)
%   MMSE_SIC detector in Massive MIMO
%   written by Amen Memmi
    modOrd = 2;
    % Create PSK modulator and demodulator System objects
    pskModulator   = comm.PSKModulator(...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitInput',         true);
    pskDemodulator = comm.PSKDemodulator( ...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitOutput',        true);
    estMMSE = zeros(N*modOrd, 1);
    orderVec = 1:N;
    k = N+1;
    % Start MMSE nulling loop
        for n = 1:N
            % Shrink H to remove the effect of the last decoded symbol
            H = H(:, [1:k-1,k+1:end]);
            % Shrink order vector correspondingly
            orderVec = orderVec(1, [1:k-1,k+1:end]);
            % Select the next symbol to be decoded
            G = (H'*H + ((N-n+1)*sigma2/sigmas2)*eye(N-n+1)) \ eye(N-n+1); % Same as inv(H'*H + Sigma_n/Es *I ), but faster
            [~, k] = min(diag(G)); %select smallest diagonal element in the covariance matrix
            symNum = orderVec(k);

            decBits = pskDemodulator(G(k,:) * H' * r);
            estMMSE(2 * (symNum-1) + (1:modOrd)) = decBits;

            if n < N
                r = r - H(:, k) * pskModulator(decBits);
            end
        end

end