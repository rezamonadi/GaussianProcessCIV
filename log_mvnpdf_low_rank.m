% log_mvnpdf_low_rank: efficiently computes
%
%   log N(y; mu, MM' + diag(d))

function log_p = log_mvnpdf_low_rank(y, mu, M, d)

  %disp(mean(y));
  %disp(mean(mu));
  %disp(mean(M(:)));
  %disp(mean(d));
  
  log_2pi = 1.83787706640934534;

  [n, k] = size(M);
 
  y = y - (mu);
  %d = ones(size(d)) * .001; %here
  d_inv = 1 ./ d;
  D_inv_y = d_inv .* y;
  % fprintf('S(d_inv):(%d,%d)\n', size(d_inv));
  % fprintf('S(y): (%d,%d)\n', size(y));
  % fprintf('S(M): (%d,%d)\n', size(M));
  D_inv_M = d_inv.*M;
  % use Woodbury identity, define
  %   B = (I + M' D^-1 M),
  % then
  %   K^-1 = D^-1 - D^-1 M B^-1 M' D^-1
  % fprintf('S(D_inv_M): (%d,%d)\n', size(D_inv_M));
  B = M' * D_inv_M;
  B(1:(k + 1):end) = B(1:(k + 1):end) + 1;
  
  L = chol(B);
  % C = B^-1 M' D^-1
  C = L \ (L' \ D_inv_M');
  % size(D_inv_M)
  % size(C)
  K_inv_y = D_inv_y - D_inv_M * (C * y);
  
  log_det_K = sum(log(d)) + 2 * sum(log(diag(L)));

  log_p = -0.5 * (y' * K_inv_y + log_det_K + n * log_2pi);
end
