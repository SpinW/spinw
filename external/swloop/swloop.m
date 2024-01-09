% Calculates the inner loop of spinw/spinwave.m in parallel
%
% [omega, Sab, warn1, orthWarn0] = swavel3(param, hklExt, ...
%   ABCD, idxAll, ham_diag, dR, RR, S0, zed, FF, ...
%   bqdR, bqABCD, idxBq, bq_ham_d, ham_MF)
%
% This is a MEX-file for MATLAB.
%
% This code uses the Eigen matrix library for linear algebra
% C++ threads to calculate the inner loop of spinwave.m
% As this is different to the LAPACK/BLAS used by the Matlab
% code there will be numerical differences in the output
% You should double check that the Matlab and Mex code gives
% consistent results before running a long calculation.
%
% Original Author: M. D. Le  [duc.le@stfc.ac.uk]


% Equivalent Matlab code to C++ code follows:
%
% function [omega, Sab, warn1, orthWarn0] = swavel3(param, hklExt, ...
%     ABCD, idxAll, ham_diag, dR, RR, S0, zed, FF, bqdR, bqABCD, idxBq, bq_ham_d, ham_MF)
%
% nHkl = size(hklExt,2);
% nMagExt = param.nMagExt;
%
% % Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
% omega = zeros(2*nMagExt, nHkl);
%
% % empty Sab
% Sab = zeros(3,3,2*nMagExt,nHkl);
%
% orthWarn0 = false;
% warn1 = false;
%
% % Could replace the code below with a generator if Matlab has one...
% idx0 = 1:nHkl;
% nHklT = nHkl / param.nTwin;
% if param.incomm
%     nHkl0 = nHkl / 3 / param.nTwin;
%     for tt = 1:param.nTwin
%         t0 = (tt-1)*nHklT + 1;
%         idx0(t0:(t0+nHklT-1)) = repmat((t0+nHkl0):(t0+2*nHkl0-1), [1 3]);
%     end
% end
%
% % diagonal of the boson commutator matrix
% gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
% % boson commutator matrix
% gComm  = diag(gCommd);
%
% sqrtS0 = sqrt(S0 / 2);
%
% for jj = 1:nHkl
%     ExpF = exp(1i*(dR' * hklExt(:,jj)))';
%     ham = accumarray(idxAll, ABCD.*repmat(ExpF, [1 3]), [1 1]*2*nMagExt)' + diag(ham_diag);
%
%     if param.bq
%         ExpF = exp(1i*(bqdR' * hklExt(:,jj)))';
%         h_bq = accumarray(idxBq, bqABCD.*repmat(ExpF, [1 3]), [1 1]*2*nMagExt)' + diag(bq_ham_d);
%         ham = ham + h_bq;
%     end
%     if param.field
%         ham = ham + ham_MF{ceil(jj / nHklT)};
%     end
%
%     ham = (ham + ham') / 2;
%
%     if param.hermit
%         [K, posDef]  = chol(ham);
%         if posDef > 0
%             try
%                 % get tolerance from smallest negative eigenvalue
%                 tol0 = eig(ham);
%                 tol0 = sort(real(tol0));
%                 tol0 = abs(tol0(1));
%                 % TODO determine the right tolerance value
%                 tol0 = tol0*sqrt(nMagExt*2)*4;
%                 if tol0>param.omega_tol
%                     error('spinw:spinwave:NonPosDefHamiltonian','Very baaaad!');
%                 end
%                 try
%                     K = chol(ham+eye(2*nMagExt)*tol0);
%                 catch
%                     K = chol(ham+eye(2*nMagExt)*param.omega_tol);
%                 end
%                 warn1 = true;
%             catch PD
%                 error('spinw:spinwave:NonPosDefHamiltonian',...
%                     ['Hamiltonian matrix is not positive definite, probably'...
%                     ' the magnetic structure is wrong! For approximate'...
%                     ' diagonalization try the param.hermit=false option']);
%             end
%         end
%
%         K2 = K*gComm*K';
%         K2 = 1/2*(K2+K2');
%         % Hermitian K2 will give orthogonal eigenvectors
%         [U, D] = eig(K2);
%         D      = diag(D);
%
%         % sort modes accordign to the real part of the energy
%         [~, idx] = sort(real(D),'descend');
%         U = U(:,idx);
%         % omega dispersion
%         omega(:,jj) = D(idx);
%
%         % the inverse of the para-unitary transformation V
%         V = inv(K)*U*diag(sqrt(gCommd.*omega(:,jj))); %#ok<MINV>
%     else
%         gham = gComm * ham;
%         gham = gham + diag([(-nMagExt+0.5):nMagExt].*1e-12);
%         [V, D, orthWarn] = eigorth(gham,param.omega_tol);
%         orthWarn0 = orthWarn || orthWarn0;
%         % multiplication with g removed to get negative and positive energies as well
%         omega(:,jj) = D;
%         M = diag(gComm * V' * gComm * V);
%         V = V * diag(sqrt(1./M));
%     end
%
%     % TODO saveV / saveH
%
%     % Calculates correlation functions.
%     ExpF = exp(-1i * sum(repmat(hklExt(:,idx0(jj)),[1 nMagExt 1]) .* RR)) .* sqrtS0;
%     if param.formfact
%         ExpF = ExpF .* FF(:,jj)';
%     end
%     zedExpF = zeros(2*nMagExt, 3);
%     for i1 = 1:3
%         zedExpF(:, i1) = transpose([zed(i1,:) .* ExpF, conj(zed(i1,:)) .* ExpF]);
%     end
%     VExp = zeros(3, 2*nMagExt, 1);
%     for i1 = 1:3
%        VExp(i1,:,:) = V' * zedExpF(:, i1);
%     end
%     for i1 = 1:(2*nMagExt)
%        Sab(:,:,i1,jj) = VExp(:,i1) * VExp(:,i1)';
%     end
% end
