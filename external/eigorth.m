function [V,D, warnOut] = eigorth(M, tol, useMex)
% orthogonal eigenvectors of defective eigenvalues
%
% [V, D, {Idx}] = EIGORTH(M, tol, useMex)
%
% Input:
%
% M         Stack of square matrices, dimensions are [nMat nMat nStack].
% tol       Tolerance, within two eigenvalue are regarded equal:
%               abs(real(diff(D))) < max([real(D) 1e-5])*tol
%           Default is 1e-5.
% useMex    If true a mex file will be used to speed up the calculation.
%           Default is false.
%
% Output:
%
% V         Matrix stack, every column contains an eigenvector of M in
%           every stack, dimensions are the same as M.
% D         Stack of vectors, that contains the eigenvalues of M,
%           dimensions are [nMat nStack].
% idx       Permutation indices of the eigenvalues and eigenvectors for
%           every stack, dimensions are [nMat nStack].
%
% See also eig, sortmode.
%

if nargin == 0
    help eigorth
    return
end

if nargout>2
    warnOut = false;
end

if nargin < 2
    tol = 1e-5;
end

if nargin < 3
    useMex = false;
end

nMat   = size(M,1);
nStack = size(M,3);

% Use OpenMP parallelised mex file if it exists
if nStack>1 && any(useMex)
    % eigenvalues are already orthogonalised by eig_omp
    [V, D] = eig_omp(M,'orth','sort','descend');
    return
end

V = zeros(nMat,nMat,nStack);
D = zeros(nMat,nStack);

for ii = 1:nStack
    [V(:,:,ii), Dtemp] = eig(M(:,:,ii));
    D(:,ii) = diag(Dtemp);
end


for jj = 1:nStack
    % sort eigenvalues-eigenvectors according to the real part of D
    [~, idxD] = sort(real(D(:,jj)),'descend');
    D(:,jj) = D(idxD,jj);
    Vi      = V(:,idxD,jj);
    % selects degenerate eigenvalues
    degIdx = [0 abs(real(diff(D(:,jj)'))) < max([real(D(:,jj)') 1e-5])*tol];
    % degIdx contains islands with increasing integers
    degVal = cumsum(diff([0 degIdx])==1);
    degVal = [degVal(2:end) degVal(end)];
    degIdx = (degIdx | [degIdx(2:end) 0]).*degVal;
    
    try
        % orthogonalise them
        for ii = 1:max(degIdx)
            Vi(:,degIdx==ii) = orth(Vi(:,degIdx==ii));
        end
        % permute back Vi
        %Vi(:,idxD) = Vi;
        
        V(:,:,jj) = Vi;
    catch %#ok<CTCH>
        % if there are not enough orthogonal vector, skip the orthogonalisation
        if nargout>2
            % catch warning
            warnOut = true;
        else
            warning('eigorth:NoOrth','Eigenvectors of defective eigenvalues cannot be orthogonalised!');
        end
    end
end

end