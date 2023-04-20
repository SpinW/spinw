function C = mmat(A,B,dim)
% Simple matrix multiplication of multidimensional arrays.
%
% Input:
% A, B      Multidimensional input arrays.
% dim       Contains two number, that selects two dimensions.
%
% The multiplication is a standard matrix multiplication. The two matrices
% are selected by dim:
%   AB = A(... dim(1) ... dim(2) ...) * B(... dim(1) ... dim(2) ...)
% The necessary condition that the multiplication can be performed:
%   size(A,dim(2)) = size(B,dim(1))
%
% Singleton dimensions in both A and B matrices are supported.
%
% The default value for dim is dim = [1 2].
%
% Examples:
% For 2D matrices mmat is identical to the Matlab built-in multiplication:
% A = [1 2 3 4];
% B = [1;2;3;4];
% C = mmat(A,B)
%
% C will be 30.
%
% For multidimensional arrays:
% A = repmat([1 2 3 4],[1 1 5]);
% B = [1 2 3 4]';
% C = mmat(A,B)
% C will be an array with dimensions of 1x1x5 and every element is 30.
%

if nargin == 0
    help mmat
    return
end

if (nargin < 3)
    dim = [1 2];
end
isdefaultdim = all(dim == [1 2]);

if numel(dim)~=2
    error('mmat:WrongInput','dim has to be a two element array!');
end

if size(A,dim(2)) ~= size(B,dim(1))
    error('mmat:WrongInput','Wrong input matrix sizes!');
end

nDA = ndims(A);
nDB = ndims(B);
nD = max(nDA,nDB);

if nD == 2 && isdefaultdim
    C = A * B;
    return
end

nA = [size(A),ones(1,nD-nDA)]; nA = nA(dim); 
nB = [size(B),ones(1,nD-nDB)]; nB = nB(dim);

% Array is double float which is 8 bytes per element, but this needs to be doubled
% (16 bytes/element) as sum / bsxfun do not operate in-place and need a temp matrix
neededMem = (numel(A)*nB(2) + numel(B)*nA(1)) * 16;
if sw_freemem < neededMem && isdefaultdim
    % Not enough memory to expand matrix to use bsxfun; use slow loop instead
    szA = size(A); if numel(szA) > 3, A = reshape(A, [szA(1:2) prod(szA(3:end))]); end
    szB = size(B); if numel(szB) > 3, B = reshape(B, [szB(1:2) prod(szB(3:end))]); end
    if numel(szA) == 2 && numel(szB) > 2
        output_shape = [szA(1) szB(2:end)];
        C = zeros([szA(1) szB(2) size(B,3)]);
        for ii = 1:size(B,3)
            C(:,:,ii) = A * B(:,:,ii);
        end
    elseif numel(szA) > 2 && numel(szB) == 2
        output_shape = [szA(1) szB(2) szA(3:end)];
        C = zeros([szA(1) szB(2) size(A,3)]);
        for ii = 1:size(A,3)
            C(:,:,ii) = A(:,:,ii) * B;
        end
    elseif size(A, 3) == size(B, 3)
        output_shape = [szA(1) szB(2:end)];
        C = zeros([szA(1) szB(2) size(A,3)]);
        for ii = 1:size(A,3)
            C(:,:,ii) = A(:,:,ii) * B(:,:,ii);
        end
    else
        error('Extra dimensions do not agree');
    end
    szC = size(C);
    if numel(szC) ~= numel(output_shape) || ~all(szC == output_shape)
        C = reshape(C, output_shape);
    end
    return
end

% form A matrix
% (nA1) x (nA2) x nB2
A = repmat(A,[ones(1,nD) nB(2)]);
% form B matrix
% nA1 x (nB1) x (nB2)
idx = 1:nD+1; idx([dim end]) = idx([end dim]);
repv = ones(1,nD+1); repv(dim(1)) = nA(1);

B = repmat(permute(B,idx),repv);

% multiply with expanding along singleton dimensions
C = sumsym(bsxfunsym(@times,A,B),dim(2));


idx2 = 1:nD+1; idx2([dim end]) = idx2([dim(1) end dim(2)]);

% permute back the final result to the right size
C = permute(C,idx2);

end
