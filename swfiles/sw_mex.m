function sw_mex(varargin)
% compiles and tests the mex files
% 
% ### Syntax
% 
% `sw_mex(Name,Value)`
% 
% ### Description
% 
% `sw_mex(Name,Value)` compiles and tests the generated mex files. The
% compiled mex files will speed up the [spinw.spinwave] function. The
% expected speedup is larger for smaller magnetic unit cells. Once the mex
% files are compiled, use the `swpref.setpref('usemex',true)` command to
% switch on using mex files in [spinw.spinwave].
% 
% ### Name-Value Pair Arguments
%
% `'compile'`
% : If `false`, mex files will not be compiled. Default is
%   `true`.
% 
% `'test'`
% : If `true`, the compiled .mex files will be tested. Default is
%   `false`.
% 
% `'swtest'`
% : If `true`, 3 spin wave calculation will run with and without .mex
%   files and the results will be compared. Default is `false`.
% 
% ### See Also
% 
% [swpref]
%

inpForm.fname  = {'test' 'compile' 'swtest' };
inpForm.defval = {false  true      false    };
inpForm.size   = {[1 1]  [1 1]     [1 1]    };

param = sw_readparam(inpForm, varargin{:});

if param.compile
    % save current folder
    aDir = pwd;
    eig_omp_dir = [sw_rootdir filesep 'external' filesep 'eig_omp'];
    chol_omp_dir = [sw_rootdir filesep 'external' filesep 'chol_omp'];
    mtimesx_dir = [sw_rootdir filesep 'external' filesep 'mtimesx'];
    % compile the mex files
    if ispc
        cd(eig_omp_dir);
        mex('-v','-largeArrayDims','eig_omp.cpp','-lmwlapack',...
            'COMPFLAGS=$COMPFLAGS /openmp','LINKFLAGS=$LINKFLAGS /nodefaultlib:vcomp "$MATLABROOT\bin\win64\libiomp5md.lib"')
        cd(chol_omp_dir);
        mex('-v','-largeArrayDims','chol_omp.cpp','-lmwlapack',...
            '-lmwblas','COMPFLAGS=$COMPFLAGS /openmp',...
            'LINKFLAGS=$LINKFLAGS /nodefaultlib:vcomp "$MATLABROOT\bin\win64\libiomp5md.lib"')
        cd(mtimesx_dir);
        mex('-v','-largeArrayDims','sw_mtimesx.c','-lmwblas','COMPFLAGS=$COMPFLAGS /openmp',...
            'LINKFLAGS=$LINKFLAGS /nodefaultlib:vcomp "$MATLABROOT\bin\win64\libiomp5md.lib"')
    elseif ismac
        % add =libiomp5 after -fopenmp?
        cd(eig_omp_dir);
        mex('-v','-largeArrayDims','eig_omp.cpp', 'CXX_FLAGS="-Xclang -fopenmp -pthread"', ...
            'LDFLAGS="$LDFLAGS -L$MATLABROOT/sys/os/maci64 -liomp5 -lmwblas -lmwlapack"', ...
            'CXXOPTIMFLAGS="$CXXOPTIMFLAGS -Xclang -fopenmp"', "-I/usr/local/opt/libomp/include");
        cd(chol_omp_dir);
        mex('-v','-largeArrayDims','chol_omp.cpp', 'CXX_FLAGS="-Xclang -fopenmp -pthread"', ...
            'LDFLAGS="$LDFLAGS -L$MATLABROOT/sys/os/maci64 -liomp5 -lmwblas -lmwlapack"', ...
            'CXXOPTIMFLAGS="$CXXOPTIMFLAGS -Xclang -fopenmp"', "-I/usr/local/opt/libomp/include");
        cd(mtimesx_dir);
        mex('-v','-largeArrayDims','sw_mtimesx.c', 'CXX_FLAGS="-Xclang -fopenmp -pthread"', ...
            'LDFLAGS="$LDFLAGS -L$MATLABROOT/sys/os/maci64 -liomp5 -lmwblas"', ...
            'CXXOPTIMFLAGS="$CXXOPTIMFLAGS -Xclang -fopenmp"', "-I/usr/local/opt/libomp/include");
    else
        % linux?
        cd(eig_omp_dir);
        mex('-v','-largeArrayDims','eig_omp.cpp','-lmwlapack','CXXFLAGS=$CXXFLAGS -fopenmp -pthread','LDFLAGS=$LDFLAGS -fopenmp')
        cd(chol_omp_dir);
        mex('-v','-largeArrayDims','chol_omp.cpp','-lmwlapack','-lmwblas','CXXFLAGS=$CXXFLAGS -fopenmp -pthread','LDFLAGS=$LDFLAGS -fopenmp')
        cd(mtimesx_dir);
        mex('-v','-largeArrayDims','sw_mtimesx.c','-lmwblas','CXXFLAGS=$CXXFLAGS -fopenmp -pthread','LDFLAGS=$LDFLAGS -fopenmp')
    end
    % return back to original folder
    cd(aDir);
end

if param.test
    
    % Parameters
    n  = 8;   % Matrix size == n x n
    nb = 20000;  % Number of matrices nb
    
    % Test matrices (real nonsymmetric, complex nonhermitian, real symmetric, complex hermitian)
    dge = rand(n);
    zge = rand(n)+1i*rand(n);
    dsy = rand(n); dsy = triu(dsy)+triu(dsy,1)';
    zhe = rand(n)+1i*rand(n); zhe = triu(zhe,1)+triu(zhe,1)'+diag(real(diag(zhe)));
    
    % Functionality tests - should not give any errors, just want to check that correct number of outputs given.
    % {
    [V, E]=eig_omp(zhe);
    D = eig_omp(zhe);
    eig_omp(zhe);
    [V, E]=eig_omp(dsy);
    D = eig_omp(dsy);
    eig_omp(dsy);
    [V, E]=eig_omp(zge);
    D = eig_omp(zge);
    eig_omp(zge);
    [V, E]=eig_omp(dge);
    D = eig_omp(dge);
    eig_omp(dge);
    %}
    
    % Numerical test vs built-in eig()
    % {
    %disp(sprintf('\nTest vs builtin:  \teig_omp Eigenvalues\teig Eigenvalues  \tEigenvectors differences'));
    %[V1,E1]=eig_omp(dsy); [V2,E2]=eig(dsy); disp(sprintf('Real symmetric   \t%10.5f\t\t%10.5f\t\t%10.5f\n',[diag(E1),diag(E2),sum(abs(V1)-abs(V2),2)]'));
    %[V1,E1]=eig_omp(zhe); [V2,E2]=eig(zhe); disp(sprintf('Complex hermitian\t%10.5f\t\t%10.5f\t\t%10.5f\n',[diag(E1),diag(E2),sum(abs(V1)-abs(V2),2)]'));
    %[V1,E1]=eig_omp(dge); [V2,E2]=eig(dge); disp(sprintf('Real general     \t%10.5f +i%10.5f\t%10.5f +i%10.5f\t%10.5f\n',[real(diag(E1)),imag(diag(E1)),real(diag(E2)),imag(diag(E2)),sum(abs(V1)-abs(V2),2)]'));
    %[V1,E1]=eig_omp(zge); [V2,E2]=eig(zge); disp(sprintf('Complex general  \t%10.5f +i%10.5f\t%10.5f +i%10.5f\t%10.5f\n',[real(diag(E1)),imag(diag(E1)),real(diag(E2)),imag(diag(E2)),sum(abs(V1)-abs(V2),2)]'));
    fprintf('\nTest vs builtin:  \tEigenvalues diffs  \tEigenvectors differences\n');
    fprintf('------------------------------------------------------------------------\n');
    [V1,E1]=eig_omp(dsy); [V2,E2]=eig(dsy); fprintf('Real symmetric   \t% 10.5g\t\t% 10.5g\n',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2))));
    [V1,E1]=eig_omp(zhe); [V2,E2]=eig(zhe); fprintf('Complex hermitian\t% 10.5g\t\t% 10.5g\n',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2))));
    [V1,E1]=eig_omp(dge); [V2,E2]=eig(dge); fprintf('Real general     \t% 10.5g\t\t% 10.5g\n',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2))));
    [V1,E1]=eig_omp(zge); [V2,E2]=eig(zge); fprintf('Complex general  \t% 10.5g\t\t% 10.5g\n',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2))));
    fprintf('\nConsistency tests:\tV*E*inv(V)-A\tV*E*V''-A\tV*V''-I\t\tV''*V-I\n');
    fprintf('------------------------------------------------------------------------------\n');
    [V,E]=eig_omp(dsy); fprintf('Real symmetric   \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g\n',sum(sum(V*E/V-dsy)),sum(sum(V*E*V'-dsy)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n))));
    [V,E]=eig_omp(zhe); fprintf('Complex hermitian\t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g\n',sum(sum(V*E/V-zhe)),sum(sum(V*E*V'-zhe)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n))));
    [V,E]=eig_omp(dge); fprintf('Real general     \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g\n',sum(sum(V*E/V-dge)),sum(sum(V*E*V'-dge)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n))));
    [V,E]=eig_omp(zge); fprintf('Complex general  \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g\n',sum(sum(V*E/V-zge)),sum(sum(V*E*V'-zge)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n))));
    %}
    
    matfn{1} = @(a,b) triu(a)+triu(a,1)';
    matfn{2} = @(a,b) triu(a,1)+1i*triu(b,1)+triu(a,1)'-1i*triu(b,1)'+diag(diag(a));
    matfn{3} = @(a,b) a;
    matfn{4} = @(a,b) a+1i*b;
    
    % Numerical test of sorting functionality
    % {
    clear V1; clear V2; clear E1; clear E2; clear B;
    nt = 10; if(nt>nb); nt=nb; end
    for tp=1:4
        for ii=1:nt; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
        for ii=1:nt
            [vv,D]=eig(B(:,:,ii)); [~,isr]=sort(real(diag(D)),'descend');
            V1(:,:,ii)=vv(:,isr); E1(:,ii)=diag(D); E1(:,ii)=E1(isr,ii);
        end
        [V2,E2]=eig_omp(B,'sort',-1);
        sortres(:,tp) = [sum(sum(sum(abs(V1)-abs(V2)))) sum(sum(abs(E1-E2)))];
    end
    fprintf('\nSort Test:  \t\tEigenvecDiffs   EigenvalDiffs\n');
    fprintf('-----------------------------------------------------\n');
    fprintf('Real symmetric   \t% 12.5g\t% 12.5g\n',sortres(:,1));
    fprintf('Complex hermitian\t% 12.5g\t% 12.5g\n',sortres(:,2));
    fprintf('Real general     \t% 12.5g\t% 12.5g\n',sortres(:,3));
    fprintf('Complex general  \t% 12.5g\t% 12.5g\n',sortres(:,4));
    %}
    
    % Numerical test of orthogonalisation functionality
    % {
    clear V1; clear V2; clear E1; clear E2; clear B;
    for tp=1:4
        for ii=1:nb; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
        %for ii=1:nb;
        %  [vv,D]=eigorth(B(:,:,ii),eps);  [~,isr]=sort(real(D));
        %  V1(:,:,ii)=vv(:,isr); E1(:,ii) = D(isr);
        %end;
        tic; [V1,E1]=eigorth(B,eps); t1(tp)=toc;
        tic; [V2,E2]=eig_omp(B,'orth','sort'); t2(tp)=toc;
        for ii=1:nb; [~,isr]=sort(real(E1(:,ii))); V1(:,:,ii)=V1(:,isr,ii); E1(:,ii) = E1(isr,ii); end
        orthres(:,tp) = [sum(sum(sum(abs(V1)-abs(V2)))) sum(sum(abs(E1-E2)))];
    end
    fprintf('\nOrth Test:  \t\tEigenvecDiffs   EigenvalDiffs   Time(eigorth)   Time(eig_omp)         Speedup\n');
    fprintf('-----------------------------------------------------------------------------------------------------\n');
    fprintf('Real symmetric   \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',orthres(:,1),t1(1),t2(1),t1(1)/t2(1));
    fprintf('Complex hermitian\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',orthres(:,2),t1(2),t2(2),t1(2)/t2(2));
    fprintf('Real general     \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',orthres(:,3),t1(3),t2(3),t1(3)/t2(3));
    fprintf('Complex general  \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',orthres(:,4),t1(4),t2(4),t1(4)/t2(4));
    %}
    
    % Speed tests
    % {
    clear E1; clear E2; clear B; clear ch; clear ch2;
    for tp=1:4
        for ii=1:nb; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
        tic; [V,E]=eig_omp(B); t1(tp)=toc;
        for ii=1:nb
            ch(ii)=sum(sum(V(:,:,ii)*diag(E(:,ii))*V(:,:,ii)'-B(:,:,ii))); ch2(ii)=sum(sum(V(:,:,ii)'*V(:,:,ii)-eye(size(V(:,:,ii)))));
        end; check1(tp,:)=[sum(ch) sum(ch2)];
        tic; for ii=1:size(E,2); [V(:,:,ii),ee]=eig(B(:,:,ii)); E2(:,ii)=diag(ee); end; t2(tp)=toc;
        for ii=1:nb
            ch(ii)=sum(sum(V(:,:,ii)*diag(E(:,ii))*V(:,:,ii)'-B(:,:,ii))); ch2(ii)=sum(sum(V(:,:,ii)'*V(:,:,ii)-eye(size(V(:,:,ii)))));
        end; check2(tp,:)=[sum(ch) sum(ch2)];
        check3(tp)=sum(sum(E-E2));
    end
    fprintf('\nSpeed Test:  \t\tTime(eig_omp)\t   Time(eig)\t   Speedup    Eigenvalue diffs   Sum(V*E*V''-A)\tSum(V''*V-I)\n');
    fprintf('-------------------------------------------------------------------------------------------------------------------\n');
    fprintf('Real symmetric   \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t1(1),t2(1),t2(1)/t1(1),abs(check3(1)),abs(check1(1,1)),abs(check1(1,2)));
    fprintf('Complex hermitian\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t1(2),t2(2),t2(2)/t1(2),abs(check3(2)),abs(check1(2,1)),abs(check1(2,2)));
    fprintf('Real general     \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t1(3),t2(3),t2(3)/t1(3),abs(check3(3)),abs(check1(3,1)),abs(check1(3,2)));
    fprintf('Complex general  \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t1(4),t2(4),t2(4)/t1(4),abs(check3(4)),abs(check1(4,1)),abs(check1(4,2)));
    %}
    
    % Test of Cholesky factorisation routines
    % {
    chlfn{1} = @(a,b) triu(a)+eye(size(a,1))*size(a,1);
    chlfn{2} = @(a,b) triu(a)+1i*triu(b,1)+eye(size(a,1))*size(a,1);
    chlfn{3} = @(a,b) tril(a)+eye(size(a,1))*size(a,1);
    chlfn{4} = @(a,b) tril(a)+1i*tril(b,-1)+eye(size(a,1))*size(a,1);
    chlfn{5} = @(a,b) triu(a);
    chlfn{6} = @(a,b) triu(a)+1i*triu(b,1);
    chlfn{7} = @(a,b) tril(a);
    chlfn{8} = @(a,b) tril(a)+1i*tril(b,-1);
    
    clear mm; clear K1;
    for tp=1:4
        % These matrices should be positive definite (with large diagonal elements)
        for ii=1:nb; mm(:,:,ii) = chlfn{tp}(rand(n),rand(n)); end
        tic; for ii=1:nb; if(tp>2); K1(:,:,ii) = chol(mm(:,:,ii),'lower'); else; K1(:,:,ii) = chol(mm(:,:,ii)); end; end; t1(tp)=toc;
        tic; if(tp>2); K2 = chol_omp(mm,'lower'); else; K2 = chol_omp(mm); end; t2(tp)=toc;
        ch(tp) = sum(sum(sum(abs(K1-K2))));
    end
    for tp=5:8
        % These matrices should not be positive definite
        for ii=1:nb; mm(:,:,ii) = chlfn{tp}(rand(n),rand(n)); end
        tic; for ii=1:nb; if(tp>6); [~,p1(ii)] = chol(mm(:,:,ii),'lower'); else; [~,p1(ii)] = chol(mm(:,:,ii)); end; end; t1(tp)=toc;
        tic; if(tp>6); [~,p2] = chol_omp(mm,'lower'); else; [~,p2] = chol_omp(mm); end; t2(tp)=toc;
        ch(tp) = sum(sum(p1-p2));
    end
    fprintf('\nCholesky Test:  \tTime(chol_omp)\t   Time(chol)\t   Speedup    CholFactor diffs   PositiveDef diffs\n');
    fprintf('------------------------------------------------------------------------------------------------------------------\n');
    fprintf('Real upper       \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(1)+t2(5),t1(1)+t1(5),(t1(1)+t1(5))/(t2(1)+t2(5)),ch(1),ch(5));
    fprintf('Real lower       \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(2)+t2(6),t1(2)+t1(6),(t1(2)+t1(6))/(t2(2)+t2(6)),ch(2),ch(6));
    fprintf('Complex upper    \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(3)+t2(7),t1(3)+t1(7),(t1(3)+t1(7))/(t2(3)+t2(7)),ch(3),ch(7));
    fprintf('Complex lower    \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(4)+t2(8),t1(4)+t1(8),(t1(4)+t1(8))/(t2(4)+t2(8)),ch(4),ch(8));
    
    if(mod(n,2)==1); n=n+1; end
    clear mm; clear K1; clear inv1;
    g=diag([ones(1,n/2) -ones(1,n/2)]);
    for tp=1:4
        % These matrices should be positive definite (with large diagonal elements)
        for ii=1:nb; mm(:,:,ii) = chlfn{tp}(rand(n),rand(n)); end
        tic; for ii=1:nb
            if(tp>2); K = chol(mm(:,:,ii),'lower'); else; K = chol(mm(:,:,ii)); end
            inv1(:,:,ii)=inv(K); K1(:,:,ii)=K*g*K';
        end; t1(tp)=toc;
        tic; if(tp>2); [K2,inv2] = chol_omp(mm,'Colpa','lower'); else; [K2,inv2] = chol_omp(mm,'Colpa'); end; t2(tp)=toc;
        ch(tp) = sum(sum(sum(abs(K1-K2))));
        ch(tp+4) = sum(sum(sum(abs(inv1-inv2))));
    end
    fprintf('\nColpa Test:     \tTime(chol_omp)\t   Time(chol+eig)  Speedup    K*g*K'' diffs\t inv(K) diffs\n');
    fprintf('------------------------------------------------------------------------------------------------------------------\n');
    fprintf('Real upper       \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(1),t1(1),t1(1)/t2(1),ch(1),ch(5));
    fprintf('Real lower       \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(2),t1(2),t1(2)/t2(2),ch(2),ch(6));
    fprintf('Complex upper    \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(3),t1(3),t1(3)/t2(3),ch(3),ch(7));
    fprintf('Complex lower    \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\n',t2(4),t1(4),t1(4)/t2(4),ch(4),ch(8));
    
end

if param.swtest
    pref = swpref;
    
    % Antiferromagnetic square lattice (tutorial 4, Cu+1 S=1) - small system [3 spins, n=6]
    AFsq = spinw;
    AFsq.genlattice('lat_const',[3 3 10],'angled',[90 90 90],'sym',0)
    AFsq.addatom('r',[0 0 0],'S', 1,'label','Cu1','color','b');
    AFsq.gencoupling('maxDistance',9)
    AFsq.addmatrix('label','J1','value',1,'color','red');      AFsq.addcoupling('mat','J1','bond',1)
    AFsq.addmatrix('label','J2','value',-0.1,'color','green'); AFsq.addcoupling('mat','J2','bond',2)
    AFsq.addmatrix('label','D','value',diag([-0.2 0 0]));      AFsq.addaniso('D');
    AFsq.genmagstr('mode','helical','k',[1/2 1/2 0],'n',[0 0 1], 'S',[1; 0; 0],'nExt',[2 2 1]);
    AFsq.optmagsteep
    %plot(AFsq,'range',[2 2 0.5],'zoom',-1)
    % Runs test
    hkl = {[1/4 3/4 0] [1/2 1/2 0] [1/2 0 0] [3/4 1/4 0] [1 0 0] [3/2 0 0] 50000};
    nm = 15;   % Ensure same number of slices for all tests
    pref.usemex = true;
    tic; linespec_herm_mex      = AFsq.spinwave(hkl,'hermit',true,'optmem',nm,'fid',0);  t1=toc;
    pref.usemex = false;
    tic; linespec_herm_nomex    = AFsq.spinwave(hkl,'hermit',true,'optmem',nm,'fid',0); t2=toc;
    pref.usemex = true;
    tic; linespec_nonherm_mex   = AFsq.spinwave(hkl,'hermit',false,'optmem',nm,'fid',0); t3=toc;
    pref.usemex = false;
    tic; linespec_nonherm_nomex = AFsq.spinwave(hkl,'hermit',false,'optmem',nm,'fid',0);t4=toc;
    
    fprintf('             %16s  %16s  %16s  %16s\n','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex');
    fprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f\n',t1*5,t2*5,t3*5,t4*5);
    
    % Times
    %   ndlt811:         43.7453   72.0950  406.1366  611.1799
    %   eryenyo:         24.1470   85.5668  382.7930  563.6103
    %   ndl01wkc26243:   36.2106  107.8983  527.3301  796.2985
    %   simon_i9_32Gb:  113.8068  149.5165  159.9780  260.2178
    %   simonWork:      173.9701  231.1420  216.0472  455.3942
    
    % KCu3As2O7(OD)3 kagome (Nilsen PRB 89 140412) - Tutorial 18, medium sized [18 spins, n=36]
    J   = -2; Jp  = -1; Jab = 0.75; Ja  = -J/.66 - Jab; Jip = 0.01;
    hK = spinw;
    hK.genlattice('lat_const',[10.2 5.94 7.81],'angled',[90 117.7 90],'sym','C 2/m');
    hK.addatom('r',[0   0   0],'S',1/2,'label','MCu2','color','b');
    hK.addatom('r',[1/4 1/4 0],'S',1/2,'label','MCu2','color','k');
    hK.gencoupling;
    hK.addmatrix('label','J',  'color','r',   'value',J);   hK.addcoupling('mat','J','bond',1);
    hK.addmatrix('label','J''','color','g',   'value',Jp);  hK.addcoupling('mat','J''','bond',2);
    hK.addmatrix('label','Ja', 'color','b',   'value',Ja);  hK.addcoupling('mat','Ja','bond',3);
    hK.addmatrix('label','Jab','color','cyan','value',Jab); hK.addcoupling('mat','Jab','bond',5);
    hK.addmatrix('label','Jip','color','gray','value',Jip); hK.addcoupling('mat','Jip','bond',10);
    hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',[0.77 0 0.115],'next',[1 1 1]);
    optpar.func = @gm_planar;
    optpar.nRun = 5;
    optpar.xmin = [    zeros(1,6), 0.5 0 0.0, 0 0];
    optpar.xmax = [2*pi*ones(1,6), 1.0 0 0.5, 0 0];
    magoptOut = hK.optmagstr(optpar);
    kOpt = hK.mag_str.k(:)';
    hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',kOpt,'next',[1 1 1]);
    %plot(hK,'range',[2 2 0.3],'sSpin',2)
    
    % Runs test
    hkl = {[0 0 0] [1 0 0] 50000};
    nm = 10;  % Ensure same number of slices for all tests
    pref.usemex = true;
    tic; linespec_herm_mex = hK.spinwave(hkl,'hermit',true,'optmem',nm); t1=toc;
    pref.usemex = false;
    tic; linespec_herm_nomex = hK.spinwave(hkl,'hermit',true,'optmem',nm); t2=toc;
    pref.usemex = true;
    tic; linespec_nonherm_mex = hK.spinwave(hkl,'hermit',false,'optmem',nm); t3=toc;
    pref.usemex = false;
    tic; linespec_nonherm_nomex = hK.spinwave(hkl,'hermit',false,'optmem',nm); t4=toc;
    fprintf('             %16s  %16s  %16s  %16s\n','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex');
    fprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f\n',t1,t2,t3,t4);
    
    %{
[~,nSuperlat] = rat(hK.mag_str.k,0.01);
hK.genmagstr('mode','helical','next',nSuperlat)
hK.mag_str.k = [0 0 0];
out = optmagsteep(hK,'nRun',200); hK = out.obj;

hkl = {[0 0 0] [1 0 0] 40};
nm = 10;  % Ensure same number of slices for all tests
tic; linespec_herm_mex = hK.spinwave(hkl,'hermit',true,'useMex',true,'optmem',nm); t5=toc;
tic; linespec_herm_nomex = hK.spinwave(hkl,'hermit',true,'useMex',false,'optmem',nm); t6=toc;
tic; linespec_nonherm_mex = hK.spinwave(hkl,'hermit',false,'useMex',true,'optmem',nm); t7=toc;
tic; linespec_nonherm_nomex = hK.spinwave(hkl,'hermit',false,'useMex',false,'optmem',nm); t8=toc;
disp(sprintf('Supercell    % 16.6f  % 16.6f  % 16.6f    % 16.6f',t5,t6,t7,t8));
    %}
    
    % Times
    %   ndlt811:         26.4661   40.5760   60.6112   63.3341
    %   eryenyo:          7.6710   36.0402   56.3213   73.5568
    %   ndl01wkc26243:   17.2135   42.8147   74.9062  106.8541
    %   simon_i9_32Gb:   49.5562   55.9125   54.1221   65.5144
    %   simonWork:       68.3912   84.6127   75.3424  100.9470
    
    % Bi4Fe5O13F - large(ish) system [80 spins, n=160]
    ff=11.6*2.5; Jc1 = 34/ff; Jc2 = 20/ff; Jab1 = 45/ff; Jab2 = 74/ff; Jd = 191/ff;
    bfof = spinw();
    bfof.genlattice('lat_const',[8.29950 8.29950 18.05730],'angled',[90 90 90],'sym','P 42/m b c');
    bfof.addatom('r',[0.5  0.   0.0800],'S',2.5,'color',[0 0 255],'label','Fe1');
    bfof.addatom('r',[0.8515 0.8388 0], 'S',2.5,'color',[255 0 0],'label','Fe2');
    bfof.addatom('r',[0.5  0.   0.25  ],'S',2.5,'color',[0 0 128],'label','Fe1_3');
    bfof.gencoupling('maxBond',99,'maxDistance',10);
    bfof.addmatrix('value',Jc1,'label','Jc1','color','r');         bfof.addcoupling('mat','Jc1','bond',1);
    bfof.addmatrix('value',Jc2,'label','Jc2','color',[128 0 0]);   bfof.addcoupling('mat','Jc2','bond',2);
    bfof.addmatrix('value',Jab1,'label','Jab1','color','b');       bfof.addcoupling('mat','Jab1','bond',3);
    bfof.addmatrix('value',Jab2,'label','Jab2','color',[0 255 0]); bfof.addcoupling('mat','Jab2','bond',4);
    bfof.addmatrix('value',Jd,'label','Jd','color','k');           bfof.addcoupling('mat','Jd','bond',5);
    bfof.addmatrix('value',diag([0 0 0.2]),'label','D');           bfof.addaniso('D');
    S2a = -[4.05 -0.35 0]; S2b = -[-0.35 4.05 0]; S1a = [2.18 -2.53 0]; S1b = [2.53 2.18 0];
    S = [S1b; S1a; S1a; S1b; S1b; S1a; S1a; S1b; S2a; -S2a; -S2b; S2b; S2b; -S2b; S2a; -S2a; -S1b; -S1a; -S1b; -S1a];
    Sv=[S; -S; -S; S];
    bfof.genmagstr('mode','direct','S',Sv','nExt',[2 2 1]);
    %bfof.plot('range',[0 2; 0 2; 0.0 1.0]); view([45 90]); set(gcf,'Tag','');
    out = optmagsteep(bfof,'nRun',200);
    %plot(out.obj,'range',[0 2; 0 2; 0.0 1.0]); view([45 90]);
    bfof = out.obj;
    
    % Runs test
    nm = 60;  % For laptops with 16GB memory, to have the same number of slices
    %nm = 6;   % For the workstation with 50GB.
    hkl={[0 0 0] [1 1 0] [1 1 1] [0 0 1] 2000};
    pref.usemex = true;
    tic; linespec_herm_mex      = bfof.spinwave(hkl,'hermit',true,'optmem',nm); t1=toc;
    pref.usemex = false;
    tic; linespec_herm_nomex    = bfof.spinwave(hkl,'hermit',true,'optmem',nm); t2=toc;
    pref.usemex = true;
    tic; linespec_nonherm_mex   = bfof.spinwave(hkl,'hermit',false,'optmem',nm); t3=toc;
    pref.usemex = false;
    tic; linespec_nonherm_nomex = bfof.spinwave(hkl,'hermit',false,'optmem',nm); t4=toc;
    
    fprintf('             %16s  %16s  %16s  %16s\n','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex');
    fprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f\n',t1,t2,t3,t4);
    % Times
    %   ndlt811:        123.6812  139.6476  197.1951  363.2818
    %   eryenyo:        127.4535  135.8621  220.9112  389.4686
    %   ndl01wkc26243:   67.8451  117.7316  146.2180  474.1894
    %   simon_i9_32Gb:  74.31700   84.7511  253.0013  342.0254
    %   simonWork:      75.40301  107.8941  216.4226  467.2630

end
end
