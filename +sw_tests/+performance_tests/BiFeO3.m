function BiFeO3()
    % setup
    bfo = spinw;
    bfo.genlattice('lat_const', [5.58 5.58 13.86], 'angled', [90 90 120], ...
                   'sym', 'R 3 c');
    bfo.addatom('r', [0 0 0.2212], 'S', 2.5, 'label', 'MFe3');
    bfo.gencoupling('maxDistance', 20)
    bfo.addmatrix('label', 'J1', 'value', 0.69029);
    bfo.addmatrix('label', 'J2', 'value', 0.2)
    bfo.addcoupling('mat', 'J1', 'bond', 1)
    bfo.addcoupling('mat', 'J2', 'bond', 2)
    bfo.addmatrix('label', 'D', 'value', [1 -1 0]*0.185)
    bfo.addcoupling('mat', 'D', 'bond', 2)
    bfo.optmagk('kbase', [1; 1; 0], 'seed', 1);
    bfo.optmagsteep('random',false,'TolX', 1e-12);
    
    % test parameters
    % add omega tol for imag eigenvalues (ignored if hermit=0)
    spinwave_args_common = {{[-1/2 0 0], [0 0 0], [1/2 1/2 0], 30}, ...
                             'sortMode', false, 'hermit', 0, ...
                             'omega_tol', 0.05};
    egrid_args = {'component','Sperp','Evect',0:0.1:5, 'imagChk', 0};
    inst_args = {'dE',0.1};
    
    % do a spin wave calculation for incommensurate case ([nExt=[1,1,1]) 
    % and commensurate nExt=0.01 - which corresponds to a supercell of 
    % [11 11 1]
    for do_supercell = 0:1
        if do_supercell
            nExt = 0.01;
            bfo.genmagstr('nExt', nExt)
            mexs = 1;  % requires ~100GB RAM with no mex
        else
            nExt = [1,1,1];
            mexs = 0:1;
        end
        for mex = mexs
            swpref.setpref('usemex', logical(mex))
            if mex && ~do_supercell
                hermits = 0:1;
            else
                % chol throws error as matrix is not positive-definite
                hermits = 0; 
            end
            for hermit = hermits
                for optmem = 0:5:10
                    spinwave_args = [spinwave_args_common, ...
                                 'hermit', hermit, 'optmem', optmem];
                    test_name = [mfilename(), '_mex_', ...
                                 num2str(swpref.getpref('usemex').val), ...
                                 '_nExt_', regexprep(num2str(nExt), ' +', '_'), ...
                                 '_hermit_', num2str(hermit), ...
                                 '_optmem_', num2str(optmem)];
                    sw_tests.utilities.profile_spinwave(test_name, bfo, ...
                                                        spinwave_args,...
                                                        egrid_args, ...
                                                        inst_args, ...
                                                        0:1);
                end
            end
        end
    end
    


end