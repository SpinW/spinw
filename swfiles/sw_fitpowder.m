classdef sw_fitpowder < handle & matlab.mixin.SetGet
    % class to fit powders
    
    %
    % [spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html)
    %

    properties (SetObservable)
       swobj
       data
    end
        
    methods
        function obj = sw_fitpowder(swobj, data)
            % constructor
            obj.swobj = swobj;
            obj.data = data;
        end

        function result = fit(obj, fit_struct)
            result = obj.swobj.fitpow(obj.data, fit_struct);
        end
    end
end