classdef HookToEdgeFinder < handle   
    properties
        hookIndices
        connectionType
    end
    
    methods
        function obj = HookToEdgeFinder(edgeCodes, connectionType)
             if nargin > 0          
                obj.hookIndices = edgeCodes;
                
                obj.connectionType = [];
                if nargin > 1
                    obj.connectionType = connectionType;
                end
             end
        end
        
        function I = compute(obj, hookIndex)  
            [i, j] = find(obj.hookIndices == hookIndex);
            I = i;
            
            if ~isempty(obj.connectionType)
                T = obj.connectionType(i);
                
                switchMask = j == 2;
                
                oneMask = switchMask & T == 1;
                twoMask = switchMask & T == 2;
                
                T(oneMask) = 2;
                T(twoMask) = 1;
                
                I = [I T];
            end
        end
    end
end

