classdef Stack < handle
    %INTEGERSTACK Summary of this class goes here
    %   Detailed explanation goes here
    properties
        data
        height
        currentIndex
    end
    
    methods
        function obj = Stack(size, width)
            w = 1;
            
            if nargin > 1
                w = width;
            end
            
            obj.data = zeros(size, w);
            obj.height = size;
            obj.currentIndex = 0;
        end
        
        function push(obj, v)
            if (obj.currentIndex < obj.height)
                obj.currentIndex = obj.currentIndex + 1;
                obj.data(obj.currentIndex, :) = v;
            else
                error('Stack is full');
            end          
        end
        
        function v = pop(obj)
            if (obj.currentIndex > 0)
                v = obj.data(obj.currentIndex, :);
                obj.currentIndex = obj.currentIndex - 1;
            else
                error('Stack is empty');
            end
        end
        
        function c = count(obj)
            c = obj.currentIndex;
        end
    end
end