classdef ObjectArray
   properties
      Value
   end
   methods
      function obj = ObjectArray(F)
         if nargin ~= 0
            % obj = Hook(width, pos2d, rotZAngleRadians)
            numHooks = size(F,1);
            obj(numHooks, 1) = ObjectArray;
            for i = 1 : numHooks
                obj(i, 1).Value = Hook(F(i,1), F(i, 2:3)', F(i, 4));
            end
         end
      end
   end
end