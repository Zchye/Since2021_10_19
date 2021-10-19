classdef BSType
%YXC
%This class enumerates all possible BS types and initializes some BS properties related to the types. The instance of this class may act as input to the constructor of the class hNRGNB.m. 

%Example
%We want to specify a BS as a Macrocell BS. In the script we type:
%   foo = BSType.Macro;
%Check the initialized property:
%   foo.height
%the result will be 25.

properties
height %meter
end

enumeration
    Macro (0)
    Micro (1)
end

methods
    function obj = BSType(val) %Constructor
        switch val
        case 0
            obj.height=25;
        case 1
            obj.height=10;
        otherwise
            obj.height=0;
        end

    end
end

end