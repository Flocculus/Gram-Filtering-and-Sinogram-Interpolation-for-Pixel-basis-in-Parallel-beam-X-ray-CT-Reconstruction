classdef GenSheppLoganPhantom
  
  
  methods
	function p = eval(obj, xStart, xStep, xSize)
        p = phantom(xSize(1), xSize(2));
	end
  end
  
end