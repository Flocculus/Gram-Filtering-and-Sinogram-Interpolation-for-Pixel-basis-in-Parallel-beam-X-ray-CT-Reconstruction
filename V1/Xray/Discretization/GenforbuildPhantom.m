classdef GenforbuildPhantom
  
  
  methods
	function p = eval(obj, xStart, xStep, xSize)
        p = forbild_gen(xSize(1), xSize(2), 1, 1);
	end
  end
  
end