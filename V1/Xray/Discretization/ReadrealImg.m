classdef ReadrealImg
    properties
        Img;
    end  
  
  methods
    function obj = ReadrealImg()
        load('./realimg/realimg.mat')
        tmp = realimg;
        obj.Img = (tmp - min(min(tmp)))./(max(max(tmp)) - min(min(tmp)));
        
    end
	function p = eval(obj, xStart, xStep, xSize)
        rate = floor(512/xSize(1));
        p = obj.Img(1:rate:xSize(1)*rate,1:rate:xSize(1)*rate);
	end
  end
  
end