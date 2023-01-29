function En= makes_en(vt)

En=sum(vt);
En([1:find(~En,1,'first') find(~En,1,'last'):end])=0;       %remove clipped events

end

