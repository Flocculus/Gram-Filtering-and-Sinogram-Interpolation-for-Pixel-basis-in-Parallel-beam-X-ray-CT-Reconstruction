function y = enum(names)

y.names = names;
for i = 1:length(names)
	y.(names{i}) = i;
end