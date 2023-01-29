function F=get_distribution(X)
[F, ~ ,lim]= histi(int32(X));             %distribution and its range
F=[zeros(1,lim(1)-1) F];                %ensure distribution starts at 1
end
