function cax = caxRange( image, minThresh, maxThresh )

[N,x] = histcounts( image, 100000 );
x = 0.5*(x(2:end)+x(1:end-1));
Nc = cumsum(N) / sum(N);

cax = [ x(find(Nc>=minThresh,1,'first')), ...
        x(find(Nc>=maxThresh,1,'first'))];
end