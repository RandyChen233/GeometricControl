function ts = arrangeT(waypts,T)
    x = waypts(:,2:end) - waypts(:,1:end-1);
    y = waypts(:,3:end) - waypts(:,2:end-1);
    y(:,end+1) = 0;
    dist = sum(x.^2,1).^0.5 + sum(y.^2,1).^0.5;
    k = T/sum(dist);
    ts = [0 cumsum(dist*k)];
end