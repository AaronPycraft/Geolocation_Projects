% https://www.movable-type.co.uk/scripts/latlong.html
% uses great circle method to compute distance
function [dm] = ll2km(ll1, ll2)
    lat1 = ll1(1);
    lon1 = ll1(2);
    lat2 = ll2(1);
    lon2 = ll2(2);
    %%
    deg2rad = pi/180;
    R = 6371e3; % metres
    %%
    phi1 = lat1 * deg2rad;
    phi2 = lat2 * deg2rad;
    deltaPhi    = (lat2-lat1) * deg2rad;
    deltaLambda = (lon2-lon1) * deg2rad;

    a = sin(deltaPhi/2) * sin(deltaPhi/2) ...
        + cos(phi1)*cos(phi2)*sin(deltaLambda/2)*sin(deltaLambda/2);
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    dm = R * c;
end