function aveTemp = layerTemp(top,bottom,wtrT,depths,bthA,bthD)

%----Author: Jordan S Read, 2009----
% adjusted 14 August 2013, from layerDensity to calculate volumetric
% average temperature for the layer
%
%Finds the average temperature thermal layer, bounded by "top" and "bottom"
% where distances are measured from the surface of the water column.
%
%Input:
%   -wtrT: temperature values (celsius)
%   -depths: corresponding depth values (m)
%   -top: surface to top of metalimnion
%   -bottom: surface to bottom of metalimnion
%   -bthA: bathymetric areas (m2)
%   -bthD: corresponding bathymetric depth values (m)
%
%Output:
%   -aveTemp: volumetric average layer temperature (oC)


if top > bottom
    error('bottom depth must be greater than top')
end
if ne(length(wtrT),length(depths))
    error(['water temperature array must be the same '...
        'length as the depth array'])
elseif nargin < 4
    error('not enough input arguments')
elseif any(isnan([wtrT depths bthA bthD]))
    error('input arguments must be numbers')
end

% if bathymetry has negative values, intepolate to 0
if lt(min(bthD),0)
    useI = ge(bthD,0);
    if ~eq(bthD,0)
        depT = [0 bthD(useI)];
    else
        depT = bthD(useI);
    end
    bthA = interp1(bthD,bthA,depT);
    bthD = depT;
end

dz = 0.1;   % (m)

numD = length(wtrT);
if max(bthD) > depths(numD)
    wtrT(numD+1) = wtrT(numD);
    depths(numD+1) = max(bthD);
elseif max(bthD) < depths(numD)
    bthD = [bthD depths(numD)];
    bthA = [bthA 0];
end
if min(bthD)<depths(1)
    wtrT = [wtrT(1) wtrT];
    depths = [min(bthD) depths];
end

[~,Io] = min(depths);
Ao = bthA(Io);
if eq(Ao,0)
    error('Surface area cannot be zero, check *.bth file')
end

%interpolates the bathymetry data
layerD = top:dz:bottom;
layerT  = interp1(depths,wtrT,layerD);
layerA  = interp1(bthD,bthA,layerD);
layerV = layerA.*dz;

aveTemp = sum(layerT.*layerV)/(sum(layerV));


end
