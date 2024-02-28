function fEncke = fEnckeDirect(r,rBody)
    nr = norm(r);
    nb = norm(rBody);
    
    q = (nr^2)/(nb^2)-2*dot(r,rBody)/(nb^2);
    fEncke = ((1+q)^1.5)-1;
end