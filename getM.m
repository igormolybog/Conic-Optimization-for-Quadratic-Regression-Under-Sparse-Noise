function M = getM(V, ndx)
    len = size(V, 1);
    [ndxx,ndxy] = find(ndx);
    M = sparse(len, len);
    e = speye(len);
    for t = 1:size(ndxx,1)
        i = ndxx(t);
        j = ndxy(t);
        vv = V(i)*conj(V(i)) + V(j)*conj(V(j));
        M11 = 1 - V(i)*conj(V(i))/vv;
        M22 = 1 - V(j)*conj(V(j))/vv;
        M12 = - V(i)*conj(V(j))/vv;
        M21 =  - V(j)*conj(V(i))/vv;
        Mij = sparse([M11 M12;
                        M21 M22]);
        M = M + [e(:,i) e(:,j)]*Mij*[e(:,i) e(:,j)].';
    end
end

