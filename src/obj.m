function result = obj(P,F_til,F_0)
    result = norm(F_til * (diag(P) * F_til + F_0));
end