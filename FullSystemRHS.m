function dpdt = FullSystemRHS(t,p)

global W

if size(p,1) == 1
    p = p.';
end

dpdt = W*p ;

return

end