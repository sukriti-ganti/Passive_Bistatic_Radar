function s = setfield_local(s, pairs)
%SETFIELD_LOCAL  Set several struct fields from a flat {name, value, ...} list.
%
%   S = SETFIELD_LOCAL(S, {'a', 1, 'b', 2}) returns S with S.a = 1, S.b = 2.
%
%   A convenience for building the many small parameter structs that the
%   link budget sweeps need, without writing four assignment lines each time.

for i = 1:2:numel(pairs)
    s.(pairs{i}) = pairs{i + 1};
end
end
