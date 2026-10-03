function write_csv(path, header, rows)
%WRITE_CSV  Write a comma separated file with a header line.
%
%   WRITE_CSV(PATH, HEADER, ROWS)
%
%   HEADER is a cell array of column names.  ROWS is a cell array with one
%   row per record; numeric entries are written with %.6g and character
%   entries verbatim.  Written with fopen and fprintf so that nothing here
%   depends on a toolbox, and so that the files stay diffable in git - the
%   PRD asks for results/ to be committed as CSV.

d = fileparts(path);
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end

fid = fopen(path, 'w');
if fid < 0
    error('write_csv:open', 'Could not open %s for writing.', path);
end

fprintf(fid, '%s', header{1});
fprintf(fid, ',%s', header{2:end});
fprintf(fid, '\n');

for i = 1:size(rows, 1)
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if j > 1, fprintf(fid, ','); end
        if ischar(v)
            fprintf(fid, '%s', v);
        else
            fprintf(fid, '%.6g', v);
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('  results -> %s (%d rows)\n', path, size(rows, 1));
end
