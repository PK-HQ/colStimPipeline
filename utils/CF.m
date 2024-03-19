function CF()
try
    % find all handles of axes (graphs)
    axh = findall(groot,'type','axes');
    % get handles of parent figures containing graphs
    fxh = get(axh,'parent');
    % close figures containg axes
    close(fxh{:});
catch
end
end