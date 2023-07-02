function [w]=suppressWarning()
w = warning('query','last');
id=w.identifier;
warning('off',id);
end