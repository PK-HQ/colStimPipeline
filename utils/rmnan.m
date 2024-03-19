function array=rmnan(array)
array=array(~isnan(array));
end