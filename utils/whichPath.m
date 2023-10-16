function mainPath=whichPath()
pcID=getenv('COMPUTERNAME');
switch pcID
  case {'SEIDEMANN2PANAL'} %SEA
    boxDirStr='pktan\';
    dataDirStr='PK\';
  case {'LA-CPSA10646WD'} %ARC RM4
    boxDirStr='esexpt\';
    dataDirStr='PK\';
  case {'LA-CPSA07019WD'} %ARC RM5
    boxDirStr='esexpt\';
    dataDirStr='PK\';
end
mainPath.box=['C:\Users\' boxDirStr 'Box\'];
mainPath.data=['Y:\Users\' dataDirStr];
end
