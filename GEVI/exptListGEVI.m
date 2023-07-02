%% Biasing sessions with entry no, points to session folder/run. Typically run 0 contains orientation map, run >=1 contains biasing experiment
entryNo=1;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230407';
datastruct(entryNo).run='0';
datastruct(entryNo).series='Contrast';
datastruct(entryNo).modality='GEVI';

entryNo=2;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230407';
datastruct(entryNo).run='6';
datastruct(entryNo).series='Contrast';
datastruct(entryNo).modality='GEVI';

entryNo=3; % MERGE SUPERBLOCK
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230407';
datastruct(entryNo).run=[1 3 4 5 8]; % 4 4 10 20 20
datastruct(entryNo).series='Temporal';
datastruct(entryNo).modality='GEVI';

entryNo=4;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230426';
datastruct(entryNo).run='0';
datastruct(entryNo).series='Contrast';
datastruct(entryNo).modality='GEVI';

entryNo=5;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230426';
datastruct(entryNo).run='1';
datastruct(entryNo).series='Temporal';
datastruct(entryNo).modality='GEVI';

entryNo=6;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230502';
datastruct(entryNo).run='0';
datastruct(entryNo).series='Temporal';
datastruct(entryNo).modality='GEVI';

entryNo=7;
datastruct(entryNo).monkeyNo='28';
datastruct(entryNo).monkey='Chip';
datastruct(entryNo).date='20230502';
datastruct(entryNo).run='1';
datastruct(entryNo).series='Contrast';
datastruct(entryNo).modality='GEVI';

