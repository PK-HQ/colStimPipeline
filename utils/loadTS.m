function TS = loadTS(datastructEntry, dataType, baselineCondition)
    if baselineCondition
        TS = load(datastructEntry.(dataType), 'TS');
    else
        TS(1) = load(datastructEntry.(dataType), 'TS');
        if ~isempty(datastructEntry.baselineIntg)
            TS(2) = load(datastructEntry.baselineTS, 'TS');
        end
    end
end
