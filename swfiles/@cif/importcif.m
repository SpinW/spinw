function [cifdat, source, isfile] = importcif(~, dataStr)
% imports .cif data from file, web or string

% if exist(dataStr,'file') == 2
%     % get the name of the file
%     fid = fopen(dataStr);
%     source = fopen(fid);
%     fclose(fid);
%     
%     cifStr = regexp(fileread(dataStr), ['(?:' sprintf('\n') ')+'], 'split');
%     isfile = true;
%     
% elseif numel(dataStr) < 200
%     % try to load it from the web
%     try
%         %cifStr = char(webread(dataStr)');
%         cifStr = urlread(dataStr);
%     catch
%         error('cif:importcif:WrongInput','The requested data cannot be found!')
%     end
%     cifStr = strsplit(cifStr,'\n');
%     
%     source = dataStr;
%     isfile = false;
% else
%     % take it as a string storing the .cif file content
%     cifStr = strsplit(char(dataStr(:))','\n');
%     %cifStr = mat2cell(cifStr,ones(size(cifStr,1),1));
%     source = '';
%     isfile = false;
% end

[cifStr,info] = ndbase.source(dataStr);
isfile = info.isfile;
source = info.source;

% split the string
cifStr = strsplit(cifStr,'\n');

% unite broken lines
bLine = double(cellfun(@(x)numel(x)>0 && x(1)==';',cifStr));

% first and last line of the series of broken lines
firstL = true;
firstBL = [];
lastBL  = [];
for ii = 1:numel(bLine)
    if bLine(ii)
        if firstL
            firstBL(end+1) = ii;
            firstL = false;
        else
            lastBL(end+1) = ii;
            firstL = true;
        end
    end
end

cifStr2 = {};
BL = false;
%idx = 1;
isfirst = true;

for ii = 1:numel(cifStr)
    if any(ii==firstBL)
        BL = true;
        isfirst = true;
    end
    
    if ~BL
        cifStr2{end+1} = cifStr{ii};
        isfirst = true;
    else
        if isfirst
            cifStr2{end} = [strtrim(cifStr2{end}) ' '' ' strtrim(cifStr{ii}(2:end))];
        else
            if all(ii~=lastBL)
                cifStr2{end} = [cifStr2{end} '\n'  cifStr{ii}];            
            else
                cifStr2{end} = [cifStr2{end} ''''];
            end
        end
        % end comment sign
        %cifStr2{end} = [cifStr2{end} char([32 26])  cifStr{ii}(2:end)]; %' ?'

        isfirst = false;
    end
    
    if any(ii==lastBL)
        BL = false;
        isfirst = true;
    end
end

cifStr = cifStr2;

strout = {};
for ii = 1:numel(cifStr)
    strout{ii} = parseline(cifStr{ii});
end

% remove comments
for ii = 1:numel(strout)
    comIdx = strcmp('comment',strout{ii}(1,:));
    strout{ii} = strout{ii}(:,~comIdx); %#ok<*AGROW>
end

% remove empty lines
emptyIdx = cellfun(@(x)isempty(x),strout);
strout(emptyIdx) = [];

LP = false;
lpVar = cell(0,0);
lpVal = cell(0,0);

cifdat.name = '';
cifdat.val  = '';

savedata = false;

for ii = 1:numel(strout) 
    strLine = strout{ii};
    switch strLine{1,1}
        case 'loop'
            if LP
                % save data
                savedata = true;
            else
                LP = true;
            end
        case 'variable'
            if (~isempty(lpVal)) || (size(strLine,2)>1)
                LP = false;
                % save data
                savedata = true;
            end
            
            if LP
                lpVar = [lpVar strLine(2,1)];
            elseif (size(strLine,2)> 1) && ismember(strLine(1,2),{'number' 'string'})
                cifdat(end+1).name = strLine{2,1};
                cifdat(end).val    = strLine{2,2};
                cifdat(end).type   = strLine{1,2};
            else
                cifdat(end+1).name = strLine{2,1};
                cifdat(end).val    = '';
                %warning('Wrong stuff comes after a variable in line: %d, %s',ii,strLine{2,1});
            end
            
        case 'number'
            if LP
                
                if isempty(lpVal)
                    lpVal  = [strLine(2,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                    lpType = [strLine(1,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                else
                    filled = ~cellfun(@(x)isempty(x),lpVal(end,:));
                    if all(filled)
                        lpVal(end+1,:) = [strLine(2,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                    else
                        emptyIdx = find(~filled);
                        lpVal(end,emptyIdx(1:numel(strLine(2,:)))) = strLine(2,:);
                    end
                end
            end
        case 'string'
            if LP
                
                if isempty(lpVal)
                    lpVal  = [strLine(2,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                    lpType = [strLine(1,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                else
                    filled = ~cellfun(@(x)isempty(x),lpVal(end,:));
                    if all(filled)
                        lpVal(end+1,:) = [strLine(2,:) cell(1,numel(lpVar)-numel(strLine(2,:)))];
                    else
                        emptyIdx = find(~filled);
                        lpVal(end,emptyIdx(1:numel(strLine(2,:)))) = strLine(2,:);
                        
                    end
                end
            end
            
        otherwise
            if ~isempty(lpVal)
                % save data
                savedata = true;
            end
            LP = false;
    end
    if savedata || (ii == numel(strout))
        for jj = 1:numel(lpVar)
            cifdat(end+1).name = lpVar{jj};
            if jj > size(lpVal,2)
                cifdat(end).type   = 'string';
                cifdat(end).val    = '';
            else
                cifdat(end).type   = lpType{1,jj};
                cifdat(end).val    = lpVal(:,jj);
                if strcmp(lpType{1,jj},'number')
                    cifdat(end).val = [cifdat(end).val{:}]';
                end
            end
        end
        lpVar = cell(0,0);
        lpVal = cell(0,0);
        
        savedata= false;
    end
end

cifdat = cifdat(2:end);

end



function strout = parseline(strin)
% parse cif file line
strout = cell(2,0);

while ~isempty(strin)
    % remove whitespace
    strin = strtrim(strin);
    
    if isempty(strin)
        return;
    end
    
    % COMMENT
    if strin(1) == '#'
        endcomIdx = numel(strin)+1;
        strout{1,end+1} = 'comment';
        strout{2,end}   = strin(2:(endcomIdx-1));
        strin = strin(endcomIdx:end);
        continue;
    end
    
    % VARIABLE
    if strin(1) == '_'
        whiteIdx    = find(isspace(strin),1,'first');
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        strout{1,end+1} = 'variable';
        strout{2,end} = strin(2:(whiteIdx-1));
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
    % LOOP
    if (numel(strin) >4) && strcmpi(strin(1:5),'loop_')
        strout{1,end+1} = 'loop';
        strout{2,end} = [];
        return;
    end
    
    % STRING
    if strin(1) == ''''
        strIdx = find(strin(2:end)=='''',1,'first')+1;
        if isempty(strIdx)
            error('importcif:WrongCif','String is not closed!');
        end
        
        strout{1,end+1} = 'string';
        strout{2,end} = strin(2:(strIdx-1));
        strin = strin((strIdx+1):end);
        continue;
    end
    
    % NUMBER
    if ismember(strin(1), [num2cell('0':'9') {'-' '.'}])
       
        bracket1Idx = find(strin(1:end)=='(',1,'first');
        bracket2Idx = find(strin(1:end)==')',1,'first');
        whiteIdx    = find(isspace(strin),1,'first');
        if isempty(bracket1Idx)
            bracket1Idx = numel(strin)+1;
        end
        if isempty(bracket2Idx)
            bracket2Idx = numel(strin)+1;
        end
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        
        if whiteIdx <= bracket1Idx
            strout{1,end+1} = 'number';
            numOut = str2double(strin(1:(whiteIdx-1)));
            if ~isempty(numOut)
                strout{2,end} = numOut;
            else
                strout{2,end} = 0;
            end
            strin = strin((whiteIdx+1):end);
            continue;
        else
            strout{1,end+1} = 'number';
            numOut = str2double(strin(1:(bracket1Idx-1)));
            if ~isempty(numOut)
                strout{2,end} = numOut;
            else
                strout{2,end} = 0;
            end
            strin = strin((bracket2Idx+1):end);
            continue;
            
        end
    end
    
    % EMPTY NUMBER
    if strin(1) == '?'
        whiteIdx    = find(isspace(strin),1,'first');
        strout{1,end+1} = 'number';
        strout{2,end} = NaN;
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
    % END COMMENT WITHOUT BEGIN
    if strcmp(strin(1),char(63)) %'?'
        strin = strin(2:end);
        continue;
    end
    
    % STRING
    if (isstrprop(strin(1), 'alpha')) || true
        whiteIdx    = find(isspace(strin),1,'first');
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        
        strout{1,end+1} = 'string';
        strout{2,end} = strin(1:(whiteIdx-1));
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
end

end
