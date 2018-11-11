function latex = latexVarianceTable(input)
% An easy to use function that generates a LaTeX table from a given MATLAB
% input struct containing numeric values. The LaTeX code is printed in the
% command window for quick copy&paste and given back as a cell array.
%
% Author:       Eli Duenisch
% Contributor:  Jasmine Hao
% Date:         Nov 9, 2018
% License:      This code is licensed using BSD 2 to maximize your freedom of using it :)
% ----------------------------------------------------------------------------------
%  Copyright (c) 2016, Eli Duenisch
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%  
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
%  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ----------------------------------------------------------------------------------
%
% Input:
% input    struct containing your data and optional fields (details described below)
%
% Output:
% latex    cell array containing LaTex code
%
% Example and explanation of the input struct fields:
%
% % numeric values you want to tabulate:
% % this field has to be a matrix or MATLAB table datatype
% % missing values have to be NaN
% % in this example we use an array
% input.data = [1.12345 2.12345 3.12345; ...
%               4.12345 5.12345 6.12345; ...
%               7.12345 NaN 9.12345; ...
%               10.12345 11.12345 12.12345];
% input.variance 
% % Optional fields (if not set default values will be used):
%
% % Set the position of the table in the LaTex document using h, t, p, b, H or !
% input.tablePositioning = 'h';
%
% % Set column labels (use empty string for no label):
% input.tableColLabels = {'col1','col2','col3'};
% % Set row labels (use empty string for no label):
% input.tableRowLabels = {'row1','row2','','row4'};
%
% % Switch transposing/pivoting your table:
% input.transposeTable = 0;
%
% % Determine whether input.dataFormat is applied column or row based:
% input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'column' is used
%
% % Formatting-string to set the precision of the table values:
% % For using different formats in different rows use a cell array like
% % {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% % where myFormatString_ are formatting-strings and numberOfValues_ are the
% % number of table columns or rows that the preceding formatting-string applies.
% % Please make sure the sum of numberOfValues_ matches the number of columns or
% % rows in input.tableData!
% %
% % input.dataFormat = {'%.3f'}; % uses three digit precision floating point for all data values
% input.dataFormat = {'%.3f',2,'%.1f',1}; % three digits precision for first two columns, one digit for the last
%
% % Define how NaN values in input.tableData should be printed in the LaTex table:
% input.dataNanString = '-';
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% input.tableColumnAlignment = 'c';
%
% % Switch table borders on/off:
% input.tableBorders = 1;
%
% % Switch table booktabs on/off:
% input.booktabs = 1;
%
% % LaTex table caption:
% input.tableCaption = 'MyTableCaption';
%
% % LaTex table label:
% input.tableLabel = 'MyTableLabel';
%
% % Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
%
% % % Now call the function to generate LaTex code:
% latex = latexTable(input);

%%%%%%%%%%%%%%%%%%%%%%%%%% Default settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These settings are used if the corresponding optional inputs are not given.
%
% Placement of the table in LaTex document
if isfield(input,'tablePlacement') && (length(input.tablePlacement)>0)
    input.tablePlacement = ['[',input.tablePlacement,']'];
else
    input.tablePlacement = '';
end
% Pivoting of the input data switched off per default:
if ~isfield(input,'transposeTable'),input.transposeTable = 0;end
% Default mode for applying input.tableDataFormat:
if ~isfield(input,'dataFormatMode'),input.dataFormatMode = 'column';end
% Sets the default display format of numeric values in the LaTeX table to '%.4f'
% (4 digits floating point precision).
if ~isfield(input,'dataFormat'),input.dataFormat = {'%.4f'};end
% Define what should happen with NaN values in input.tableData:
if ~isfield(input,'dataNanString'),input.dataNanString = '-';end
% Specify the alignment of the columns:
% 'l' for left-justified, 'c' for centered, 'r' for right-justified
if ~isfield(input,'tableColumnAlignment'),input.tableColumnAlignment = 'r';end
% Specify whether the table has borders:
% 0 for no borders, 1 for borders
if ~isfield(input,'tableBorders'),input.tableBorders = 0;
% Specify whether to use booktabs formatting or regular table formatting:
if ~isfield(input,'booktabs')
    input.booktabs = 1;
else
    if input.booktabs
        input.tableBorders = 0;
    end
end
% Other optional fields:
if ~isfield(input,'tableCaption'),input.tableCaption = 'MyTableCaption';end
if ~isfield(input,'tableLabel'),input.tableLabel = 'MyTableLabel';end
if ~isfield(input,'makeCompleteLatexDocument'),input.makeCompleteLatexDocument = 0;end
if ~size(input.data) == size(input.variance), input.variance = zeros(size(input.data)); end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process table datatype
if isa(input.data,'table')
  if(~isempty(input.data.Properties.RowNames))
    input.tableRowLabels = input.data.Properties.RowNames';
  end
  if(~isempty(input.data.Properties.VariableNames))
    input.tableColLabels = input.data.Properties.VariableNames';
  end
    input.data = table2array(input.data);
end

% get size of data
numberDataRows = size(input.data,1);
numberDataCols = size(input.data,2);

% obtain cell array for the table data and labels
colLabelsExist = isfield(input,'tableColLabels');
rowLabelsExist = isfield(input,'tableRowLabels');
cellSize = [2 *numberDataRows+colLabelsExist,numberDataCols+rowLabelsExist];
C = cell(cellSize);
C((2 * (0:numberDataRows-1) + 1+colLabelsExist),1+colLabelsExist:end ) = num2cell(input.data);
C((2 * (1:numberDataRows) +colLabelsExist),1+colLabelsExist:end ) = cellfun(@(x) {sprintf("(%.4f)",x)}, num2cell(input.variance));

if rowLabelsExist
    C((2 * (0:numberDataRows-1) + 1+colLabelsExist),1)=input.tableRowLabels';
end
if colLabelsExist
    C(1,1+rowLabelsExist:end)=input.tableColLabels;
end

% obtain cell array for the format
lengthDataFormat = length(input.dataFormat);
% if lengthDataFormat==1
%     tmp = repmat(input.dataFormat(1),2 * numberDataRows,numberDataCols);
% else
%     dataFormatList={};
%     for i=1:2:lengthDataFormat
%         dataFormatList(end+1:end+input.dataFormat{i+1},1) = repmat(input.dataFormat(i),input.dataFormat{i+1},1);
%     end
%     if strcmp(input.dataFormatMode,'column')
%         tmp = repmat(dataFormatList',numberDataRows,1);
%     end
%     if strcmp(input.dataFormatMode,'row')
%         tmp = repmat(dataFormatList,1,numberDataCols);
%     end
% end
% if ~isequal(size(tmp),size(input.data))
%     error(['Please check your values in input.dataFormat:'...
%         'The sum of the numbers of fields must match the number of columns OR rows '...
%         '(depending on input.dataFormatMode)!']);
% end
dataFormatArray = cell(cellSize);
dataFormatArray(1+colLabelsExist:end,1+rowLabelsExist:end) = input.dataFormat;

% transpose table (if this switched on)
if input.transposeTable
    C = C';
    dataFormatArray = dataFormatArray';
end

% make table header lines:
hLine = '\hline';
if input.tableBorders
    header = ['\begin{tabular}','{l|',repmat([input.tableColumnAlignment],1,size(C,2)-1),'}'];
else
    header = ['\begin{tabular}','{l|',repmat(input.tableColumnAlignment,1,size(C,2)-1),'}'];
end
latex = {['\begin{table}[!htbp]',input.tablePlacement];'\centering';['\caption{',input.tableCaption,'}'];header};

% generate table
if input.booktabs
    latex(end+1) = {hLine};
end    

for i=1:size(C,1)
    if i==2 && input.booktabs
        latex(end+1) = {hLine};
    end
    if input.tableBorders
        latex(end+1) = {hLine};
    end
    rowStr = '';
    for j=1:size(C,2)
        dataValue = C{i,j};
        if iscell(dataValue) || isstring(dataValue)
          dataValue = dataValue{:};
        elseif isnan(dataValue)
          dataValue = input.dataNanString;
        elseif isnumeric(dataValue)
          dataValue = num2str(dataValue,dataFormatArray{i,j});
        end
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    latex(end+1) = {[rowStr,' \\']};
end

if input.booktabs
    latex(end+1) = {hLine};
end   


% make footer lines for table:
tableFooter = {'\end{tabular}'; ...
    ['\label{table:',input.tableLabel,'}'];'\end{table}'};
if input.tableBorders
    latex = [latex;{hLine};tableFooter];
else
    latex = [latex;tableFooter];
end

% add code if a complete latex document should be created:
if input.makeCompleteLatexDocument
    % document header
    latexHeader = {'\documentclass[letterpaper,10pt]{article}'};
    if input.booktabs
        latexHeader(end+1) = {'\usepackage{booktabs}'};
    end 
    latexHeader(end+1) = {'\begin{document}'};
    % document footer
    latexFooter = {'\end{document}'};
    latex = [latexHeader';latex;latexFooter];
end

% print latex code to console:
disp(char(latex));

end
