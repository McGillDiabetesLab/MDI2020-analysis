function plotTable(ui, cellTable, varargin)
YPad = 0.05;
XPad = 0.05;
YOffset = 0.05;
XOffset = 0.05;

TitleLength = 0.55;
XTitleLength = 0.35;

TableFont = 12;
XTitleFont = 13;
TitleFont = 15;
for nVar = 1:2:length(varargin)
    switch (lower(varargin{nVar}))
        case 'ypad'
            YPad = varargin{nVar+1};
        case 'xpad'
            XPad = varargin{nVar+1};
        case 'yoffset'
            YOffset = varargin{nVar+1};
        case 'xoffset'
            XOffset = varargin{nVar+1};
        case 'titlelength'
            TitleLength = varargin{nVar+1};
        case 'xtitlelength'
            XTitleLength = varargin{nVar+1};
        case 'titlefont'
            TitleFont = varargin{nVar+1};
        case 'xtitlefont'
            XTitleFont = varargin{nVar+1};
        case 'tablefont'
            TableFont = varargin{nVar+1};
        otherwise
            error('EstUtils::plotTable:WrongInput Unkown option %s', varargin{nVar});
    end
end

nbrElement = size(cellTable, 1) + 1;
YSize = (1 - (2 * YOffset + (nbrElement - 1) * YPad)) / nbrElement;

for lvl = 1:nbrElement - 1
    if isempty(cellTable{lvl, 1})
        continue;
    end
    
    if strcmpi(cellTable{lvl, 1}, 'title')
        nbrOfColumns = size(cellTable{lvl, 2}, 1);
        XSize = (1 - (2 * XOffset + nbrOfColumns * XPad + XTitleLength)) / nbrOfColumns;
        for colIdx = 1:nbrOfColumns
            temp = uicontrol(ui, 'Style', 'text');
            temp.Units = 'normalized';
            temp.String = cellTable{lvl, 2}(colIdx, :);
            temp.FontWeight = 'bold';
            temp.FontSize = XTitleFont;
            temp.Position = [XOffset + XTitleLength + colIdx * XPad + (colIdx - 1) * XSize, YOffset + (lvl - 1) * (YSize + YPad), XSize, YSize];
        end        
    else
        temp = uicontrol(ui, 'Style', 'text');
        temp.Units = 'normalized';
        temp.String = cellTable{lvl, 1};
        temp.FontSize = XTitleFont;
        temp.FontWeight = 'bold';
        temp.Position = [XOffset, YOffset + (lvl - 1) * (YSize + YPad), XTitleLength, YSize];
        
        nbrOfColumns = size(cellTable{lvl, 2}, 1);
        XSize = (1 - (2 * XOffset + nbrOfColumns * XPad + XTitleLength)) / nbrOfColumns;
        for colIdx = 1:nbrOfColumns
            temp = uicontrol(ui, 'Style', 'edit');
            temp.Units = 'normalized';
            temp.String = cellTable{lvl, 2}(colIdx, :);
            temp.FontSize = TableFont;
            temp.Position = [XOffset + XTitleLength + colIdx * XPad + (colIdx - 1) * XSize, YOffset + (lvl - 1) * (YSize + YPad), XSize, YSize];
        end
    end
end

lvl = lvl + 1;
temp = uicontrol(ui, 'Style', 'text');
temp.Units = 'normalized';
temp.String = ui.Tag;
temp.FontSize = TitleFont;
temp.FontWeight = 'bold';
temp.BackgroundColor = [0.3, 0.3, 0.3];
temp.ForegroundColor = [0.9, 0.9, 0.9];
temp.Position = [(1 - TitleLength) / 2, YOffset + (lvl - 1) * (YSize + YPad), TitleLength, YSize];
end
