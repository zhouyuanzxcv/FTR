function position = get_figure_position(idx)

if nargin < 1
    idx = 1;
end

if idx == 1 % 1:1
    position = [599 422 272 250]; % 1681 × 1622 for 3C
elseif idx == 2 % 1:2
    position = [264,404,620,277];
elseif idx == 3 % 1:3
    position = [264,404,968,277];
elseif strcmp(idx,'fig5')
    position = [827   194   251   798];
elseif strcmp(idx, 'fig51')
    position = [367   264   251   604];
elseif strcmp(idx,'figS12')
    position = [636,302,408,432];
elseif idx == 4
    position = [680   802   478   190];
elseif idx == 5
    position = [599   483   219   195];
elseif idx == 6
    position = [691   679   842   265];
elseif idx == 7
    position = [599   475   236   197];
elseif idx == 8
    position = [785   391   282   125];
elseif idx == 9
    position = [450   353   259   243];
elseif idx == 10
    position = [1006  410  335  276];
elseif idx == 11
    position = [605   567   287   116];
elseif idx == 12
    position = [269   297   819   540];
end

end