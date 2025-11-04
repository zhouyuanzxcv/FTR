function num_bioms = get_num_bioms_from_ver(biom_ver, dataset_name)
if strcmp(dataset_name, 'NACC')
    if biom_ver == "HM"
        num_bioms = 32;
    elseif biom_ver == "HS"
        num_bioms = 64;
    elseif biom_ver == "LM"
        num_bioms = 7;
    elseif biom_ver == "LS"
        num_bioms = 14;
    end
else
    if biom_ver == "HM"
        num_bioms = 41;
    elseif biom_ver == "HS"
        num_bioms = 82;
    elseif biom_ver == "LM"
        num_bioms = 13;
    elseif biom_ver == "LS"
        num_bioms = 26;
    end
end

end
