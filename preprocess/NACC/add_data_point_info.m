function mri2 = add_data_point_info(mri1)
% add diagnosis, MMSE
mri1 = renamevars(mri1, 'years', 'years_from_bl');
ref_datetime = datetime(2000,1,1);
years1 = years(mri1.scandt - ref_datetime);
mri1 = addvars(mri1, years1, 'after', 'naccid');
mri1 = renamevars(mri1, 'years1', 'years');
mri1 = renamevars(mri1, 'naccid', 'RID');


clinical = readtable('./data/nacc_clinical.csv',VariableNamingRule='preserve');
inds = ismember(clinical.NACCID, unique(mri1.RID));
clinical = clinical(inds, {'NACCID','NACCUDSD','NACCMMSE','CDRSUM','VISITYR','VISITMO','VISITDAY'});
visit_dt = datetime(clinical.VISITYR, clinical.VISITMO, clinical.VISITDAY);
years1 = years(visit_dt - ref_datetime);
clinical = addvars(clinical, years1, 'after', 'NACCID');
clinical = renamevars(clinical, 'years1', 'years');
clinical = renamevars(clinical, 'NACCID', 'RID');
clinical = removevars(clinical, {'VISITYR','VISITMO','VISITDAY'});

% convert RID to numbers
mri1.RID = cellfun(@(x) str2num(x(5:end)), mri1.RID);
clinical.RID = cellfun(@(x) str2num(x(5:end)), clinical.RID);

% join the 2 tables
mri2 = join_by_RID_years(mri1, clinical, 'find_right_for_each_left');
mri2 = renamevars(mri2, {'NACCUDSD','NACCMMSE'}, {'diagnosis','MMSE'});

inds = mri2.MMSE == -4 | mri2.MMSE == 88 | ...
    mri2.MMSE == 95 | mri2.MMSE == 96 | ...
    mri2.MMSE == 97 | mri2.MMSE == 98;
mri2{inds, 'MMSE'} = NaN;

mri2.diagnosis = convert_NACCUDSD_to_dx(mri2.diagnosis);

mri2 = removevars(mri2, {'scandt','years'});
mri2 = renamevars(mri2, 'years_from_bl', 'years');
end
