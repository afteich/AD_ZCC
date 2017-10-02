
% This script analyzes ARACNe/MARINa data to filter MRs based on their
% support for genes in a given ontology category.  This script was used to
% generate a list of MRs that support synaptic function.  However, this
% basic outline could be used identify MRs that support any other ontology
% group.  In order to use this script, one needs the following; 1) A
% normalized gene expression profile for all samples in the analysis, 2) A
% list of MRs from MARINa that you are interested in, and 3) The 
% interactome that dictates the targets of each MR as well as the MI and 
% Spearman's correlation between the MR and each target.

% NOTE - xlswrite will NOT work on a Mac - Excel currently does not offer
% support for this function on Macs.  So you have to save your matricies
% into .txt files using extra code if you run this on a Mac (see below).

% First we load some data
disp('Importing dataset')
[data textdata] = importfile('Alzheimer_HGU133plus2_affy_NDAD.exp');
% Importfile is a homemade function included with this program; it spits
% out a data file with numbers and a textdata file with character strings.
% Both matricies are the size of the entire original matrix.  The string
% matrix has blanks where there are numbers.  The data matrix has a
% numerical value for each spot (it makes up a number where there are
% characters).
probes = textdata(2:end,1);
gene_names = textdata(2:end,2);
samples = textdata(1,3:end);
clear textdata
disp('Reading network ...')
adjfile = ['Alzheimer_HGU133plus2_Affy_NDAD_bs100_p9e-10_5col_bonf.txt'];
[tf tg tmp] = textread(adjfile,'%s%s%[^\n]','delimiter','\t');
corr_and_pvalue = dlmread(adjfile,'\t',0,3);
pvalue = dlmread(adjfile,'\t',0,4);
TF_Names = textread('TF_Names.txt','%s');
TG_Names = textread('TG_Names.txt','%s');

% Load your syanptic target list (we curated this list from Ingenuity)
Restricted_Synaptic_Targets = textread('Restricted_Synapse_List.txt','%s');

% Load your list of candidate MRs the first column has probe number.  The
% next six columns represent EC, Hipp, MTG, PC, SFG, and VCX; a 1 is
% present if the MR for that row is an MR in that region, a 0 if not
[AB_a,AB_b,AB_c,AB_d,AB_e,AB_f,AB_g] = textread('AllBrain_Control_vs_Affected_MRcandidateProbes.txt','%s%f%f%f%f%f%f');
MRwRegions = cell(length(AB_a),7);
MRwRegions(:,1) = AB_a;

% Load numeric values into cell array
for i = 1:length(AB_a)
    MRwRegions{i,2} = AB_b(i);
    MRwRegions{i,3} = AB_c(i);
    MRwRegions{i,4} = AB_d(i);
    MRwRegions{i,5} = AB_e(i);
    MRwRegions{i,6} = AB_f(i);
    MRwRegions{i,7} = AB_g(i);
end



% Get a list of all the variations of my synaptic targets
% Because there are multiple probes of a given gene, my list of
% Syn_Probe_List and Syn_Gene_List will be longer than my initial list
% "Restricted_Synaptic_Targets".  I will simultaneously generate a "data" matrix
% that consists only of synaptic genes, and the probe and name list will
% correspond to the rows of the data_Syn matrix
index = 1;
for i = 1:length(Restricted_Synaptic_Targets)
    a = strcmp(Restricted_Synaptic_Targets(i),gene_names);
    b = find(a>0);
    for j = 1:length(b)
        Syn_Probe_List(index,1) = probes(b(j));
        Syn_Gene_List(index,1) = gene_names(b(j));
        data_Syn(index,:) = data(b(j),:);
        index = index+1;
    end
end

% The "data" file has the following column indices (rows are genes):
% EC Affected values are indices 1 to 9
EC_Aff = [1:9];
% EC Control values are indices 10 to 20
EC_Ctl = [10:20];
% EC NDAD values are indices 21 to 26; we aren't using these
% Hippocampus Affected values are indices 27 to 36
HIPP_Aff = [27:36];
% Hippocampus Control values are indices 37 to 49
HIPP_Ctl = [37:49];
% Hippocampus NDAD values are indices 50 to 55; we aren't using these
% MTG affected values are indices 56 to 71
MTG_Aff = [56:71];
% MTG Control values are indices 72 to 83
MTG_Ctl = [72:83];
% MTG NDAD values are indices 84 to 89; we aren't using these
% PC (posterior cingulate) affected values are indices 90 to 98
PC_Aff = [90:98];
% PC control values are indices 99 to 111
PC_Ctl = [99:111];
% PC NDAD values are indices 112 to 116; we aren't using these
% SFG affected values are indices 117 to 139
SFG_Aff = [117:139];
% SFG control values are indices 140 to 149
SFG_Ctl = [140:149];
% SFG NDAD values are indices 150 to 155; we aren't using these
% VCX affected values are indices 156 to 174
VCX_Aff = [156:174];
% VCX control values are indices 175 to 186
VCX_Ctl = [175:186];
% VCX NDAD values are indices 187 to 191; we aren't using these


% THE COLUMNS OF THE NEXT MATRIX (AllBrain_TFwithTargets_CONvsAFF) are as
% follows: probe number of TF, TF name, probe number of target, target
% name, the correlation of the target with the TF, the sign of the correlation,
% and the difference in expression of that target between AD and control
% (only differences that pass an FDR threshold are included - see below).


% First generate your matrix for Control vs. Affected
index = 1;
for i = 1:length(MRwRegions) % get the number of rows, or TFs
    
    
    % Generate a matrix of gene expression differences specific for the
    % regions where this MR is ranked as significant (this is found in the
    % MRwRegions matrix that has been loaded).
    AllBrain_Control = data_Syn(:,nonzeros([EC_Ctl.*MRwRegions{i,2},HIPP_Ctl.*MRwRegions{i,3},MTG_Ctl.*MRwRegions{i,4},PC_Ctl.*MRwRegions{i,5},SFG_Ctl.*MRwRegions{i,6},VCX_Ctl.*MRwRegions{i,7}])');
    [xdim_control,ydim_control] = size(AllBrain_Control);

    AllBrain_Affected = data_Syn(:,nonzeros([EC_Aff.*MRwRegions{i,2},HIPP_Aff.*MRwRegions{i,3},MTG_Aff.*MRwRegions{i,4},PC_Aff.*MRwRegions{i,5},SFG_Aff.*MRwRegions{i,6},VCX_Aff.*MRwRegions{i,7}])');
    [xdim_affected,ydim_affected] = size(AllBrain_Affected);

    % Sum and average
    AllBrain_Control = sum(AllBrain_Control,2)/ydim_control; 
    AllBrain_Affected = sum(AllBrain_Affected,2)/ydim_affected;

    % Normalize by control value
    AllBrain_ControlvsAffected = (AllBrain_Affected - AllBrain_Control)./AllBrain_Control;

    % We want to know if the difference for each gene is statistically significant.  
    % If not, then we should not count the difference.  So we'll do a
    % t-test on each gene (with FDR correction).
    % We will then generate a vector with 1's for a gene if FDR<0.05 and 0
    % otherwise.  We'll then multiply this vector of 1's and 0's by the 
    % vector of differences for each comparison, and we'll proceed from
    % there with only statistically significant differences.
    FDR_threshold = 0.05;

    for j = 1:length(data_Syn) % data_Syn, AllBrain_Control, AllBrain_Affected, and AllBrain_ControlvsAffected all have same length
        temp_AllBrain_control = data_Syn(j,nonzeros([EC_Ctl.*MRwRegions{i,2},HIPP_Ctl.*MRwRegions{i,3},MTG_Ctl.*MRwRegions{i,4},PC_Ctl.*MRwRegions{i,5},SFG_Ctl.*MRwRegions{i,6},VCX_Ctl.*MRwRegions{i,7}])');
        temp_AllBrain_affected = data_Syn(j,nonzeros([EC_Aff.*MRwRegions{i,2},HIPP_Aff.*MRwRegions{i,3},MTG_Aff.*MRwRegions{i,4},PC_Aff.*MRwRegions{i,5},SFG_Aff.*MRwRegions{i,6},VCX_Aff.*MRwRegions{i,7}])');
        [temp_a,temp_b] = ttest2(temp_AllBrain_control,temp_AllBrain_affected);
        Significance_AllBrain_ControlvsAffected(j,1:2) = [temp_a,temp_b];

    end

    Significance_AllBrain_ControlvsAffected(:,3) = mafdr(Significance_AllBrain_ControlvsAffected(:,2));
    for j = 1:length(Significance_AllBrain_ControlvsAffected)
        if  Significance_AllBrain_ControlvsAffected(j,3) < FDR_threshold
            Significance_AllBrain_ControlvsAffected(j,4) = 1;
        else
            Significance_AllBrain_ControlvsAffected(i,4) = 0;
        end
    end


    AllBrain_ControlvsAffected = AllBrain_ControlvsAffected.*Significance_AllBrain_ControlvsAffected(:,4); % filter out values less than FDR threshold


    
    
    % Now use this matrix of expression values specific for your MR (AllBrain_ControlvsAffected) to
    % figure out the rest of the entries
    
    e = strcmp(MRwRegions(i,1),tf);
    f = find(e>0);
    clear tg_temp;
    tg_temp = tg(f);
    corr_pvalue_temp =  corr_and_pvalue(f,1);
    for j = 1:length(tg_temp)
       if (sum(strcmp(tg_temp(j),Syn_Probe_List))>0) % This means the probe is a synaptic target
           AllBrain_TFwithTargets_CONvsAFF(index,1) = MRwRegions(i,1);
           AllBrain_TFwithTargets_CONvsAFF(index,3) = tg_temp(j);
           AllBrain_TFwithTargets_CONvsAFF{index,5} = corr_pvalue_temp(j);
           AllBrain_TFwithTargets_CONvsAFF{index,6} = sign(corr_pvalue_temp(j));
           g = strcmp(tg_temp(j),Syn_Probe_List);
           h = find(g>0);
           AllBrain_TFwithTargets_CONvsAFF{index,7} = AllBrain_ControlvsAffected(h); % this is a modified data_Syn, which has the same order as Syn_Probe_List
           index = index+1;
       end
    end
end


for i = 1:length(AllBrain_TFwithTargets_CONvsAFF)
    e = strcmp(AllBrain_TFwithTargets_CONvsAFF(i,1),probes);
    f = find(e>0);
    AllBrain_TFwithTargets_CONvsAFF(i,2) = gene_names(f);
    e = strcmp(AllBrain_TFwithTargets_CONvsAFF(i,3),probes);
    f = find(e>0);
    AllBrain_TFwithTargets_CONvsAFF(i,4) = gene_names(f);
end


% Save your data (this only works for Windows)
%xlswrite('AllBrain_TFwithTargets_CONvsAFF.xlsx', AllBrain_TFwithTargets_CONvsAFF);


% If you're on a Mac, use the code below to save your data
[nrows,ncols]= size(AllBrain_TFwithTargets_CONvsAFF);
filename = 'AllBrain_TFwithTargets_CONvsAFF.txt';
fid = fopen(filename, 'w+');
for row=1:nrows
    fprintf(fid, '%s\t %s\t %s\t %s\t %f\t %f\t %f\n', AllBrain_TFwithTargets_CONvsAFF{row,:});
end
fclose(fid);



% Remember: the columns of AllBrain_TFwithTargets_CONvsAFF are as
% follows: probe number of TF, TF name, probe number of target, target
% name, the correlation of the target with the TF, the sign of the correlation,
% and the difference in expression of that target between AD and control
% (only differences that pass an FDR threshold are included - see below).


% Now start generating the results matrix SUMMARY_AllBrain_TFwithTargets_CONvsAFF
GENE_CHANGE_SUM = 0;
GENE_COUNT = 0;
GENE_CHANGE_SUM_WITH_CORR = 0;
GENE_CHANGE_SUM_WITH_SIGN = 0;
GENE_SIGN_DIFFERENCE = 0;
index = 1;
for i = 1:length(AllBrain_TFwithTargets_CONvsAFF)
    
    GENE_CHANGE_SUM = GENE_CHANGE_SUM + AllBrain_TFwithTargets_CONvsAFF{i,7};
    GENE_CHANGE_SUM_WITH_CORR = GENE_CHANGE_SUM_WITH_CORR + AllBrain_TFwithTargets_CONvsAFF{i,5}*AllBrain_TFwithTargets_CONvsAFF{i,7};
    GENE_CHANGE_SUM_WITH_SIGN = GENE_CHANGE_SUM_WITH_SIGN + AllBrain_TFwithTargets_CONvsAFF{i,6}*AllBrain_TFwithTargets_CONvsAFF{i,7};
    GENE_SIGN_DIFFERENCE = GENE_SIGN_DIFFERENCE + AllBrain_TFwithTargets_CONvsAFF{i,6};
    GENE_COUNT = GENE_COUNT + 1;
    
    if i == length(AllBrain_TFwithTargets_CONvsAFF) % you're at end of list 
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF(index,1)=AllBrain_TFwithTargets_CONvsAFF(i,1);
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF(index,2)=AllBrain_TFwithTargets_CONvsAFF(i,2);
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,3} = GENE_CHANGE_SUM;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,4} = GENE_COUNT;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,5} = GENE_CHANGE_SUM_WITH_CORR;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,6} = GENE_CHANGE_SUM_WITH_SIGN;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,7} = GENE_SIGN_DIFFERENCE;
        
        % All done!
        % This contingency comes first because otherwise the loop trips at
        % the end 
        
        
    elseif strcmp(AllBrain_TFwithTargets_CONvsAFF{i,1},AllBrain_TFwithTargets_CONvsAFF{i+1,1})==0 % the next line is a new MR
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF(index,1)=AllBrain_TFwithTargets_CONvsAFF(i,1);
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF(index,2)=AllBrain_TFwithTargets_CONvsAFF(i,2);
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,3} = GENE_CHANGE_SUM;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,4} = GENE_COUNT;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,5} = GENE_CHANGE_SUM_WITH_CORR;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,6} = GENE_CHANGE_SUM_WITH_SIGN;
        SUMMARY_AllBrain_TFwithTargets_CONvsAFF{index,7} = GENE_SIGN_DIFFERENCE;
        
        GENE_CHANGE_SUM = 0;
        GENE_COUNT = 0;
        GENE_CHANGE_SUM_WITH_CORR = 0;
        GENE_CHANGE_SUM_WITH_SIGN = 0;
        GENE_SIGN_DIFFERENCE = 0;
        index = index + 1;
        
        
        
    end
end
        

% This adds on the end a total number of genes in each interactome and
% performs Fisher exact test
for i = 1:length(SUMMARY_AllBrain_TFwithTargets_CONvsAFF(:,1))
    temp_TF = SUMMARY_AllBrain_TFwithTargets_CONvsAFF(i,2);
    temp_index1 = strcmp(TF_Names(:),temp_TF);
    temp_index2 = find(temp_index1>0);
    temp_TG = TG_Names(temp_index2);
    TG_unique = unique(temp_TG);
    TG_number = length(TG_unique);
    SUMMARY_AllBrain_TFwithTargets_CONvsAFF{i,8} = TG_number;
    
    TF_syn_genes = 0;
    for j = 1:length(TG_unique)
        if (sum(strcmp(TG_unique(j),Restricted_Synaptic_Targets))>0) % This means the probe is a synaptic target
            TF_syn_genes = TF_syn_genes + 1;
        end
    end
    SUMMARY_AllBrain_TFwithTargets_CONvsAFF{i,9} = TF_syn_genes;
    
    
    total_syn = length(unique(Syn_Gene_List));
    total_genes = length(unique(gene_names)) - total_syn;
    TF_total_genes = TG_number - TF_syn_genes;
    
    x = table([total_syn;TF_syn_genes],[total_genes;TF_total_genes],'VariableNames',{'Syn','NoSyn'},'RowNames',{'AllGenes','TFGenes'});
    % [h,p,stats] = fishertest(x, 'Tail', 'left'); % Fisher Exact Test
    p = hygecdf(TF_syn_genes,length(unique(gene_names)),total_syn,TG_number,'upper'); % Hypergeometric Test
    SUMMARY_AllBrain_TFwithTargets_CONvsAFF{i,10} = p;
    
end




% Save your data (Use the xlswrite function only for Windows)
%xlswrite('SUMMARY_AllBrain_TFwithTargets_CONvsAFF.xlsx', SUMMARY_AllBrain_TFwithTargets_CONvsAFF);



% If you're on a Mac, use the below code to save your data
[nrows,ncols]= size(SUMMARY_AllBrain_TFwithTargets_CONvsAFF);
filename = 'SUMMARY_AllBrain_TFwithTargets_CONvsAFF.txt';
fid = fopen(filename, 'w+');
for row=1:nrows
    fprintf(fid, '%s\t %s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n', SUMMARY_AllBrain_TFwithTargets_CONvsAFF{row,:});
end
fclose(fid);









