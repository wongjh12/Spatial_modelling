*** things to note ****
check microbeNames and abundancelist if they match!
change all S_ to wtv value

% upload human model and get core rxns
load('Recon3E_O2Fixed.mat')
model = modelO2Fixed

% get core rxns for biomass obj
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model,90);
core1 = find(grRatio<0.5);

% get core rxns for lactate obj
model.c(find(model.c))=0;
model.c(find(ismember(model.rxns,"EX_lac_D(e)")))=1;
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model,90);
core2 = find(grRatio<0.5);

% union of both core rxns
core = union(core1,core2)
%%%%% saved as starthere.mat


%%%%%%%% Start run from here

% build individual models
load('starthere.mat')
data = readtable("DataS_.csv");
expdata = data(:, 3:width(data));
genes = cellstr(string(data.ENTREZID));
names = data.Properties.VariableNames(:,3:width(data)); %(3:19)

for i=1:width(expdata)
	levels = table2array(expdata(:,i));
	reaction_levels = gene_to_reaction_levels(model,genes, levels , @min, @(x,y)(x+y));
	reaction_levels(isnan(reaction_levels))=0;
	reactionWeights = max(reaction_levels)-reaction_levels;
	tissueModel = fastCoreWeighted_mod(core, model, reactionWeights);	
	savename = names{i}
	save(savename, 'tissueModel')
end

% build comm models
% some lb of exchange rxns need to be changed
load('communityworkspace.mat', 'lbchange')
load('communityworkspace.mat', 'excrxnchange')

for i=1:str2double(extractBetween(names{end},"_00","_")) % 7 comm models
	modellist={};
	for ii=1:length(names)
		if find(contains(names(ii),['S12_00'+string(i)]))
		modelx = load(string(names(ii)) +'.mat', 'tissueModel');
		modeladd = modelx.tissueModel;
		modellist=[modellist;modeladd];
		end
	end

	for iii=1:length(modellist)
		id = find(contains(names,['S12_00'+string(i)]));
		biomasstemp = [names{id(iii)},'_biomass[c]'];
		modellist{iii,1} = addMetabolite(modellist{iii,1},biomasstemp);
		x= find(ismember(modellist{iii,1}.mets,biomasstemp));
		y= find(ismember(modellist{iii,1}.rxns,'biomass_reaction'));
		modellist{iii,1}.S(x,y)=1;
		microbeNames(i,iii) = cellstr(names(id(iii)));
	end

	commModelname = sprintf('commModel%d',i);
	[commModel] = createMultipleSpeciesModel(modellist);
	abundances = ones(length(modellist),1);
	commModellist{i} = addMicrobeCommunityBiomass(commModel,  string(microbeNames(i,1:length(find(~cellfun(@isempty,microbeNames(i,:))))))', abundances)
	commModellist{i}.c=zeros(length(commModellist{i}.rxns),1);
	commModellist{i}.c(find(ismember(commModellist{i}.rxns,'communityBiomass')))=1;

	% remove indiv biomass and set the lb of indiv ex_biomass to 0
	for iv=1:length(modellist)
		commModellist{i}.lb(find(ismember(commModellist{i}.rxns,['EX_'+string(microbeNames{i,iv})+'_biomass[c]'])))=0;
		commModellist{i}.ub(find(ismember(commModellist{i}.rxns,['EX_'+string(microbeNames{i,iv})+'_biomass[c]'])))=0;
		commModellist{i}.lb(find(ismember(commModellist{i}.rxns,['model1_IEX_'+string(microbeNames{i,iv})+'_biomass[c]tr'])))=0;
		commModellist{i}.lb(find(ismember(commModellist{i}.rxns,['model2_IEX_'+string(microbeNames{i,iv})+'_biomass[c]tr'])))=0;
		commModellist{i}.lb(find(ismember(commModellist{i}.rxns,['model3_IEX_'+string(microbeNames{i,iv})+'_biomass[c]tr'])))=0;
	end
	commModellist{i}.lb(find(ismember(commModellist{i}.rxns,['EX_microbeBiomass[fe]'])))=0;

	for v=1:length(lbchange)
		commModellist{i}.lb(find(ismember(commModellist{i}.rxns,excrxnchange(v))))=lbchange(v);
	end

	commModellist{i}.c(commModellist{i}.c==1)=0;
	commModellist{i}.c(find(ismember(commModellist{i}.rxns,'EX_microbeBiomass[fe]')))=1;
	sol=optimizeCbModel(commModellist{i});
	commfluxlist{i}=sol.f;
end


%%%%%%%%%% ABUNDANCES calculate by adding up all the nuclei parts tgt thn taking eg tumor / (cd3+tumor)
% format: following microbeNames
abundancelist = table2array(readtable("abundanceS_.csv"));

for i = 1:length(commModellist)
	ind = find(commModellist{i}.S(:,find(ismember(commModellist{i}.rxns,'communityBiomass'))))
	ind= ind(1:length(ind)-1,1);
	commModellist{i}.S(ind,find(ismember(commModellist{i}.rxns,'communityBiomass'))) = -nonzeros(abundancelist(i,:))
end

% union to get reactions col
U=commModellist{1}.rxns; %cell array containing arrays
for i=2:length(commModellist)
	temp = commModellist{i}
	U=union(U,temp.rxns);
end
U2 = cell2table([[commModellist{1}.rxns],[commModellist{1}.subSystems]] ,"VariableNames",["Reactions" "Subsystem"])

for i=2:length(commModellist)
	tmp = cell2table([[commModellist{i}.rxns],[commModellist{i}.subSystems]],"VariableNames",["Reactions" "Subsystem"])
	U2 = union(U2,tmp);
end

%%%%%% get fluxes
for i=1:length(commModellist)
	model = commModellist{i}

	% remove obj to do metabolicEP
	model.c(model.c==1)=0;
	model.c(find(ismember(model.rxns,'EX_microbeBiomass[fe]')))=1;
	sol = optimizeCbModel(model)
	ind = find(ismember(model.rxns,'EX_microbeBiomass[fe]'))
	model.lb(ind) = 0.1*sol.f
	model.c(model.c==1)=0;

	modelEP = pre_processing(model);
	% Set sampling parameters
	exp_i = 0;
	av_exp = 0;
	va_exp = 0;
	Beta=1e8;
	damping=0.9;
	precision=1e-5;
	minvar=1e-50;
	maxvar=1e50;
	maxit=1e3;

	% metabolicEP
	[mu_free_un, s_free_un, a_free_un, d_free_un, av_free_un, va_free_un, cov_free_un, t_EP_free_un] = ...
	MetabolicEP(full(modelEP.S), modelEP.b, modelEP.lb, modelEP.ub, Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);
	%modelind = find(ismember(U,modelEP.rxns)); %ismember cannot cause will
	%sort
	modelind=[];
	for ii =1:length(modelEP.rxns)
		modelind(ii,1) = find(strcmp(U,modelEP.rxns{ii}));
	end
	commetmu = cell(length(U),1);
	commetmu(modelind)=num2cell(mu_free_un);
	newmatrixmet_mu(i) = {commetmu'};
end

resultTable = cell2table(vertcat(newmatrixmet_mu{:})');
finalTable = [U2,resultTable]
writetable(finalTable,'muresults.csv')

save('results.mat','commModellist','commfluxlist','microbeNames','abundancelist')
