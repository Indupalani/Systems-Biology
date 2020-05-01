function model = modelParser(model,genes,head,node_parent,parent,printLevel)
% This function look into the GPRparser and builds auxiliary nodes for each
% layer of the GPR rules.
% For example:
%     g1 & g2 & g3 is one layer in the model,   RXN = g1 + g2 + g3
%     (g1 | g2) & g3 is a two layer model:      RXN = (g1 | g2) + g3
%                                               (g1 | g2) = g1
%                                               (g1 | g2) = g2



%  Name nodes according to position and layer inside the GPR rule
node = strcat({[node_parent,'_']},cellfun(@num2str,num2cell(1:length(head.children)),'UniformOutput', false));

% model = addMetabolite(model, node);
model.mets = vertcat(model.mets,node');
[~, idx_new_mets] = ismember(node, model.mets);

idx_parent_mets = find(strcmp(parent, model.mets));

%  Expand model with every node, using OR/AND rules
if isa(head,'OrNode')
    idx_new_rxn = length(model.rxns) + (1:length(head.children));
    model.rxns(idx_new_rxn) = node;
    model.S(idx_new_mets,idx_new_rxn) = -eye(length(head.children));
    model.S(idx_parent_mets,idx_new_rxn) = +1;
else
    idx_new_rxn = length(model.rxns)+1;
    model.rxns{idx_new_rxn} = ['AND_' parent];
    model.S(idx_new_mets,idx_new_rxn) = -1;
    model.S(idx_parent_mets,idx_new_rxn) = +1;
end


for i=1:length(head.children)
    child = head.children(i);
    % check if the child is final layer
    if isa(child,'LiteralNode')
        gen = genes(str2double(child.id));
        idx_new_rxn = length(model.rxns)+1;
        idx_met_gen = strcmp(gen,model.mets);
        idx_met_node = strcmp(node(i),model.mets);
        model.rxns{idx_new_rxn} = [gen{1},'_',node{i}];
        model.S(idx_met_node,idx_new_rxn) = +1;
        model.S(idx_met_gen,idx_new_rxn) = -1;
    else
        model = modelParser(model,genes,child,node{i},node{i},printLevel);
    end
end

end


function model = xyz(model,genes,head,node_parent,parent,printLevel)
% This function look into the GPRparser and builds auxiliary nodes for each
% layer of the GPR rules.
% For example:
%     g1 & g2 & g3 is one layer in the model,   RXN = g1 + g2 + g3
%     (g1 | g2) & g3 is a two layer model:      RXN = (g1 | g2) + g3
%                                               (g1 | g2) = g1
%                                               (g1 | g2) = g2



%  Name nodes according to position and layer inside the GPR rule
node = strcat({[node_parent,'_']},cellfun(@num2str,num2cell(1:length(head.children)),'UniformOutput', false));

% model = addMetabolite(model, node);
model.mets = vertcat(model.mets,node');
[~, idx_new_mets] = ismember(node, model.mets);

idx_parent_mets = find(strcmp(parent, model.mets));

%  Expand model with every node, using OR/AND rules
if isa(head,'OrNode')
    idx_new_rxn = length(model.rxns) + (1:length(head.children));
    model.rxns(idx_new_rxn) = node;
    model.S(idx_new_mets,idx_new_rxn) = -eye(length(head.children));
    model.S(idx_parent_mets,idx_new_rxn) = +1;
else
    idx_new_rxn = length(model.rxns)+1;
    model.rxns{idx_new_rxn} = ['AND_' parent];
    model.S(idx_new_mets,idx_new_rxn) = -1;
    model.S(idx_parent_mets,idx_new_rxn) = +1;
end


for i=1:length(head.children)
    child = head.children(i);
    % check if the child is final layer
    if isa(child,'LiteralNode')
        gen = genes(str2double(child.id));
        idx_new_rxn = length(model.rxns)+1;
        idx_met_gen = strcmp(gen,model.mets);
        idx_met_node = strcmp(node(i),model.mets);
        model.rxns{idx_new_rxn} = [gen{1},'_',node{i}];
        model.S(idx_met_node,idx_new_rxn) = +1;
        model.S(idx_met_gen,idx_new_rxn) = -1;
    else
        model = modelParser(model,genes,child,node{i},node{i},printLevel);
    end
end

end