function [dados_treino, dados_teste, resumo] = preparacao_dados_clinvar(inputFile, outputMat, maxRows)
%   Pré-processa o ClinVar para o fluxo do sistema.
%   Extrai Gene e Consequence do INFO (extrair_campos_clinvar.m),
%   constrói o perfil da variante (build_clinvar_variant_string.m),
%   separa treino/teste e devolve uma tabela pronta para os módulos.

    if nargin < 1 || isempty(inputFile)
        inputFile = 'clinvar.csv';
    end
    if nargin < 2 || isempty(outputMat)
        outputMat = 'dados_clinvar_processados.mat';
    end
    if nargin < 3
        maxRows = [];
    end

    fprintf('A carregar ClinVar de %s...\n', inputFile);

    opts = detectImportOptions(inputFile, 'TextType', 'string');
    opts.SelectedVariableNames = selectUsefulColumns(opts.VariableNames);

    if ~isempty(maxRows)
        opts.DataLines = [2, maxRows + 1];
    end

    dados = readtable(inputFile, opts);

    chromCol = findColumnName(dados.Properties.VariableNames, ["CHROM", "Chrom", "chrom"]);
    refCol = findColumnName(dados.Properties.VariableNames, ["REF", "Ref", "ref"]);
    altCol = findColumnName(dados.Properties.VariableNames, ["ALT", "Alt", "alt"]);
    infoCol = findColumnName(dados.Properties.VariableNames, ["INFO", "Info", "info"]);
    clnsigCol = findColumnName(dados.Properties.VariableNames, ["CLNSIG", "Clnsig", "clinical_significance"]);

    for idx = 1:height(dados)
        linha = dados(idx, :);
        [gene, consequence] = extrair_campos_clinvar(getFieldValue(linha, infoCol));
        dados.Gene(idx, 1) = gene; %#ok<AGROW>
        dados.Consequence(idx, 1) = consequence; %#ok<AGROW>
        dados.Perfil(idx, 1) = build_clinvar_variant_string(linha); %#ok<AGROW>
        dados.variant_key(idx, 1) = composeVariantKey(linha, chromCol, refCol, altCol); %#ok<AGROW>
        if strlength(clnsigCol) > 0
            clnsigValue = string(linha.(clnsigCol));
            dados.CLNSIG(idx, 1) = clnsigValue; %#ok<AGROW>
            dados.is_high_risk(idx, 1) = any(contains(lower(clnsigValue), ["pathogenic", "likely_pathogenic", "risk_factor", "conflicting_interpretations_of_pathogenicity"])); %#ok<AGROW>
        else
            dados.CLNSIG(idx, 1) = ""; %#ok<AGROW>
            dados.is_high_risk(idx, 1) = false; %#ok<AGROW>
        end
    end

    varsToCheck = {};
    if strlength(chromCol) > 0, varsToCheck{end + 1} = char(chromCol); end %#ok<AGROW>
    if strlength(refCol) > 0, varsToCheck{end + 1} = char(refCol); end %#ok<AGROW>
    if strlength(altCol) > 0, varsToCheck{end + 1} = char(altCol); end %#ok<AGROW>
    if strlength(infoCol) > 0, varsToCheck{end + 1} = char(infoCol); end %#ok<AGROW>
    if ~isempty(varsToCheck)
        dados = rmmissing(dados, 'DataVariables', varsToCheck);
    end

    dados = removevars(dados, intersect(string(dados.Properties.VariableNames), ["INFO"]));

    rng(0);
    total = height(dados);
    idx = randperm(total);
    nTreino = max(1, round(0.8 * total));
    dados_treino = dados(idx(1:nTreino), :);
    dados_teste = dados(idx(nTreino + 1:end), :);

    resumo = struct();
    resumo.total = total;
    resumo.treino = height(dados_treino);
    resumo.teste = height(dados_teste);
    resumo.colunas = string(dados.Properties.VariableNames);
    resumo.qtd_high_risk = sum(dados.is_high_risk);

    dados_treino = reorderForFlow(dados_treino);
    dados_teste = reorderForFlow(dados_teste);

    save(outputMat, 'dados_treino', 'dados_teste', 'resumo');

    fprintf('ClinVar preparado com sucesso:\n');
    fprintf('  - Total de amostras processadas: %d\n', total);
    fprintf('  - Treino: %d\n', resumo.treino);
    fprintf('  - Teste: %d\n', resumo.teste);
    fprintf('  - Variantes de alto risco: %d\n', resumo.qtd_high_risk);
    fprintf('  - Guardado em: %s\n', outputMat);
end

function selected = selectUsefulColumns(variableNames)
    desired = ["CHROM", "REF", "ALT", "INFO", "CLNSIG", "Varyant_Tipi", "Variant_Tipi", "Target", "CHROM_Encoded", "VT_Encoded"];
    variableNames = string(variableNames);
    selected = strings(0, 1);

    for i = 1:numel(desired)
        idx = find(lower(variableNames) == lower(desired(i)), 1);
        if ~isempty(idx)
            selected(end + 1, 1) = variableNames(idx); %#ok<AGROW>
        end
    end
end

function fieldName = findColumnName(variableNames, aliases)
    fieldName = "";
    variableNames = string(variableNames);
    aliases = string(aliases);

    for i = 1:numel(aliases)
        idx = find(lower(variableNames) == lower(aliases(i)), 1);
        if ~isempty(idx)
            fieldName = variableNames(idx);
            return;
        end
    end
end

function value = getFieldValue(row, fieldName)
    value = "";
    if strlength(fieldName) > 0
        value = string(row.(fieldName));
    end
end

function key = composeVariantKey(row, chromCol, refCol, altCol)
    if strlength(chromCol) == 0 || strlength(refCol) == 0 || strlength(altCol) == 0
        key = "";
        return;
    end

    key = string(row.(chromCol)) + ":" + string(row.(refCol)) + ">" + string(row.(altCol));
end

function tabela = reorderForFlow(tabela)
    preferred = ["Gene", "CHROM", "REF", "ALT", "Consequence", "CLNSIG", "variant_key", "Perfil", "is_high_risk"];
    available = string(tabela.Properties.VariableNames);
    order = strings(0, 1);

    for i = 1:numel(preferred)
        idx = find(lower(available) == lower(preferred(i)), 1);
        if ~isempty(idx)
            order(end + 1, 1) = available(idx); %#ok<AGROW>
        end
    end

    if ~isempty(order)
        remaining = setdiff(available, order, 'stable');
        tabela = tabela(:, [cellstr(order); cellstr(remaining(:))]);
    end
end