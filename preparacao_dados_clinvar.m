function [dados_treino, dados_teste, resumo] = preparacao_dados_clinvar(ficheiro_entrada, ficheiro_mat, max_linhas)
% Pré-processa o ClinVar para o fluxo do sistema.

    if nargin < 1 || isempty(ficheiro_entrada)
        ficheiro_entrada = 'clinvar.csv';
    end
    if nargin < 2 || isempty(ficheiro_mat)
        ficheiro_mat = 'dados_clinvar_processados.mat';
    end
    if nargin < 3
        max_linhas = [];
    end

    fprintf('A carregar ClinVar de %s...\n', ficheiro_entrada);

    opcoes = detectImportOptions(ficheiro_entrada, 'TextType', 'string');
    opcoes.SelectedVariableNames = selecionar_colunas_uteis(opcoes.VariableNames);

    if ~isempty(max_linhas)
        opcoes.DataLines = [2, max_linhas + 1];
    end

    dados = readtable(ficheiro_entrada, opcoes);

    col_chrom = procurar_nome_coluna(dados.Properties.VariableNames, ["CHROM", "Chrom", "chrom"]);
    col_ref = procurar_nome_coluna(dados.Properties.VariableNames, ["REF", "Ref", "ref"]);
    col_alt = procurar_nome_coluna(dados.Properties.VariableNames, ["ALT", "Alt", "alt"]);
    col_info = procurar_nome_coluna(dados.Properties.VariableNames, ["INFO", "Info", "info"]);
    col_clnsig = procurar_nome_coluna(dados.Properties.VariableNames, ["CLNSIG", "Clnsig", "clinical_significance"]);

    for idx = 1:height(dados)
        linha = dados(idx, :);
        [gene, consequencia] = extrair_campos_clinvar(obter_valor_campo(linha, col_info));
        dados.Gene(idx, 1) = gene;
        dados.Consequence(idx, 1) = consequencia;
        dados.Perfil(idx, 1) = build_clinvar_variant_string(linha);
        dados.variant_key(idx, 1) = compor_chave_variante(linha, col_chrom, col_ref, col_alt);
        if strlength(col_clnsig) > 0
            valor_clnsig = string(linha.(col_clnsig));
            dados.CLNSIG(idx, 1) = valor_clnsig;
            dados.is_high_risk(idx, 1) = any(contains(lower(valor_clnsig), ["pathogenic", "likely_pathogenic", "risk_factor", "conflicting_interpretations_of_pathogenicity"]));
        else
            dados.CLNSIG(idx, 1) = "";
            dados.is_high_risk(idx, 1) = false;
        end
    end

    colunas_a_verificar = {};
    if strlength(col_chrom) > 0, colunas_a_verificar{end + 1} = char(col_chrom); end
    if strlength(col_ref) > 0, colunas_a_verificar{end + 1} = char(col_ref); end
    if strlength(col_alt) > 0, colunas_a_verificar{end + 1} = char(col_alt); end
    if strlength(col_info) > 0, colunas_a_verificar{end + 1} = char(col_info); end
    if ~isempty(colunas_a_verificar)
        dados = rmmissing(dados, 'DataVariables', colunas_a_verificar);
    end

    dados = removevars(dados, intersect(string(dados.Properties.VariableNames), ["INFO"]));

    rng(0);
    total = height(dados);
    indices = randperm(total);
    n_treino = max(1, round(0.8 * total));
    dados_treino = dados(indices(1:n_treino), :);
    dados_teste = dados(indices(n_treino + 1:end), :);

    resumo = struct();
    resumo.total = total;
    resumo.treino = height(dados_treino);
    resumo.teste = height(dados_teste);
    resumo.colunas = string(dados.Properties.VariableNames);
    resumo.qtd_high_risk = sum(dados.is_high_risk);

    dados_treino = reordenar_para_fluxo(dados_treino);
    dados_teste = reordenar_para_fluxo(dados_teste);

    save(ficheiro_mat, 'dados_treino', 'dados_teste', 'resumo');

    fprintf('ClinVar preparado com sucesso:\n');
    fprintf('  - Total de amostras processadas: %d\n', total);
    fprintf('  - Treino: %d\n', resumo.treino);
    fprintf('  - Teste: %d\n', resumo.teste);
    fprintf('  - Variantes de alto risco: %d\n', resumo.qtd_high_risk);
    fprintf('  - Guardado em: %s\n', ficheiro_mat);
end

function selecionadas = selecionar_colunas_uteis(nomes_variaveis)
    desejadas = ["CHROM", "REF", "ALT", "INFO", "CLNSIG", "Varyant_Tipi", "Variant_Tipi", "Target", "CHROM_Encoded", "VT_Encoded"];
    nomes_variaveis = string(nomes_variaveis);
    selecionadas = strings(0, 1);

    for i = 1:numel(desejadas)
        idx = find(lower(nomes_variaveis) == lower(desejadas(i)), 1);
        if ~isempty(idx)
            selecionadas(end + 1, 1) = nomes_variaveis(idx);
        end
    end
end

function nome_campo = procurar_nome_coluna(nomes_variaveis, sinonimos)
    nome_campo = "";
    nomes_variaveis = string(nomes_variaveis);
    sinonimos = string(sinonimos);

    for i = 1:numel(sinonimos)
        idx = find(lower(nomes_variaveis) == lower(sinonimos(i)), 1);
        if ~isempty(idx)
            nome_campo = nomes_variaveis(idx);
            return;
        end
    end
end

function valor = obter_valor_campo(linha, nome_campo)
    valor = "";
    if strlength(nome_campo) > 0
        valor = string(linha.(nome_campo));
    end
end

function chave = compor_chave_variante(linha, col_chrom, col_ref, col_alt)
    if strlength(col_chrom) == 0 || strlength(col_ref) == 0 || strlength(col_alt) == 0
        chave = "";
        return;
    end

    chave = string(linha.(col_chrom)) + ":" + string(linha.(col_ref)) + ">" + string(linha.(col_alt));
end

function tabela = reordenar_para_fluxo(tabela)
    preferidas = ["Gene", "CHROM", "REF", "ALT", "Consequence", "CLNSIG", "variant_key", "Perfil", "is_high_risk"];
    disponiveis = string(tabela.Properties.VariableNames);
    ordem = strings(0, 1);

    for i = 1:numel(preferidas)
        idx = find(lower(disponiveis) == lower(preferidas(i)), 1);
        if ~isempty(idx)
            ordem(end + 1, 1) = disponiveis(idx);
        end
    end

    if ~isempty(ordem)
        restantes = setdiff(disponiveis, ordem, 'stable');
        tabela = tabela(:, [cellstr(ordem); cellstr(restantes(:))]);
    end
end