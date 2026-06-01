% app_oncologia.m
% Interface Grafica do Pipeline de Triagem Clinico-Oncologica - MPEI

function app_oncologia()
    addpath('utils');
    
    dados_treino = [];
    dados_teste = [];
    assinaturas_treino = [];
    minhash = [];
    filtro_bloom = [];
    nb = [];
    
    fig = uifigure('Name', 'Sistema de Triagem Clinico-Oncologica - MPEI', ...
        'Position', [100 100 800 620]);
    
    % 1. Carregamento de Dados e Treino
    pnlControl = uipanel(fig, 'Title', '1. Carregamento de Dados e Treino', 'Position', [20 500 760 100]);
    lblStatus = uilabel(pnlControl, 'Text', 'Estado: Aguardando carregamento...', 'Position', [20 40 400 22], 'FontWeight', 'bold');
    btnLoad = uibutton(pnlControl, 'Text', 'Carregar ClinVar e Treinar Modelos', ...
        'Position', [450 30 280 35], 'ButtonPushedFcn', @(btn, event) carregarETreinar(), ...
        'BackgroundColor', [0.1 0.5 0.8], 'FontColor', [1 1 1], 'FontWeight', 'bold');
        
    % 2. Selecionar Paciente de Teste
    pnlPaciente = uipanel(fig, 'Title', '2. Selecionar Paciente de Teste', 'Position', [20 280 760 200]);
    uilabel(pnlPaciente, 'Text', 'Escolha um Paciente do Conjunto de Teste:', 'Position', [20 145 300 22]);
    ddPacientes = uidropdown(pnlPaciente, 'Position', [20 115 720 25], 'Items', {'(Carregue os dados primeiro)'}, 'Enable', 'off', ...
        'ValueChangedFcn', @(dd, event) atualizarDetalhesPaciente());
    lblDetalhes = uilabel(pnlPaciente, 'Text', 'Detalhes do paciente selecionado aparecerao aqui.', 'Position', [20 15 720 85], 'WordWrap', 'on', 'FontAngle', 'italic');
    
    % 3. Resultados do Pipeline
    pnlResultados = uipanel(fig, 'Title', '3. Resultados do Pipeline Integrado', 'Position', [20 20 760 240]);
    btnComparar = uibutton(pnlResultados, 'Text', 'Executar Pipeline Completo', ...
        'Position', [20 180 240 32], 'Enable', 'off', 'ButtonPushedFcn', @(btn, event) executarPipeline(), ...
        'BackgroundColor', [0.1 0.7 0.4], 'FontColor', [1 1 1], 'FontWeight', 'bold');
        
    % Sub-paineis esquerdos (Bloom e Naive Bayes)
    pnlBloom = uipanel(pnlResultados, 'Title', 'Filtro de Bloom (Triagem)', 'Position', [20 100 320 70]);
    lblBloomVal = uilabel(pnlBloom, 'Text', 'Aguardando execucao...', 'Position', [10 10 300 30], ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        
    pnlNB = uipanel(pnlResultados, 'Title', 'Naive Bayes (Diagnostico)', 'Position', [20 20 320 70]);
    lblNBVal = uilabel(pnlNB, 'Text', 'Aguardando execucao...', 'Position', [10 10 300 30], ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        
    % Tabela de Resultados (MinHash)
    tblResultados = uitable(pnlResultados, 'Position', [350 20 390 180]);
    tblResultados.ColumnName = {'Rank', 'Similaridade', 'Gene', 'Classe Historica'};
    tblResultados.ColumnWidth = {50, 90, 80, 140};
    
    % Callbacks
    function carregarETreinar()
        btnLoad.Enable = 'off';
        lblStatus.Text = 'A ler clinvar.csv e a preparar dados...';
        drawnow;
        
        ficheiro_mat = 'dados_clinvar_processados.mat';
        try
            if ~isfile(ficheiro_mat)
                [dados_treino, dados_teste, ~] = preparacao_dados_clinvar('clinvar.csv', ficheiro_mat, 5000);
            else
                load(ficheiro_mat, 'dados_treino', 'dados_teste');
            end
            
            dados_treino(dados_treino.CLNSIG == "" | ismissing(dados_treino.CLNSIG), :) = [];
            dados_teste(dados_teste.CLNSIG == "" | ismissing(dados_teste.CLNSIG), :) = [];
            
            lblStatus.Text = 'A treinar modelos...';
            drawnow;
            
            % 1. Inicializar e Treinar Filtro de Bloom
            n_alto_risco = sum(dados_treino.is_high_risk);
            m_bloom = max(1000, round(8 * n_alto_risco));
            k_bloom = 5;
            filtro_bloom = FiltroBloom(m_bloom, k_bloom);
            chaves_alto_risco = dados_treino.variant_key(dados_treino.is_high_risk);
            for i = 1:numel(chaves_alto_risco)
                filtro_bloom.Inserir(chaves_alto_risco(i));
            end
            
            % 2. Inicializar e Treinar MinHash
            num_hashes = 150;
            tamanho_shingle = 5;
            minhash = MinHash(num_hashes, tamanho_shingle, 42);
            n_treino = height(dados_treino);
            assinaturas_treino = zeros(n_treino, num_hashes, 'uint64');
            for i = 1:n_treino
                assinaturas_treino(i, :) = minhash.gerarAssinatura(char(dados_treino.Perfil(i)));
            end
            
            % 3. Inicializar e Treinar Naive Bayes
            nb = NaiveBayes(1);
            nb.treinar(dados_treino, 'CLNSIG', {'Gene', 'CHROM', 'REF', 'ALT', 'Consequence'});
            
            lblStatus.Text = sprintf('Dados Carregados! Treino Concluido (%d amostras)', n_treino);
            
            n_teste = height(dados_teste);
            itens_dropdown = cell(n_teste, 1);
            for i = 1:n_teste
                itens_dropdown{i} = sprintf('Paciente %d: Gene %s (Mut: %s>%s)', ...
                    i, dados_teste.Gene(i), dados_teste.REF(i), dados_teste.ALT(i));
            end
            ddPacientes.Items = itens_dropdown;
            ddPacientes.Enable = 'on';
            btnComparar.Enable = 'on';
            atualizarDetalhesPaciente();
            
        catch excepcao
            lblStatus.Text = ['Erro: ' excepcao.message];
            btnLoad.Enable = 'on';
        end
    end

    function atualizarDetalhesPaciente()
        if isempty(dados_teste), return; end
        
        idx_selecionado = ddPacientes.Value;
        idx_partes = sscanf(idx_selecionado, 'Paciente %d:');
        idx_paciente = idx_partes(1);
        
        paciente = dados_teste(idx_paciente, :);
        lblDetalhes.Text = sprintf('<b>Gene:</b> %s | <b>Cromossoma:</b> %s | <b>Mutacao:</b> %s > %s | <b>Consequencia:</b> %s<br><b>Significancia Clinica Real:</b> %s<br><b>Perfil:</b> %s', ...
            paciente.Gene, paciente.CHROM, paciente.REF, paciente.ALT, paciente.Consequence, paciente.CLNSIG, paciente.Perfil);
        lblDetalhes.Interpreter = 'html'; 
     end

    function executarPipeline()
        if isempty(dados_teste) || isempty(minhash) || isempty(filtro_bloom) || isempty(nb), return; end
        
        idx_selecionado = ddPacientes.Value;
        idx_partes = sscanf(idx_selecionado, 'Paciente %d:');
        idx_paciente = idx_partes(1);
        
        paciente = dados_teste(idx_paciente, :);
        
        % 1. Filtro de Bloom (Triagem)
        e_risco = filtro_bloom.Verificar(paciente.variant_key);
        if e_risco
            lblBloomVal.Text = 'ALERTA: Variante de Risco!';
            lblBloomVal.BackgroundColor = [1 0.3 0.3];
            lblBloomVal.FontColor = [1 1 1];
        else
            lblBloomVal.Text = 'Sem perigo conhecido no Bloom';
            lblBloomVal.BackgroundColor = [0.4 0.8 0.4];
            lblBloomVal.FontColor = [0 0 0];
        end
        
        % 2. Naive Bayes (Diagnostico)
        [classe, confianca, ~] = nb.classificar(paciente);
        lblNBVal.Text = sprintf('%s (%.1f%%)', classe, confianca * 100);
        
        if contains(lower(classe), 'pathogenic')
            lblNBVal.BackgroundColor = [1 0.6 0.2];
            lblNBVal.FontColor = [1 1 1];
        elseif contains(lower(classe), 'benign')
            lblNBVal.BackgroundColor = [0.7 0.9 0.7];
            lblNBVal.FontColor = [0 0 0];
        else
            lblNBVal.BackgroundColor = [0.9 0.9 0.9];
            lblNBVal.FontColor = [0 0 0];
        end
        
        % 3. MinHash (Semelhancas)
        perfil_novo = char(paciente.Perfil);
        assinatura_novo = minhash.gerarAssinatura(perfil_novo);
        
        n_treino = height(dados_treino);
        similaridades = zeros(n_treino, 1);
        for j = 1:n_treino
            similaridades(j) = minhash.estimarJaccard(assinatura_novo, assinaturas_treino(j, :));
        end
        
        [valores_sim, indices_top] = maxk(similaridades, 5);
        
        dados_tabela = cell(5, 4);
        for k = 1:5
            idx_historico = indices_top(k);
            paciente_historico = dados_treino(idx_historico, :);
            
            dados_tabela{k, 1} = sprintf('Top %d', k);
            dados_tabela{k, 2} = sprintf('%.1f%%', valores_sim(k) * 100);
            dados_tabela{k, 3} = char(paciente_historico.Gene);
            dados_tabela{k, 4} = char(paciente_historico.CLNSIG);
        end
        tblResultados.Data = dados_tabela;
    end
end
