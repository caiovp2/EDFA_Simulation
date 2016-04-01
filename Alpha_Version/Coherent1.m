%   Executa a simulação de camada física para um arquivo XML retornado
%       pela Ferramenta de Planejamento. 
%
%                                    by B. Coutinho
%                                      05/01/2015

%% Inicialização da Simulação e Leitura do XML
PLSim_Input_Data;





%% Variáveis auxiliares

% Variável auxiliar para o cálculo do tamanho da DCF
fatorDL_DCF = 0.0;

% Variável auxiliar para o cálculo da atenuação acumulada
fatorAL_ALFA = 0.0;

%% Buscando as informações dos Enlaces a serem analisados

% Por conta de ter de começar de um determinado enlace, teremos de rodar
% todos os enlaces pelo menos 2 vezes




        % Verificando a interface que inicia o enlace
        id_interface = SegmentoFibra.id_interface_in(SegmentoFibra.id == id_segmentofibra(1));

        % Retornar o nº de canais que propagarão neste enlace
        canais_por_enlace = tab_canal_enlace(:,Enlace.id(ii));
        canais_por_enlace = find(canais_por_enlace > 0);
        Nch = length(canais_por_enlace);
        fprintf(fid,'Nº de Canais do inicio do Enlace = %d \n',Nch);

        % Inicialização necessária do Optilux
        reset_all(Nsymb,Nt,Nch);

        % Montando a Sequência de símbolos para a transmissão
        patmatx = cell(Nsymb,2);
        patmaty = cell(Nsymb,2);
        for jj=1:Nch
            A = pattern('debruijn',randi(100,1),struct('alphabet',2));
            X_I = A';
            X_Q = horzcat(X_I(33:Nsymb),X_I(1:32));
            Y_I = horzcat(X_I(200:Nsymb),X_I(1:199));
            Y_Q = horzcat(Y_I(33:Nsymb),Y_I(1:32));
            % opção: as 2 linhas comentadas abaixo são os comandos para 
            % gerar sequências sem atrasos - necessário comentar as linhas 
            % restantes deste loop.
            %[~,patmatx{jj}] = pattern('random',1,struct('alphabet',4));
            %[~,patmaty{jj}] = pattern('random',2,struct('alphabet',4));
            patmatx{jj}(:,1) = X_I';
            patmatx{jj}(:,2) = X_Q';
            patmaty{jj}(:,1) = Y_I';
            patmaty{jj}(:,2) = Y_Q';
        end

        % Gerando o sinal do laser
        % Necessidade de buscar quais lambdas estão sendo usados neste enlace
        lambdas = zeros(1,length(canais_por_enlace));
        for jj=1:length(canais_por_enlace)
            fprintf(fid,'Id de Canal Utilizado neste Enlace: %d -> ',canais_por_enlace(jj));
            canal_utilizado = CanalOptico.lambda_utilizado{CanalOptico.id == canais_por_enlace(jj)};
            fprintf(fid,'Label do Canal = %s\n',canal_utilizado);
            lambdas(jj) = Banda.valor(strcmp(Banda.label,canal_utilizado));
        end
        
        % A divisão por 40 está relacionada com os equipamentos da PadTec
        %Pin = dBm2Watt(Enlace.pin(ii))*1e3/40; % [mW]
        Pin = 1;
        E = lasersource(Pin,lambdas,Inf, options);
        %E = lasersource(Pin,1550.12,1.6, options);
        fprintf(fid,'\n Potencia de Entrada = %f mW \n',Pin);
        
        % Buscando o nó de origem da fibra que inicia o
        % enlace para a checagem dos canais que vieram de enlaces
        % anteriores
        id_noh_origem = SegmentoFibra.id_noh_origem(SegmentoFibra.id == id_segmentofibra(1));
        fprintf(fid,'\n ID do noh origem deste Enlace = %d \n',id_noh_origem);
        % buscando as fibras que estão chegando neste noh origem
        fibras_in_noh_origem = SegmentoFibra.id(SegmentoFibra.id_noh_destino == id_noh_origem);
        enlaces_in_noh_origem = SegmentoFibra.id_enlace(fibras_in_noh_origem);
        fprintf(fid,'\n ID dos enlaces que chegam neste noh origem = %d \n',enlaces_in_noh_origem);
        canais_in_noh_origem = tab_canal_enlace(:,enlaces_in_noh_origem);
        [row,col] = find(canais_in_noh_origem < 3 & canais_in_noh_origem ~= 0);
        canais_in = row;
        if ~isempty(canais_in)
            fprintf(fid,'\n ID dos canais que chegam neste noh origem = %d \n',canais_in);
        end
        
        elecx_i = zeros(Nt*Nsymb,Nch);
        elecx_q = zeros(Nt*Nsymb,Nch);
        elecy_i = zeros(Nt*Nsymb,Nch);
        elecy_q = zeros(Nt*Nsymb,Nch);
        
        Eoptx = zeros(Nt*Nsymb,Nch);
        Eopty = zeros(Nt*Nsymb,Nch);
        
        % Modulação do sinal
        % buscar quais canais estão vindo de enlaces anteriores
        for jj = 1:Nch
            
            elecx_i(:,jj) = electricsource(patmatx{jj}(:,1),'qpsk',symbrate,'cosroll',duty,roll);
            elecx_q(:,jj) = electricsource(patmatx{jj}(:,2),'qpsk',symbrate,'cosroll',duty,roll);
            elecy_i(:,jj) = electricsource(patmaty{jj}(:,1),'qpsk',symbrate,'cosroll',duty,roll);
            elecy_q(:,jj) = electricsource(patmaty{jj}(:,2),'qpsk',symbrate,'cosroll',duty,roll);
                
            % checar se o canal está vindo do enlace anterior, checar o nº
            % da rodada também
            if  nn == 1 || tab_canal_enlace(canais_por_enlace(jj),Enlace.id(ii)) == 1 || tab_canal_enlace(canais_por_enlace(jj),Enlace.id(ii)) == 4
                fprintf(fid,'\n O Canal %d serah gerado \n',canais_por_enlace(jj));
                Eoptx(:,jj)   = qi_modulator(E(:,jj), elecx_i(:,jj), elecx_q(:,jj),struct('amplitude',[0.1 0.1]));
                Eopty(:,jj)   = qi_modulator(E(:,jj), elecy_i(:,jj), elecy_q(:,jj),struct('amplitude',[0.1 0.1]));
            else
                % buscar o canal do enlace anterior
                id_canal_anterior = find(row == canais_por_enlace(jj));
                id_enlace_anterior = enlaces_in_noh_origem(col(id_canal_anterior));
                fprintf(fid,'\n O Canal %d serah recuperado do enlace %d \n',canais_por_enlace(jj),id_enlace_anterior);
                load(fullfile('tmp',strcat('X_canal',int2str(canais_por_enlace(jj)),'_enlace',int2str(id_enlace_anterior))),'X');
                load(fullfile('tmp',strcat('Y_canal',int2str(canais_por_enlace(jj)),'_enlace',int2str(id_enlace_anterior))),'Y');
                Eoptx(:,jj) = X;
                Eopty(:,jj) = Y;
            end
        end

        % Criação do campo óptico
        create_field('unique',Eoptx,Eopty,struct('power','average'));

        % Transmissão em cada fibra do enlace
        for jj=1:length(id_segmentofibra)
            fprintf(fid,'\n\n ====  Fiber %d (Tx)\n',id_segmentofibra(jj));
            comp = SegmentoFibra.comprimento(SegmentoFibra.id == id_segmentofibra(jj));
            tx.length = comp*1e3; % length [m]
            fprintf(fid,'length  = %f [m] \n',comp*1e3);
            alfa = SegmentoFibra.alfa(SegmentoFibra.id == id_segmentofibra(jj));
            tx.alphadB = alfa/comp; % attenuation [dB/km]
            fatorAL_ALFA = fatorAL_ALFA + tx.alphadB*tx.length; % Calculo da atenuação acumulada
            fprintf(fid,'alphadB = %f  [dB/km] \n',tx.alphadB);
            tx.aeff    = 80; % [um^2]
            fprintf(fid,'aeff    = 80 effective area [um^2] \n');
            tx.n2      = 2.7e-20; % nonlinear index
            fprintf(fid,'n2      = 2.7e-20 \n');
            tx.lambda = SegmentoFibra.lambda_ref(SegmentoFibra.id == id_segmentofibra(jj)); % wavelength [nm] @ dispersion
            fprintf(fid,'lambda  = %f [nm] @ dispersion \n',SegmentoFibra.lambda_ref(SegmentoFibra.id == id_segmentofibra(jj)));
            tx.disp = SegmentoFibra.coef_dc(SegmentoFibra.id == id_segmentofibra(jj)); % dispersion [ps/nm/km] @ wavelength
            fprintf(fid,'disp    = %f [ps/nm/km] @ wavelength \n',SegmentoFibra.coef_dc(SegmentoFibra.id == id_segmentofibra(jj)));
            tx.slope = SegmentoFibra.s_0(SegmentoFibra.id == id_segmentofibra(jj)); % slope [ps/nm^2/km] @ wavelength
            fprintf(fid,'slope   = %f [ps/nm^2/km] @ wavelength \n',SegmentoFibra.s_0(SegmentoFibra.id == id_segmentofibra(jj)));
            tx.dphimax = 3e-3;    % maximum nonlinear phase rotation per step
            fprintf(fid,'dphimax = 3e-3 maximum nonlinear phase rotation per step \n');
            tx.dzmax   = 500;     % maximum SSFM step
            fprintf(fid,'dzmax   = 500 maximum SSFM step \n');
            coef_pmd = SegmentoFibra.coef_pmd(SegmentoFibra.id == id_segmentofibra(jj));
            tx.pmd = coef_pmd; % ps/sqrt(km)
            fprintf(fid,'tx.pmd     = %f ps/sqrt(km) \n',tx.pmd);
            tx.manakov = 'yes';   % enables the Manakov equation for random birefringence
            fprintf(fid,'manakov = yes enables the Manakov equation for random birefringence \n');
            tx.dgd = tx.pmd * sqrt(tx.length*1e-3) ...
                        * 1e-3 * symbrate; % DGD médio [symbols]
            tx.nplates = 80;
            fprintf(fid,'dgd = %f DGD medio [ps] \n',(tx.dgd/symbrate)*1e3);

            fatorDL_DCF = fatorDL_DCF + tx.disp*comp;
            
            tic
            if CUDA
                cu_fiber(tx,'gpsx'); % <fiber.m> do Optilux 0.1 adaptado para CUDA
            else
          %      fiber(tx,'gpsx'); % <fiber.m> do Optilux 0.1
                fiber(tx,'gp--'); % <fiber.m> do Optilux 0.1
            end
            toc
            
            fprintf(fid,'\n\n');
            
            % Buscando o ID do nó destino desta fibra
            id_noh_destino = SegmentoFibra.id_noh_destino(SegmentoFibra.id == id_segmentofibra(jj));

            % Verificando se o nó destino é um amplificador
            if ~isempty(Amplificador.id_noh(Amplificador.id_noh == id_noh_destino))
                fprintf(fid,'\n Noh destino %d eh de Amplificação \n',id_noh_destino);
                Gerbio = fatorAL_ALFA/1e3;
                fatorAL_ALFA = 0.0;
                fprintf(fid,'\n Ganho necessario para o amplificador compensar a atenuacao = %f dB \n',Gerbio);
                ampliflat(Gerbio,'gain');
            end

            % Verificando se o nó destino é DCM para usar as DCFs
            if ~isempty(DCM.id_noh(DCM.id_noh == id_noh_destino))
                fprintf(fid,'\n Noh destino %d eh DCM - Usar DCF \n',id_noh_destino);
                %%%%  Fiber 2 (compensating fiber)
                dcf.alphadB = 0.11;     % attenuation [dB/km]
                dcf.aeff    = 20;      % effective area [um^2]
                dcf.n2      = 2.7e-20; % nonlinear index
                dcf.lambda  = 1550;    % wavelength [nm] @ dispersion 
                dcf.disp    = -300;    % dispersion [ps/nm/km] @ wavelength
                dcf.slope   = -0.42;   % slope [ps/nm^2/km] @ wavelength
                dcf.dphimax = 3E-3;    % maximum nonlinear phase rotation per step
                dcf.dzmax   = 500;     % maximum SSFM step
                dcf.length = (-fatorDL_DCF)/dcf.disp*1e3; % comp. fiber length [m]
                fatorDL_DCF = 0.0;
                fprintf(fid,'Tamanho necessario da DCF = %f [m] \n',dcf.length);
                dcf.pmd = 0.1;
                dcf.dgd = dcf.pmd * sqrt(dcf.length*1e-3) ...
                        * 1e-3 * symbrate; % DGD médio [symbols]
                dcf.nplates = 10;
                dcf.manakov = 'yes';
                fprintf(fid,'dgd = %f DGD medio da DCF [ps] \n',(dcf.dgd/symbrate)*1e3);
                
                tic
                if CUDA
                    cu_fiber(dcf,'gpsx');
                else
                    fiber(dcf,'gp--');
                end
                toc
            end

            % Verificando se o nó destino é OADM para a recepção coerente
            if ~isempty(TerminalOADM.id_noh(TerminalOADM.id_noh == id_noh_destino))
                fprintf(fid,'\n Noh destino %d eh OADM \n',id_noh_destino);

                % Após o término do enlace inserir o ruído e proceder a recepção do
                % sinal
                fprintf(fid,'\n === INSERINDO O RUIDO AO SINAL NO RECEPTOR \n');
                SNRm = OSNR + 10*log10(2*12.5/symbrate);
                add_awgn(SNRm);

                % Verificar para cada canal quem será recepcionado neste OADM,
                % ou seja, finaliza neste nós OADM em questão, e quais
                % canais continuarão a propagação no próximo enlace

                % Faz um backup do GSTATE original por conta da função
                % optfilter
                BACKUP = copia_gstate();

                for kk=1:Nch

                    % Filtro do canal
                    fprintf(fid,'\n === Filtragem do Canal %d \n',canais_por_enlace(kk));
                    optfilter(kk,'supergauss',0.8*Sim.espacamento/symbrate,7);

                    % checar se o canal continuará no próximo enlace
                    if tab_canal_enlace(canais_por_enlace(kk),Enlace.id(ii)) > 2

                        % Parâmetros Originais
                        cma.N = CMA_MAX;
                        N = WINDOW_MAX;
                        
                        % Caso 1: O canal será finalizado / droppado neste OADM
                        fprintf(fid,'>>>>>>>>>>>>>>>Este canal Finaliza aqui \n!');

                        % Detecção por diversidade de fase
                        fprintf(fid,'\n === Detecção por diversidade de fase \n');
                        Irx = Det_Coer_PhaseDiversity(rec,options);

                        % Conversão A/D
                        fprintf(fid,'\n === ADC \n');
                        tic
                        [X,Y,~] = ADC(Irx,Nt,Nsymb,symbrate,samplingrate,enob,0);
                        toc
                        %[BER.estimadaQf.X, BER.estimadaQf.Y] = Ber_QF(X,Y);
                        %avgber_pol_adc = 0.5*(BER.estimadaQf.X + BER.estimadaQf.Y);
                        %fprintf(fid,'BER apos o ADC = %g\n',avgber_pol_adc);
                        [~,~,Qmin_adc] = Ber_QF(X,Y);
                        fprintf(fid,'Q minimo apos o ADC = %g [dB]\n',Qmin_adc);

                        % Gerando a Ortonormalização
                        if rec.ImpAngle ~= 0
                            fprintf(fid,'\n === Ortonormalização usando Algoritmo de Gram-Schmidt \n');
                            [X,Y] = Ortonormaliza(X,Y,0);
                        end
                        
                        % Compensação de Dispersão Cromática
                        fprintf(fid,'\n === Passando pelo Comp DC \n');
                        X_orig = X;
                        Y_orig = Y;
                        Nspan = length(id_segmentofibra);
                        tic
                        [X,Y] = Compensa_Disp_CM(X,Y,tx,299792458,Nspan,symbrate,samplingrate,0);
                        toc
                        %[BER.estimadaQf.X, BER.estimadaQf.Y] = Ber_QF(X,Y);
                        %avgber_pol_dc = 0.5*(BER.estimadaQf.X + BER.estimadaQf.Y);
                        [~,~,Qmin_dc] = Ber_QF(X,Y);
                        if (Qmin_dc <= Qmin_adc)
                            fprintf(fid,'O Compensador de DC nao foi necessario neste caso! \n');
                            X = X_orig;
                            Y = Y_orig;
                        else
                            fprintf(fid,'Q minimo apos a Compensação da DC = %g [dB]\n',Qmin_dc);
                        end

                        % Reamostra sinal para que haja 2 amostras por símbolo
                        fprintf(fid,'\n === Reamostrando o sinal para que haja 2 amostras por simbolo \n');
                        dp = 2;
                        X = ppinterp(1:length(X), X, 1:dp:length(X));
                        Y = ppinterp(1:length(Y), Y, 1:dp:length(Y));

                        % PolDemux
                        fprintf(fid,'\n === Passando pelo CMA \n');
                        % verificar se o Q-factor está dentro de um limite
                        % aceitável para este módulo
                        X_orig = X;
                        Y_orig = Y;
                        melhor_X_cma = [];
                        melhor_Y_cma = [];
                        xx = 1;
                        while cma.N > 2
                            [X,Y] = Pol_Demux(X_orig,Y_orig,cma,0);
                            %[BER.estimadaQf.X, BER.estimadaQf.Y] = Ber_QF(X,Y);
                            %avgber_pol_cma = 0.5*(BER.estimadaQf.X + BER.estimadaQf.Y);
                            [~,~,Qmin_cma] = Ber_QF(X,Y);
                            fprintf(fid,'Q minimo apos o CMA = %g [dB]\n',Qmin_cma);
                            % 1º Rodada
                            if xx == 1
                                melhor_qmin_cma = Qmin_cma;
                                first_X_cma = X;
                                first_Y_cma = Y;
                            end
                            if Qmin_cma >= melhor_qmin_cma
                                fprintf(fid,'Maior Q!! \n');
                                melhor_qmin_cma = Qmin_cma;
                                melhor_X_cma = X;
                                melhor_Y_cma = Y;
                            end
                            if Qmin_cma >= VALOR_Q_FEC
                                break;
                            else
                                cma.N = cma.N - 1;
                                fprintf(fid,'Valor do # Taps do CMA alterado para %d \n',cma.N);
                            end
                            xx = xx + 1;
                        end
                        
                        % Verificando o melhor resultado do CMA
                        if ~isempty(melhor_X_cma)
                            X = melhor_X_cma;
                            Y = melhor_Y_cma;
                        else
                            X = first_X_cma;
                            Y = first_Y_cma;
                        end

                        % Phase Estimation
                        fprintf(fid,'\n === Passando pelo Estimador de Fase \n');
                        % verificando se o Q final está dentro de um
                        % limite aceitável para este módulo
                        X_orig = X;
                        Y_orig = Y;
                        xx = 1;
                        while N >= 1
                            [X,Y] = Phase_Estimation(X_orig,Y_orig,N,0);
                            %[BER.estimadaQf.X, BER.estimadaQf.Y] = Ber_QF(X,Y);
                            %avgber_pol = 0.5*(BER.estimadaQf.X + BER.estimadaQf.Y);
                            [~,~,Qmin] = Ber_QF(X,Y);
                            fprintf(fid,'Q minimo apos o Estimador = %g [dB]\n',Qmin);
                            if xx == 1 
                                melhor_qmin_phase = Qmin;
                            end
                            if Qmin >= melhor_qmin_phase
                                fprintf(fid,'Maior Q!! \n');
                                melhor_qmin_phase = Qmin;
                            end
                            if Qmin >= VALOR_Q_FEC
                                break;
                            else
                                N = N - 1;
                                fprintf(fid,'Tamanho da Janela de simbolos alterado para %d \n',N);
                            end
                            xx = xx + 1;
                        end

                        fprintf(fid,'\n Q MINIMO (RODADA %d) DO CANAL %d NO FINAL DO ENLACE %d = %g [dB]\n',nn,canais_por_enlace(kk),Enlace.id(ii),melhor_qmin_phase);
                        q_canal_enlace(canais_por_enlace(kk),Enlace.id(ii)) = melhor_qmin_phase;
                        
                    else
                        
                        % Caso 2: O canal continuará no próximo enlace
                        fprintf(fid,'>>>>>>>>>>>>>>>Este canal continuarah no próximo enlace \n!');

                        % Armazenar o vetor do sinal deste canal para uso no
                        % próximo enlace
                        [X,Y] = retorna_canal();
                        save(fullfile('tmp',strcat('X_canal',int2str(canais_por_enlace(kk)),'_enlace',int2str(Enlace.id(ii)))),'X');
                        save(fullfile('tmp',strcat('Y_canal',int2str(canais_por_enlace(kk)),'_enlace',int2str(Enlace.id(ii)))),'Y');
                        
                    end % Interface out

                    % Restaura o GSTATE com todos os canais
                    restaura_gstate(BACKUP);

                end % Nch

            end % Terminal OADM

        end % Segmentos de Fibra

    