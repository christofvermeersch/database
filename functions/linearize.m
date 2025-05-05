function linmat = linearize(mat,supp)
    A = mat;
    n = size(supp,2);
    d = max(sum(supp,2)); 
B_size = 0;
    for di = 0:d-1
        B_size = B_size + nchoosek(n+di-1,di);
    end

    % A op de correcte manier gebruiken, om de B's te bouwen
    [k, l] = size(A{1});

    % Hoeveel B's er zijn, hangt af van het aantal parameters
    B = cell(1,n+1);
    B0_11 = A{1};
    B0_22 = kron(eye(B_size-1),-eye(l));
    B0_12 = zeros(k, length(B0_22));
    B0_21 = zeros(length(B0_22), l);
    B{1} = [B0_11, B0_12; B0_21, B0_22];
    
    aantal = ones(d,n); 
    %Voor elke graad van standaardmonomen (rij)...
    for di = 1:d
        %... ga voor elke parameter \lambda_i na, hoeveel keer deze
        %voorkomt, zet dan op deze plaatsen de juiste A op de juiste plek
        for ni = 1:n     
            aantal(di,ni) = nchoosek((n-ni)+(di-1),(di-1));
        end
    end
    dim = sum(aantal(:,1));

    % We hebben A0 niet meer nodig, maar wel een nulmatrix als 'opvulling'
    % in de B-matrices, dus die zetten we even op deze plaats
    A{1} = zeros(k,l);
    % De array eerste_rij bevat steeds de 
    first_row = cell(dim);
    A_index = 2;
    space_index = 0;
    for di = 1:d
        % Voor elke rij kijken we hoeveel 'blokken' er toegevoegd worden 
        % aan de B-matrices
        new_space = aantal(di,1);
        for ni = 1:n             
            % Daarna kijken we voor de B-matrix behorend bij parameters,
            % welke er daarvan opgevuld moet worden
            fill_in = aantal(di,ni);
            % De eerste new_space - fill_in blokken worden leeg opgevuld,
            % de laatste fill_in blokken worden met delen van A ingevuld
            grens = new_space-fill_in;
            for i = 1:grens
                first_row{ni, space_index+i} = A{1};
            end
            for i = grens+1:new_space
                first_row{ni, space_index+i} = A{A_index};
                A_index = A_index + 1;
            end
        end
        space_index = space_index + new_space;
    end

    % De array eerste_rij bevat steeds de 
    enen = eye(l,l);
    nullen = zeros(l,l);
    identity_rows = cell(n*(B_size-1),B_size);
    identity_rows(:) = {nullen};
    start_row = 1;
    start_column = 1;
    for di = 1:d-1
        % Het aantal kolommen dat we in deze stap zullen gaan toevoegen is: 
        new_columns = aantal(di,1);
        % We gaan voor elke parameter
        for ni = 1:n
            % over de nieuwe rijen, maar houden hierbij in het achterhoofd
            % dat elke nieuwe rij slechts 1 keer ingevuld kan worden. Om
            % deze reden wordt er bij de start_row steeds het aantal
            % ingevulde rijen opgeteld. Er wordt per 'degree block' ook
            % over de nieuwe kolommen geitereerd, en alles wordt van achter
            % naar voor ingevuld. Nadat dit voor 1 parameter gebeurd is,
            % schakelen we voor een nieuwe graad over naar een nieuwe
            % start_column.
            fill_in_rows = aantal(di,ni);
            for index = new_columns-1:-1:new_columns-fill_in_rows
                % Over nieuwe rijen voor deze parameter (het verschil is om
                % ervoor te zorgen dat je niet te ver schuift als je niet
                % voldoende toevoegt)
                nr = start_row + index - (new_columns-fill_in_rows);
                % en over nieuwe kolommen
                nc = start_column + index;
                % en die opvullen met enen
                identity_rows{(ni-1)*(dim-1)+nr,nc} = enen;
            end
            start_row = start_row + fill_in_rows;
        end
        start_column = start_column + new_columns;
    end
    
    % De eerste rij en het eenheidsgedeelte uit de cell arrays halen en in
    % een matrix zetten
    B_index = 2;
    for ni = 1:n
        Bi = zeros(k, dim*l);
        Bi_rij = zeros(l, dim*l);
        Bi_rest = zeros((dim-1)*l, dim*l);
        for dimi = 1:dim
            Bi(:,(dimi-1)*l+1:dimi*l) = first_row{ni,dimi};
        end
        for dimi = 1:dim-1
            for dimj = 1:dim
                Bi_rij(:,(dimj-1)*l+1:dimj*l) = identity_rows{(dim-1)*(ni-1)+dimi,dimj};
            end
            Bi_rest((dimi-1)*l+1:dimi*l,:) = Bi_rij;
            Bi_rij = zeros(l, dim*l);
        end
        B{B_index} = [Bi; Bi_rest];
        B_index = B_index + 1;
    end

    linmat = B; 
end


