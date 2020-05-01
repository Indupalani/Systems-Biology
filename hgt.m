%
% hgt.m
%
% Horizontal Gene Transfer simulation
%

maxT = 1e4; % Number of timesteps

b_total = 1e6; % cells/mL
b_avg = b_total / n;

if runs_i == 0
    fsIterations = 1;
end
fullSpread = [];

I = zeros(fsIterations, n, maxT);
for fsIterator = 1:fsIterations

    b = ones(n, 1) * b_avg;

    % Conjuagtion:
    CONJUGATION = zeros(n, n);

    for recipient = 1:n

        if rand(1) < prob_super_permissive
            bonus = 1;
        else
            bonus = 1/super_permissive_bonus;
        end

        for donor = 1:n
            CONJUGATION(recipient, donor) = gamma_conjugation * rand(1) * bonus;
        end
    end

    % Transformation:
    TRANSFORMATION = zeros(n, n);

    while sum(sum(TRANSFORMATION)) == 0 && gamma_transformation > 0
        for index = 1:n
            if rand(1) < prob_transformable
                %fprintf('Species %d is competent.\n', index);
                TRANSFORMATION(index, :) = gamma_transformation;
            end
        end
    end

    % Transduction:
    TRANSDUCTION = zeros(n, n);

    module_start = 1;
    module_size = max_module_size;
    while module_start <= n

        % handle cluster of bacteria indexed as:
        % (cluster_start) thru (cluster_start + cluster_size - 1)

        for recipient = module_start:min(n, module_start + module_size - 1)
            for donor = module_start:min(n, module_start + module_size - 1)
                if recipient >= donor
                    overlap = 1;
                else
                    overlap = (recipient - module_start + 1)/(donor - module_start + 1);
                end
                TRANSDUCTION(recipient, donor) = overlap * gamma_transduction;
            end
        end

        module_start = module_start + module_size;
        module_size = module_size - 1;
    end

    % Vesicles:
    VESICLES = zeros(n, n);

    donor_efficiency = rand(1, n) .* (rand(1, n) < prob_produce_vesicles);
    recipient_efficiency = rand(1, n) .* (rand(1, n) < prob_uptake_vesicles);

    for recipient = 1:n
        for donor = 1:n
            VESICLES(recipient, donor) = gamma_vesicles * recipient_efficiency(recipient) * donor_efficiency(donor);
        end
    end

    % Constant: sum HGT Relationships Matrix
    CONJUGATION = acceleration * CONJUGATION;
    TRANSFORMATION = acceleration * TRANSFORMATION;
    TRANSDUCTION = acceleration * TRANSDUCTION;
    VESICLES = acceleration * VESICLES;
    A = CONJUGATION + TRANSFORMATION + TRANSDUCTION + VESICLES;

    % Keep track of total amount of each mechanism performed [cells/mL]
    total_conjugation = 0;
    total_transformation = 0;
    total_transduction = 0;
    total_vesicle_mediated = 0;
    total_vertical = 0;

    % Take note when we reach full spread
    done_mark = -1;

    % Count of plasmid+ bacteria in each species (starting state)
    i = zeros(n, maxT); % cells/mL
    i(1, 1) = b(1)/10;

    for t = 2:maxT

        % compute the actual change due to HGT
        di = (b - i(:, t-1)) .* (A * i(:, t-1));
        i(:, t) = i(:, t-1) + di;

        temp = exp(acceleration * m / generation_time);
        VERTICAL = temp*i(:, t) ./ (1 + (i(:, t)./b) * (temp - 1)) - i(:, t);
        i(:, t) = temp*i(:, t) ./ (1 + (i(:, t)./b) * (temp - 1));

        % Don't allow i to surpass b (also, see if we've reached full spread):
        done_count = 0;
        for index = 1:n
            if i(index, t) > b(index)*0.99 % arbitrary threshold
                i(index, t) = b(index);
                done_count = done_count + 1;
            end
        end

        if done_mark < 0 && done_count == n
            done_mark = t;
        end

        % record the contribution of each mechanism
        total_conjugation = total_conjugation + sum( (b - i(:, t-1)) .* (CONJUGATION * i(:, t-1)) );
        total_transformation = total_transformation + sum( (b - i(:, t-1)) .* (TRANSFORMATION * i(:, t-1)) );
        total_transduction = total_transduction + sum( (b - i(:, t-1)) .* (TRANSDUCTION * i(:, t-1)) );
        total_vesicle_mediated = total_vesicle_mediated + sum( (b - i(:, t-1)) .* (VESICLES * i(:, t-1)) );
        total_vertical = total_vertical + sum(VERTICAL);
    end
    
    if done_mark > 0
        fullSpread = [fullSpread done_mark];
    end
    
    I(fsIterator, :, :) = i;
end

if runs_i > 0
    i(1, maxT) = b(1);
    theplot = subplot(runs_m, runs_n, runs_i);
    theplot.OuterPosition(1) = theplot.OuterPosition(1) - 0.0500;
    theplot.Position(2) = theplot.Position(2) + 0.0000;
    theplot.Position(4) = theplot.Position(4) - 0.0200;
    imagesc(bsxfun(@rdivide, i, b)); % normalize i by the size of the species
    set(gca, 'XTick', []);
    set(gca, 'FontSize', 12);
    %colorbar
    ylabel('Species', 'FontSize', 13);
    xlabel('Time', 'FontSize', 13);
    
    % print percentage total contribution of each mechanism
    total_all = total_conjugation + total_transformation + total_transduction + total_vesicle_mediated;% + total_vertical;
    fprintf('\n');
    fprintf('Conjugation:      %8.6f\n', total_conjugation/total_all);
    fprintf('Transformation:   %8.6f\n', total_transformation/total_all);
    fprintf('Transduction:     %8.6f\n', total_transduction/total_all);
    fprintf('Vesicle-Mediated: %8.6f\n', total_vesicle_mediated/total_all);
    fprintf('Vertical:         %8.6f\n', total_vertical/total_all);

    % print total time simulated
    fprintf('\n');
    days = maxT * acceleration / (60*24);
    if days < 365.2425
        fprintf('Time Simulated: %d days\n', days);
    else
        fprintf('Time Simulated: %d years\n', days/365.2425);
    end
    
    % Mark time axis
    xticks = [1 maxT/2 maxT];
    xticklabels = {'0' format_time_interval(maxT/2 * acceleration) format_time_interval(maxT * acceleration)};
    
    set(gca, 'XTick', xticks);
    set(gca, 'XTickLabels', xticklabels);
    
    if ~isempty(fullSpread) % Mark full spread
        mean_fullSpread = mean(fullSpread);
        
        if fsIterations > 1 % draw patch
            std_fullSpread = std(fullSpread);
            patchCoordsX = [mean_fullSpread-std_fullSpread mean_fullSpread-std_fullSpread mean_fullSpread+std_fullSpread mean_fullSpread+std_fullSpread];
            patchCoordsY = [0 n+1 n+1 0];
            patch(patchCoordsX, patchCoordsY, 'g', 'LineStyle', 'none', 'FaceAlpha', 0.8);
        else % draw line
            line([fullSpread fullSpread], get(gca,'YLim'), 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
        end
        
        if std_fullSpread >= 3e2 % if there's room for a label in the stripe
            meanLabel = text(mean_fullSpread - 0.005*maxT, n/2, format_time_interval(mean_fullSpread * acceleration), 'FontSize', 11, 'Rotation', 90, 'HorizontalAlignment', 'center');
            if fsIterations <= 1
                meanLabel.Color = 'w';
                meanLabel.VerticalAlignment = 'top';
            end
        end
        fprintf('Run %i: N=%i, mean=%d minutes, SD=%d minutes\n', runs_i, fsIterations, mean_fullSpread*acceleration, std_fullSpread*acceleration);
    end
    fprintf('\n---\n\n');
end
