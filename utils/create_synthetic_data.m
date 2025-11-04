function [outputArg1,outputArg2] = create_synthetic_data(inputArg1,inputArg2)

create_synthetic_data_sigmoid();
% create_synthetic_data_sustain();

end

function create_synthetic_data_sustain()
% The actual code is in python and is located in
% ./competing_methods/pySuStaIn-master/sim/create_synthetic_data_sustain.py

cmd1_str = 'cd ./competing_methods/pySuStaIn-master/sim';
cmd2_str = [get_python_path(), ' create_synthetic_data_sustain.py'];
system([cmd1_str, ' && ', cmd2_str]);
end

function create_synthetic_data_sigmoid()
clear;
save_dir = './input/Synthetic_data/synthetic_data_Sigmoid/';
alp_mu = 10;
alp_sgm = 2;
bta_low = 0.1;
bta_upp = 0.9;
nsamp = 3000; % total number of points
nsamp_control = 600; % number of NC points (located at 0)
% nsamp_control = 0;
nsbtyp = 3;

nbiom_list = [5,10,15,20,25,30];
noise_list = [0.01,0.1,0.2];
thres = 0.5 * ones(nsbtyp,1);
%nbiom_list = 10;
%noise_list = 0.2;

for noise = noise_list
    for nbiom = nbiom_list
        for f = 1:10
            folder = [save_dir,'Sigmoid_b',num2str(nbiom),'_n', ...
                num2str(noise),'_f',num2str(f),'/'];
            mkdir(folder);

            data = zeros(nsamp,nbiom);
            stage = rand(nsamp,1);
            stage(1:nsamp_control) = 0;
            subtype = randi(nsbtyp,nsamp,1);
            traj_data = cell(nsbtyp,nbiom,2);
            ord = zeros(nsbtyp,nbiom);

            %% create the trajectory
            for j = 1:nsbtyp
                comp = zeros(nbiom,1);
                %figure;
                hold on
                for i = 1:nbiom
                    traj_data(j,i,:) = create_traj_sigmoid(alp_mu,alp_sgm,bta_low,bta_upp);
                    comp(i) = sum(traj_data{j,i,2}< thres(j),2);
                    %plot_traj(traj_data{j,i});
                end
                %saveas(figure(j),join(['.\',folder,'\figure_',int2str(j),'.png']));
                hold off
                [~,ord(j,:)] = sort(comp);
                %disp(ord(j,:))
            end

            %% generate the data point with random noise
            for k = 1:nsamp
                j = subtype(k);
                for i = 1:nbiom
                    % the noise is scaled by the maximum of the trajectory
                    data(k,i) = normrnd(interp_traj(stage(k),traj_data(j,i,:)),noise*5);
                end
            end

            writematrix(data,[folder,'/data.csv']);
            writematrix(stage,[folder,'/stages.csv']);
            writematrix(subtype, [folder, '/subtypes.csv']);
            writematrix(ord,[folder,'/order.csv']);
            for j = 1:nsbtyp
                writematrix(cell2mat(traj_data(j,:,2)')', [folder,'/traj_',num2str(j),'.csv']);
            end
        end
    end
end

end

function traj_data = create_traj_sigmoid(alp_mu,alp_sgm,bta_low,bta_upp)
traj_x = linspace(0, 1, 1000);

alpha = normrnd(alp_mu,alp_sgm);
beta = bta_low + (bta_upp-bta_low)*rand;

traj_y = 1 ./ (1 + exp(-alpha * (traj_x - beta)));

% shift and scale s.t. the first point starts at 0 and the maximum is 5
traj_y = traj_y - traj_y(1);
traj_y = traj_y/max(traj_y) * 5;

traj_data = {traj_x, traj_y};

end


function h = plot_traj(traj_data)
stage_t = 0:0.001:1;
h = plot(stage_t,interp_traj(stage_t,traj_data), 'linewidth', 1);
end

function f = interp_traj(s,traj_data)    
f = interp1(traj_data{1},traj_data{2},s);
end


