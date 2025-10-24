
function runDirichletExamples()
    % Run the twoelement_dirichletboundary.m script
    disp('Running twoelement_dirichletboundary.m...');
    twoelement_dirichletboundary; % Ensure this script is in the current path

    % Run the dirichlet_pdetoolbox.m script
    disp('Running dirichlet_pdetoolbox.m...');
    dirichlet_pdetoolbox; % Ensure this script is in the current path

    disp('Both scripts have been executed.');
end