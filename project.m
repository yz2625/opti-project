load Cow_PartialDistance
num_pt = size(Dist,1);

omega = ones(num_pt, num_pt);
for ind1 = 1:num_pt
    for ind2 = 1:num_pt
        omega(ind1, ind2) = ind1 == ind2 || Dist(ind1, ind2) ~= 0;
    end
end
omega = sparse(omega);

% Initialization
nonzero = nnz(omega);
tau = 3000; % parameter
delta = 1.2 * (num_pt * num_pt) / nonzero; % step size
Y = zeros(num_pt, num_pt); % Y0 = 0
Y = sparse(Y);
iteration = 200; % maximum iteration
l = 5; % increment
r = 0;

err = zeros(1, iteration);
err(1) = inf;


for itr = 1:iteration
    s = r + 1;
    while true
        U = sparse(num_pt, s);
        S = sparse(s, s);
        V = sparse(s, num_pt);
        [U, S, V] = svds(Y, s);
        
        s = s + l;
        if S(s-l, s-l) <= tau
            break
        end
    end
        r = max(find(diag(S) < tau));
        S = shrink(S, tau);
        X = U * S * V';
        Y = Y + delta * ((Dist - X).*omega);
        err(itr) = norm((X - Dist).*omega,'fro')/norm(Dist.*omega,'fro');
end

X = (X + X.')/2; % symmetry; round off error
[a, ~] = size(X);
for i = 1:a
    X(i,i)=0; % zeros along the diagonal
end

plot(err);
plt = mdscale(X, 3);
figure(2);
ViewMesh(plt, trg);