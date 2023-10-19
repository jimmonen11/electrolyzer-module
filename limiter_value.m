function val = limiter_value(r,method)
%Returns flux limiter values for different methods

% deal with the situation where 'r' is a scalar rather than an array
% still have some python in here haven't updated for different schemes
if method == "UPWIND"
    val = 0;

elseif method == "CENTRAL"
    val = 1;

% elif method == 'CHARM':
%     val = np.zeros_like(r)
%     ix = np.where(r>0)
%     if len(ix)>0:
%         val[ix] = r[ix] * (3*r[ix] + 1) / (r[ix] + 1)**2

elseif method == "HCUS"
    val = zeros(1, length(r));
    for i = 1:length(val)
        val(i) = 1.5*( r(i) + abs(r(i)) ) / (r(i)+2);
    end

elseif method == "KOREN"
    val = zeros(1, length(r));
    for i = 1:length(val)
        val(i) = max(0 , min(2*r(i), min( (2+r(i))/3, 2)) );
    end
% elif method == 'MINMOD':
%     val = np.zeros_like(r)
%     for i in range(len(val)):
%         val[i] = max(0, min(1, r[i]))
% 
% elif method == 'SUPERBEE':
%     val = np.zeros_like(r)
%     for i in range(len(val)):
%         val[i] = max(0, min(2*r[i], 1), min(r[i],2) )
% 
% elif method == 'VANLEER':
%     val = np.zeros_like(r)
%     for i in range(len(val)):
%         val[i] = (r[i]+abs(r[i]))/(1+abs(r[i]))

else
    assert(False) % unsupported option
end