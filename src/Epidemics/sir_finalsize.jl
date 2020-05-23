
using Plots

@fastmath function SIR_fs(N,bet,gamm)
    final_size = zeros(N+1,1);
    final_size[2] = 1;
    for Z2 in 0:N
        @inbounds for Z1 in Z2+1:N-1
            p1 = 1 / ( 1 + gamm/(bet*(N-Z1)));
            final_size[Z1+2] = final_size[Z1+2] + final_size[Z1+1]*p1;
            final_size[Z1+1] = final_size[Z1+1]*(1-p1);
        end
    end
    return final_size;
end

N = 20;
bet = 2/(N-1);
gamm = 1.0;

final_size = SIR_fs(N,bet,gamm);
@time final_size = SIR_fs(N,bet,gamm);

bar(0:N,final_size,label="Cumulative infections")
