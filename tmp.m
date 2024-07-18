N_para_dis = []


for Nd_i = [1,5]
    for Nd_j = [1,5]
        for Nc_i = [1,2]
            for Nc_j = [1,2]
                temp = [Nd_i;Nc_i;Nd_j;Nc_j;abs(Nd_i + 5.*Nc_i - Nd_j - 5.*Nc_j)];
                N_para_dis = [N_para_dis,temp];

            end

        end
    end
end