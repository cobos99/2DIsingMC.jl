module mc2dising

    using DelimitedFiles
    using Printf

    function save_results_formatted(outfile::String, results_matrix)
        to_save_string = (x -> (@sprintf("%.8e", x))).(results_matrix)
        open(outfile, "w") do out
            writedlm(out, to_save_string, '\t', quotes=false)
        end
    end
    
    function run_mc2dising(outfile::String, initial_config::Matrix{Int}, 
                           J::Real, temperature::Real, 
                           execution_time_per_site::Int, 
                           saveeach_time_per_site::Int;
                           save_to_file_each::Int=-1,
                           print_mode::Bool=false)

        if (length(size(initial_config)) != 2) ||
           (size(initial_config)[1] != size(initial_config)[2]) ||
           (any(abs.(initial_config) .!== 1))

           throw(ArgumentError("initial_config is not a 2D square matrix with entries ±1")) 
        
        end

        L = size(initial_config)[1]
        nsites = L^2
        energy = 0.

        if save_to_file_each < 1
            save_to_file_each::Int = execution_time_per_site ÷ 10
        end

        # Compute the initial energy
        for col in 1:L
            @simd for row in 1:L
                top_neigh = mod1(row - 1, L)
                right_neigh = mod1(col + 1, L)
                @inbounds energy += -J*initial_config[row, col]*(initial_config[top_neigh, col] + initial_config[row, right_neigh])
            end
        end

        # Compute the initial magnetization
        magnetization = sum(initial_config)

        # Initialize an array to store simulation results and the last spin_config
        results = zeros(div(execution_time_per_site, saveeach_time_per_site) + 1, 3)
        results[1, :] = [0., energy, magnetization]
        last_spin_config = copy(initial_config)
        last_spin_config_filename = prod(split(outfile, ".")[1:end-1] .* ".")[1:end-1] * "_lspinconf.txt"

        # Precompute the value of the exponentials
        # for the acceptance probabilities
        exps_dict = Dict{Tuple{Int, Int}, AbstractFloat}()
        for sμ_f in -1:2:1
            for sum_sμ_nn in -4:2:4
                    ΔE = 2*sμ_f*J*sum_sμ_nn
                    exps_dict[(sμ_f, sum_sμ_nn)] = exp(-ΔE/temperature)
            end
        end

        # Run the Monte Carlo simulation
        spin_config = copy(initial_config)
        for t_per_site in 1:execution_time_per_site

            if print_mode
                pad = length(string(execution_time_per_site))
                current_time = string(t_per_site, pad=pad)
                print(@sprintf("t_per_site = %s/%i\r", current_time, execution_time_per_site))
            end

            for t_this_interval in 1:nsites
                # Pick one spin position at random
                flip_row, flip_col = rand(1:L, 2)
                og_spin_val = spin_config[flip_row, flip_col]

                # Compute the energy after flipping the spin
                top_neigh = mod1(flip_row - 1, L)
                bot_neigh = mod1(flip_row + 1, L)
                right_neigh = mod1(flip_col + 1, L)
                left_neigh = mod1(flip_col - 1, L)
                neigh_sum = spin_config[top_neigh, flip_col] + spin_config[bot_neigh, flip_col] + spin_config[flip_row, right_neigh] + spin_config[flip_row, left_neigh]
                
                # Compute the acceptance probability for the chosen transition
                if neigh_sum*og_spin_val > 0
                    acceptance_prob = exps_dict[(og_spin_val, neigh_sum)]
                else
                    acceptance_prob = 1.0
                end

                # Roll de dice to check whether to flip and
                # store the new values of the magnetization and energy
                # if flipped
                rn = rand()
                if rn < acceptance_prob
                    spin_config[flip_row, flip_col] *= -1
                    ΔE = 2*og_spin_val*J*neigh_sum
                    energy = energy + ΔE
                    magnetization = magnetization - 2*og_spin_val
                end                          
            end
            # Save the energy and magnetization to a file
            if execution_time_per_site % saveeach_time_per_site == 0
                results[div(t_per_site, saveeach_time_per_site) + 1, :] = [t_per_site, energy, magnetization]
                last_spin_config = copy(spin_config)
                if execution_time_per_site % save_to_file_each == 0
                    save_results_formatted(outfile, results)
                    writedlm(last_spin_config_filename, last_spin_config, '\t')
                end
            end
        end
        save_results_formatted(outfile, results)
        writedlm(last_spin_config_filename, last_spin_config, '\t')
        if print_mode
            print("\n")
        end
    end

    function ferromagnetic_initialization(L::Int; negspins::Bool=false)::Matrix{Int}
        return (-1)^negspins .* ones(Int, L, L)
    end

    function random_initialization(L::Int)::Matrix{Int}
        return rand([-1, 1], L, L)
    end

    function read_last_config(filename::String)::Matrix{Int}
        return readdlm(filename, '\t', Int)
    end
end
