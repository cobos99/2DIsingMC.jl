using ArgParse
using Printf
include("2disingmc.jl")
using .mc2dising

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "L"
            help = "Linear dimension of the 2D lattice to be simulated"
            arg_type = Int
            required = true
        "temp"
            help = "Temperature in Joules (T*k_B)"
            arg_type = Float64
            required = true
        "J"
            help = "Coupling constant"
            arg_type = Float64
            required = true
        "exec_time"
            help = "Simulation execution time per site"
            arg_type = Int
            required = true
        "--init", "-i"
            help = "Initialization strategy. Possible values: r (random), f (ferromagnetic), fn (negative ferromagnetic), <filename>"
            arg_type = String
            default = "r"
        "--seach", "-s"
            help = "Time per site interval between saves"
            arg_type = Int
            default = 1
    end
    
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    filename = @sprintf("MC2DIsing_%ix%i_T_%.2f_J_%.2f_texecps_%i.txt", parsed_args["L"], parsed_args["L"], parsed_args["temp"], parsed_args["J"], parsed_args["exec_time"])
    if parsed_args["init"] == "r"
        initial_config = mc2dising.random_initialization(parsed_args["L"])
    elseif parsed_args["init"] == "f"
        initial_config = mc2dising.ferromagnetic_initialization(parsed_args["L"])
    elseif parsed_args["init"] == "fn"
        initial_config = mc2dising.ferromagnetic_initialization(parsed_args["L"], true)
    else
        initial_config = mc2dising.read_last_config(parsed_args["init"])
    end
    println("Running simulation...")
    mc2dising.run_mc2dising(filename, initial_config, parsed_args["J"], parsed_args["temp"],
                  parsed_args["exec_time"], parsed_args["seach"], print_mode=true)
end

main()