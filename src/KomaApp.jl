module CLI

import ..KomaMRI: CLI_HELP, run_cli

function print_help()
    print(CLI_HELP)
    return nothing
end

function (@main)(ARGS)
    if length(ARGS) == 1 && ARGS[1] in ("-h", "--help")
        return print_help()
    end
    return run_cli(ARGS)
end

end
