using OnlineSampleStatistics
using StatsBase

const README = joinpath(@__DIR__, "..", "README.md")
const EXCLUDED_LINE_PATTERNS = [
    r"@btime",
    r"100_000_000",
    r"BenchmarkTools",
]

function extract_usage_blocks(path::AbstractString)
    lines = readlines(path)
    in_usage = false
    in_block = false
    lang = ""
    blocks = NamedTuple{(:lang, :code, :index), Tuple{String, String, Int}}[]
    current = IOBuffer()
    block_index = 0

    for line in lines
        if startswith(line, "## ")
            in_usage = (strip(line) == "## Usage")
            if !in_usage && !isempty(blocks)
                break
            end
        end

        if !in_usage
            continue
        end

        if startswith(line, "```")
            fence = strip(line)
            if !in_block
                lang = replace(fence, "```" => "")
                in_block = true
                current = IOBuffer()
                block_index += 1
            else
                code = String(take!(current))
                if strip(lang) == "julia" || strip(lang) == "julia-repl"
                    push!(blocks, (lang = strip(lang), code = code, index = block_index))
                end
                in_block = false
                lang = ""
            end
            continue
        end

        in_block && println(current, line)
    end

    return blocks
end

function repl_to_code(block::AbstractString)
    code = IOBuffer()
    for line in split(block, '\n')
        startswith(line, "julia> ") && println(code, line[8:end])
    end
    return String(take!(code))
end

function filter_heavy_lines(block::AbstractString)
    filtered = IOBuffer()
    for line in split(block, '\n')
        any(p -> occursin(p, line), EXCLUDED_LINE_PATTERNS) && continue
        println(filtered, line)
    end
    return String(take!(filtered))
end

function main()
    blocks = extract_usage_blocks(README)
    isempty(blocks) && error("No Julia/JULIA-REPL example blocks found in README Usage section")

    executed = 0
    for block in blocks
        raw_code = block.lang == "julia-repl" ? repl_to_code(block.code) : block.code
        code = filter_heavy_lines(raw_code)
        isempty(strip(code)) && continue

        try
            Base.include_string(Main, code, "README:usage:block$(block.index)")
            executed += 1
        catch err
            @error "README example block failed" block_index = block.index language = block.lang exception = (err, catch_backtrace())
            rethrow()
        end
    end

    return println("README Usage examples: OK (", executed, " block(s) executed)")
end

main()
