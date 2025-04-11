struct MyStructT{N}
    data::NTuple{N,Int}
end

struct MyStructNT{T}
    data::@NamedTuple{T}
end

struct MyStructA
    data::Vector{Int}
end
# Example usage
my_instance = MyStruct((1, 2, 3))
println(my_instance)

function UnivariateStatistic2{T,K,I}(weights::I, rawmoments...) where {T,K,I}
    nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
    new{T,K,I,typeof(rawmoments)}(weights, rawmoments)
end