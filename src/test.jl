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