using StructTypes

abstract type Vehicle end

struct Car <: Vehicle
    type::String
    make::String
    model::String
    seatingCapacity::Int
    topSpeed::Float64
end

struct Truck <: Vehicle
    type::String
    make::String
    model::String
    payloadCapacity::Float64
end

# example from StructTypes deserialization
Car = StructTypes.read("""
{
    "type" "car",
    "make" "Mercedes-Benz",
    "model" "S500",
    "seatingCapacity" 5,
    "topSpeed" 250.1
}""", Vehicle)