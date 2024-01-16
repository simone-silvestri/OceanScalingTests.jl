using JLD2

Na = (4320, 1800, 100)
Ny = (4320, 1801, 100)
Nz = (4320, 1800, 101)
Nη = (4320, 1800, 1)

T = Tuple(zeros(Na) for i in 1:12)
S = Tuple(zeros(Na) for i in 1:12)
u = Tuple(zeros(Na) for i in 1:12)
v = Tuple(zeros(Ny) for i in 1:12)
w = Tuple(zeros(Nz) for i in 1:12)
η = Tuple(zeros(Nη) for i in 1:12)

myfiles = readdir("./")
myfiles = filter(x -> length(x) >= 21, myfiles)
myfiles = filter(x -> x[1:21] == "compressed_iteration_", myfiles)

@inline remove_last_character(s) = s[1:end-1]
myfiles = remove_last_character.(myfiles)

numbers = parse.(Int, filter.(isdigit, myfiles))
numbers = sort(numbers)

monthly_days(year) = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

total_days = [sum([monthly_days(1995)[i][end] for i in 1:j]) for j in 1:12]

months = zeros(12)

for num in numbers[1:3]
    file = jldopen("compressed_iteration_" * string(num) * ".jld2")
    day  = mod(file["clock"].time / 86400, 365)
    mnth = findlast(x -> x >= day, total_days) 

    T[mnth] .+= file["T"]
    S[mnth] .+= file["S"]
    u[mnth] .+= file["u"]
    v[mnth] .+= file["v"]
    w[mnth] .+= file["w"]
    η[mnth] .+= file["η"]

    months[mnth] += 1
end

jldsave("climatology_data.jld2", T = T, S = S, u = u, v = v, w = w, η = η, months = months)

