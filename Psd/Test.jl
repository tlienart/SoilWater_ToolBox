using Winston

x = 0:0.1:10
y = sin.(x)
y2 = sin.(2sin.(2sin.(x)))
a=plot(x, y, "g^", x, y2, "b-o")
plot(x, y, "g^", x, y2, "b-o")
Winston.savefig(a, "Test.svg")