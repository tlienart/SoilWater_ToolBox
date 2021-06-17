using CairoMakie

    x = 1:10
    lines(x, x.^2, color = :black, linewidth = 2, linestyle = ".-", label = "x²",
    figure = (resolution = (700,450), fontsize = 18, backgroundcolor = "#D0DFE6FF"),
    axis = (xlabel = "x", ylabel = "x²", backgroundcolor = :white))
    axislegend("legend", position = :lt)
    limits!(0,10,0,100)
    current_figure()
    #save("lines1.png"), current_figure(), px_per_unit = 2.0)

