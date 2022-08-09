
function visualize(bs::RegularBsplines)
    # evalute bsplines
    t = range(bs.lower, bs.upper, length = 500)
    Bs = basis(t, bs)
    # plot
    Mke.lines(t, Bs[:, 1])
    for i = 2:size(Bs, 2)
        Mke.lines!(t, Bs[:, i])
    end
    Mke.current_figure()
end

visualize(bs)
Mke.save("plop.pdf", Mke.current_figure())

