#=
using Pkg
Pkg.add("PrettyTables")
Pkg.add("SpecialFunctions")
Pkg.add("Formatting")
Pkg.add("Statistics")
Pkg.add("DataFrames")
Pkg.add("Plots")
Pkg.add("Printf")
Pkg.add("SymPy")
Pkg.add("PyCall")
=#


using PrettyTables
using SpecialFunctions
using Formatting
using Statistics
using DataFrames
using Plots
using Printf
using SymPy
using PyCall



x = Sym("x")
h = Sym("h")

#Controls the number of decimal places printed. Don't think I need ths anymore (was for graphing... mostley sign charts)
#Base.show(io::IO, f::Float64) = @printf(io, "%.8f", f)

#c style printf command
#@sprintf("%3.2f",convert(Float64,x))

"""
    limittable(f, a; rows::Int=5, dir::String="", format="%10.8f")

a: a finite number\n
rows: number of rows to compute (default is 5 rows)\n
dir: a string indicating which side to take the limit from\n
format: a string that specifies c-style printf format for numbers (default is %10.8f)

|dir|meaning|
|---|-------|
|""|approach from both sides (default)|
|"-"|approach from the left|
|"+"|approach from the right|
"""
function limittable(f, a; rows::Int=5, dir::String="", format="%10.8f")

    if dir == "+"
        X = a .+ [10.0^(-i) for i in 1:rows-2]
        X = vcat([a + 1, a + 0.5], X)
        X_str = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_str = [ sprintf1(format,y) for y in Y ]
        s=""
        for (x,y) in zip(X_str, Y_str)
            s *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre =
        """<table style="text-align:center; background-color: lightcoral">
        <caption style="text-align:center">Limit table: approach from the right</caption>
            <tr>
                <th scope: col>x</th>
                <th scope: col>y</th>
            </tr>"""
        table_str_post = """</table>"""
        table_str = table_str_pre*s*table_str_post
    elseif dir == "-"
        X = a .- [10.0^(-i) for i in 1:rows-2]
        X = vcat([a - 1, a - 0.5], X)
        X_str = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_str = [ sprintf1(format,y) for y in Y ]
        s=""
        for (x,y) in zip(X_str, Y_str)
            s *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre =
        """<table style="text-align:center; background-color: lightblue">
        <caption style="text-align:center">Limit table: approach from the left</caption>
            <tr>
                <th scope: col>x</th>
                <th scope: col>y</th>
            </tr>"""
        table_str_post = """</table>"""
        table_str = table_str_pre*s*table_str_post
    else
        Xr = a .+ [10.0^(-i) for i in 1:rows-2]
        Xr = vcat([a + 1, a + 0.5], Xr)
        Xr_str = [ sprintf1(format,x) for x in Xr ]
        Yr = [N(f(z)) for z in Xr]
        Yr_str = [ sprintf1(format,y) for y in Yr ]
        Xl = a .- [10.0^(-i) for i in 1:rows-2]
        Xl = vcat([a - 1, a - 0.5], Xl)
        Xl_str = [ sprintf1(format,x) for x in Xl ]
        Yl = [N(f(z)) for z in Xl]
        Yl_str = [ sprintf1(format,y) for y in Yl ]
        sl=""
        for (x,y) in zip(Xl_str, Yl_str)
            sl *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        sr=""
        for (x,y) in zip(Xr_str, Yr_str)
            sr *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre = """
        <table>
        <caption style="text-align:center">Limit table</caption>
        <tr>
            <th style="text-align:center; background-color: lightblue; scope: colgroup">Left</th>
            <th style="text-align:center; background-color: lightcoral; scope: colgroup">Right</th>
        </tr>
        <tr>
            <td>
                <table style="text-align:center; background-color: lightblue">
                  <tr>
                    <th scope: col>x</th>
                    <th scope: col>y</th>
                  </tr>
                  """
              table_str_mid = """
                </table>
            </td>
            <td>
                <table style="text-align:center; background-color: lightcoral">
                  <tr>
                    <th scope: col>x</th>
                    <th scope: col>y</th>
                  </tr>
              """
              table_str_post = """
                </table>
            </td>
        </tr>
        </table>
        """
        table_str = table_str_pre*sl*table_str_mid*sr*table_str_post
    end
    #display("text/markdown", table_str)
    display("text/html", table_str)
end

"""
    limittable(f, a::Sym; rows::Int=5, format="%10.2f")

a: is either oo (meaning &#8734;) or -oo (meaning -&#8734;)\n
rows: number of rows to compute (default is 5 rows)\n
dir: a string indicating which side to take the limit from\n
format: a string that specifies c-style printf format for numbers (default is %10.2f)
"""
function limittable(f, a::Sym; rows::Int=5, format="%10.2f")
    if a == oo
        X = [10.0^(i+1) for i in 1:rows]
        X_str = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_str = [ sprintf1(format,y) for y in Y ]
        s=""
        for (x,y) in zip(X_str, Y_str)
            s *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre =
        """<table style="text-align:center; background-color: lightcoral">
        <caption style="text-align:center">Limit table: x approaching &#8734;</caption>
            <tr>
                <th scope: col>x</th>
                <th scope: col>y</th>
            </tr>"""
        table_str_post = """</table>"""
        table_str = table_str_pre*s*table_str_post
    elseif a == -oo
        X = [-1*10.0^(i+1) for i in 1:rows]
        X_str = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_str = [ sprintf1(format,y) for y in Y ]
        s=""
        for (x,y) in zip(X_str, Y_str)
            s *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre =
        """<table style="text-align:center; background-color: lightblue">
        <caption style="text-align:center">Limit table: x approaching -&#8734;</caption>
            <tr>
                <th scope: col>x</th>
                <th scope: col>y</th>
            </tr>"""
        table_str_post = """</table>"""
        table_str = table_str_pre*s*table_str_post
    end

    #display("text/markdown", table_str)
    display("text/html", table_str)
end


################################################################################
#Trying to use PrettyTables for output
function limittable_pt(f, a; rows::Int=5, dir::String="", format="%10.8f", backend=:html)

    if dir == "+"
        X = a .+ [10.0^(-i) for i in 1:rows-2]
        X = vcat([a + 1, a + 0.5], X)
        X_data = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_data = [ sprintf1(format,y) for y in Y ]
        data = hcat(X_data, Y_data)
        header = ["x", "y"]
        even_row = HTMLHighlighter((data,i,j) -> (i%2==0), HTMLDecoration(color = "red"))
        pretty_table(data, header, backend = backend, highlighters = (even_row))
    elseif dir == "-"
        X = a .- [10.0^(-i) for i in 1:rows-2]
        X = vcat([a - 1, a - 0.5], X)
        X_str = [ sprintf1(format,x) for x in X ]
        Y = [N(f(z)) for z in X]
        Y_str = [ sprintf1(format,y) for y in Y ]
        s=""
        for (x,y) in zip(X_str, Y_str)
            s *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre =
        """<table style="text-align:center; background-color: lightblue">
        <caption style="text-align:center">Limit table: approach from the left</caption>
            <tr>
                <th scope: col>x</th>
                <th scope: col>y</th>
            </tr>"""
        table_str_post = """</table>"""
        table_str = table_str_pre*s*table_str_post
    else
        Xr = a .+ [10.0^(-i) for i in 1:rows-2]
        Xr = vcat([a + 1, a + 0.5], Xr)
        Xr_str = [ sprintf1(format,x) for x in Xr ]
        Yr = [N(f(z)) for z in Xr]
        Yr_str = [ sprintf1(format,y) for y in Yr ]
        Xl = a .- [10.0^(-i) for i in 1:rows-2]
        Xl = vcat([a - 1, a - 0.5], Xl)
        Xl_str = [ sprintf1(format,x) for x in Xl ]
        Yl = [N(f(z)) for z in Xl]
        Yl_str = [ sprintf1(format,y) for y in Yl ]
        sl=""
        for (x,y) in zip(Xl_str, Yl_str)
            sl *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        sr=""
        for (x,y) in zip(Xr_str, Yr_str)
            sr *= "<tr>
              <td>$x</td>
              <td>$y</td>
            </tr>"
        end
        table_str_pre = """
        <table>
        <caption style="text-align:center">Limit table</caption>
        <tr>
            <th style="text-align:center; background-color: lightblue; scope: colgroup">Left</th>
            <th style="text-align:center; background-color: lightcoral; scope: colgroup">Right</th>
        </tr>
        <tr>
            <td>
                <table style="text-align:center; background-color: lightblue">
                  <tr>
                    <th scope: col>x</th>
                    <th scope: col>y</th>
                  </tr>
                  """
              table_str_mid = """
                </table>
            </td>
            <td>
                <table style="text-align:center; background-color: lightcoral">
                  <tr>
                    <th scope: col>x</th>
                    <th scope: col>y</th>
                  </tr>
              """
              table_str_post = """
                </table>
            </td>
        </tr>
        </table>
        """
        table_str = table_str_pre*sl*table_str_mid*sr*table_str_post
    end

end
################################################################################


function lim(f, var, c; dir = "")
    lhl = limit(f, var, c, dir="-")
    rhl = limit(f, var, c, dir="+")
    if dir == ""
        rhl == lhl ? rhl : missing
    elseif dir == "+"
        rhl
    else
        lhl
    end
end

#Graphing stuff
function criticalpoints(f)
    crpt_num = solve(f)
    #crpt_num = N.(crpt_num)
    crpt_den = solve(simplify(1/f)) #without simplify failure
    #crpt_den = N.(crpt_den)
    sort(crpt_num), sort(crpt_den)
end

function getsign(f, val)
    f(val) > 0 ? "+" : "-"
end

function getsigns(f, crpt, values)
    signs = String[]

    m = Real[]

    push!(m, mean([values[1], crpt[1]]))
    push!(signs, getsign(f, m[1]))

    for i in 1:length(crpt)-1
        push!(m,mean([crpt[i], crpt[i+1]]))
        push!(signs, getsign(f, m[i+1]))
    end

    push!(m, mean([values[2], crpt[end]]))
    push!(signs, getsign(f, last(m)))
    m, signs
end

"""
    signchart(f, values, label="", horiz_jog = 0.2; dotverticaljog = 0.12, marksize = 6, imageFormat = :svg)
"""
function signchart(f, values, label="", horiz_jog = 0.2; size=(500, 200), dotverticaljog = 0, marksize = 8, tickfontsize = 20, imageFormat = :svg)
    crpt_num, crpt_den = criticalpoints(f) #get the critical pts
    crpt = sort(vcat(crpt_num, crpt_den)) #put them into a list and sort

    #If there are no crpts just plot the sign of the function at 0
    if length(crpt) == 0
        mdpts = 0
        signs = f(0) > 0 ? "+" : "-"
    #else get the midpts between crpts and the sign there
    else
        mdpts, signs = getsigns(f, crpt, values)
    end

    y = zeros(length(mdpts))
    gr()
    p = plot(mdpts, y, yaxis = false, label = label, legend = false, framestyle = :grid, fmt = imageFormat)
    p = yaxis!(label)
    p = yticks!([0],[""])
    p = plot!(size=size)
    p = plot!(ylims=(-.05,.6))
    p = plot!(xlims=(values[1]-horiz_jog, values[2]))
    p = hline!([0], linecolor = :black)
    #new
    ticklabels = [ @sprintf("%3.2f",convert(Float64,x)) for x in crpt ]
    crpt_float = [convert(Float64,x) for x in crpt]
    p = plot!(xticks=(crpt_float, ticklabels), tickfontsize=tickfontsize)
    #new
    #p = plot!(xlims=(values[1]-horiz_jog, values[2]), xticks=crpt, tickfontsize=tickfontsize)
    p = annotate!([(mdpts[i], .5, text(signs[i], :center, 30)) for i in 1:length(mdpts)])
    p = vline!(crpt, line = :dash, linecolor = :red)
    if !isempty(crpt_num)
        p = plot!(crpt_num, zeros(length(crpt_num)) .+ dotverticaljog, seriestype = :scatter, markercolor = :black, markersize = marksize)
    end
    if !isempty(crpt_den)
        p = plot!(crpt_den, zeros(length(crpt_den)) .+ dotverticaljog, seriestype = :scatter, markercolor = :white, markersize = marksize)
    end
    p
end
#=
function signchart(f, values, label="", horiz_jog = 0.2; size=(600, 200), dotverticaljog = 0, marksize = 8, tickfontsize = 20, imageFormat = :svg)
    crpt_num, crpt_den = criticalpoints(f) #get the critical pts
    crpt = sort(vcat(crpt_num, crpt_den)) #put them into a list and sort

    #If there are no crpts just plot the sign of the function at 0
    if length(crpt) == 0
        mdpts = 0
        signs = f(0) > 0 ? "+" : "-"
    #else get the midpts between crpts and the sign there
    else
        mdpts, signs = getsigns(f, crpt, values)
    end

    y = zeros(length(mdpts))
    gr()
    p = plot(mdpts, y, yaxis = false, label = label, legend = false, framestyle = :origin, fmt = imageFormat)
    p = yaxis!(label)
    p = plot!(size=size)
    p = plot!(ylims=(-.05,.6))
    #new
    ticklabels = [ @sprintf("%3.2f",x) for x in crpt ]
    p = plot!(xlims=(values[1]-horiz_jog, values[2]), xticks=(crpt, ticklabels), tickfontsize=tickfontsize)
    #new
    #p = plot!(xlims=(values[1]-horiz_jog, values[2]), xticks=crpt, tickfontsize=tickfontsize)
    p = annotate!([(mdpts[i], .5, text(signs[i], :center, 30)) for i in 1:length(mdpts)])
    p = vline!(crpt, line = :dash)
    if !isempty(crpt_num)
        p = plot!(crpt_num, zeros(length(crpt_num)) .+ dotverticaljog, seriestype = :scatter, markercolor = :black, markersize = marksize)
    end
    if !isempty(crpt_den)
        p = plot!(crpt_den, zeros(length(crpt_den)) .+ dotverticaljog, seriestype = :scatter, markercolor = :white, markersize = marksize)
    end
    p
end
=#

"""
    functionplot(values, f, horiz_ticks, label = ""; imageFormat = :svg)
"""
function functionplot(values, f, horiz_ticks, label = ""; imageFormat = :svg)
    #imageFormat can be :svg, :png, ...
    #p = plot(values, f, legend = :outertopright, framestyle = :origin, xticks=horiz_ticks, fmt = imageFormat)
    p = plot(values, f, legend = false, framestyle = :origin, xticks=horiz_ticks, fmt = imageFormat)
    p = yaxis!(label)
    if length(solve(simplify(1/f))) != 0
        p = vline!(solve(simplify(1/f)), line = :dash)
    end
    #Draw horizontal asymptote if it exists
    lf_lim=limit(f, x, -oo)
    if lf_lim == -oo || lf_lim == oo || lf_lim == -oo*1im || lf_lim == oo*1im
        lf_infty_bool = true
    else
        lf_infty_bool = false
    end
    rt_lim=limit(f, x, oo)
    if rt_lim == -oo || rt_lim == oo || lf_lim == -oo*1im || lf_lim == oo*1im
        rt_infty_bool = true
    else
        rt_infty_bool = false
    end

    #println("lf_infty_bool = $lf_infty_bool and rt_infty_bool = $rt_infty_bool")

    if !lf_infty_bool && !rt_infty_bool
        if lf_lim == rt_lim
            #println("lf_lim = rt_lim")
            p = hline!([lf_lim], line = :dash)
        else
            #println("lf_lim ≠ rt_lim")
            p = hline!([lf_lim], line = :dash)
            p = hline!([rt_lim],  line = :dash)
        end
    elseif !lf_infty_bool
        #println("lf_lim = $lf_lim")
        p = hline!([lf_lim], line = :dash)
    elseif !rt_infty_bool
        #println("rt_lim = $rt_lim")
        p = hline!([rt_lim],  line = :dash)
    end
    p
end

"""
    plot_function_sign_chart(f, fcn_width, horiz_ticks; labels = ["", "", ""], stacked::Bool=true, heights=[0.85 ,0.05, 0.05, 0.05], size=(900,800))
"""
function plot_function_sign_chart(f, fcn_width, horiz_ticks; labels = ["", "", ""], stacked::Bool=true, heights=[0.85 ,0.05, 0.05, 0.05], size=(500,800))
    s = signchart(q(x), fcn_width, labels[1])

    f′(x) = diff(f(x))
    s′ = signchart(f′(x), fcn_width, labels[2])

    f′′(x) = diff(f′(x))
    s′′ = signchart(f′′(x), fcn_width, labels[3])

    #p = functionplot(fcn_width[1]:.1:fcn_width[2], f(x), horiz_ticks[1]:horiz_ticks[2], labels[1])
    p = functionplot(fcn_width[1]:.1:fcn_width[2], f(x), horiz_ticks, labels[1])

    if stacked
        lay_out = grid(4,1, heights=heights)
    else
        lay_out = @layout[
            a{0.5w} grid(3,1)
        ]
    end
    plot(p, s, s′, s′′, layout = lay_out, size=size)
end

#summary code
abstract type Interval end
struct OpenOpenInterval <: Interval
    left
    right
end
struct OpenClosedInterval <: Interval
    left
    right
end
struct ClosedOpenInterval <: Interval
    left
    right
end
struct ClosedClosedInterval <: Interval
    left
    right
end

"""
    summary(f, values, domain::String = "(-∞, ∞)", label = ["y", "y′", "y′′"], fig_width = 200; dotverticaljog=0.05, marksize=10)

Takes a function f and outputs the:\n
y-intercept,\n
function sign chart, derivative sign chart, second derivative sign chart,\n
local max,\n
local min,\n
inflection points,\n
behaviour at infinites/endpoints.

The inputs are:\n
f: is the function to summarize,\n
values: min and max values for the sign charts,\n
domain: the domain of the function entered as a string (default "(-∞, ∞)"),\n
labels: is applied to the sign charts (default ["y", "y′", "y′′"]),\n
fig_width: gives the width of each sign chart (default 200),\n
dotverticaljog: changes the height of the points on the sign chart (default 0),\n
marksize: changes the diamater of the points on the sign chart (default 8),\n
tickfontsize: changes the font size of the tick marks (default 20),\n
format: displays the sighn charts either side-by-side or stacked on top of each other (default :side, other :stack).
"""
function summary(f, values, domain::String = "(-∞, ∞)", labels = ["y", "y′", "y′′"], fig_width = 200; dotverticaljog=0, marksize=8, tickfontsize = 20, format = :side, format_num="%3.2f")
    #Need to adjust values based on domain
    gr()
    interval = convert_to_interval(domain)
    y_intercept = f(0)
    y_int = "y-intercept: $(sprintf1(format_num, convert(Float64, y_intercept)))"

    if !isdir("figs")
        mkdir("figs")
    end
    f_name = string(f)
    f_name = replace(f_name, "*" => "")
    f_name = replace(f_name, "/" => "_")
    f′ = diff(f)
    f′′ = diff(f′)

    s = signchart(f, values, labels[1], dotverticaljog=dotverticaljog, marksize=marksize, tickfontsize=tickfontsize)
    filename = "./figs/sc"*f_name
    if isfile(filename*".png")
        rm(filename*".png")
    end
    savefig(filename)
    s′ = signchart(f′, values, labels[2], dotverticaljog=dotverticaljog, marksize=marksize, tickfontsize=tickfontsize)
    filename = "./figs/scp"*f_name
    if isfile(filename*".png")
        rm(filename*".png")
    end
    savefig(filename)
    s′′ = signchart(f′′, values, labels[3], dotverticaljog=dotverticaljog, marksize=marksize, tickfontsize=tickfontsize)
    filename = "./figs/scpp"*f_name
    if isfile(filename*".png")
        rm(filename*".png")
    end
    savefig(filename)

    f_min, f_max = extrema(f, domain)
    min_app = extrema_approx(f_min, true)
    max_app = extrema_approx(f_max, true)
    if f_min ==[]
        min_str = "local minimum: none"
    else
        min_str = replace("local minimum: $min_app", ")(" => "), (")
    end
    if f_min ==[]
        max_str = "local maximum: none"
    else
        max_str = replace("local maximum = $max_app", ")(" => "), (")
    end

    infpt = inflection_points(f, domain)
    infpt_app = extrema_approx(infpt, true)
    if infpt == []
        infp_str = "inflection point: none"
    else
        infp_str = replace("inflection point: $infpt_app", ")(" => "), (")
    end


    left_end, right_end = end_behavior(f, interval)

    if format == :side
        display("text/markdown", y_int*"<br>"
            #*"![](sc.png) ![](scp.png) ![](scpp.png)"
            *"""
            <table>
            <tr>
            <td><img src="./figs/sc$f_name.png" alt="drawing" width="$fig_width"/></td>
            <td><img src="./figs/scp$f_name.png" alt="drawing" width="$fig_width"/></td>
            <td><img src="./figs/scpp$f_name.png" alt="drawing" width="$fig_width"/></td>
            </tr>
            </table>
            """
            *"<br>"*min_str*"<br>"*max_str*"<br>"*infp_str*"<br>"*left_end*"<br>"*right_end)
    elseif format == :stack
        display("text/html", y_int*"<br>"
            #*"![](sc.png) ![](scp.png) ![](scpp.png)"
            *"""
            <img src="./figs/sc$f_name.png" alt="drawing" width="$fig_width"/><br>
            <img src="./figs/scp$f_name.png" alt="drawing" width="$fig_width"/><br>
            <img src="./figs/scpp$f_name.png" alt="drawing" width="$fig_width"/><br>
            """
            *"<br>"*min_str*"<br>"*max_str*"<br>"*infp_str*"<br>"*left_end*"<br>"*right_end)
    end

end

function convert_to_interval(domain::String)
    pattern = r"\s*(\[|\()\s*(-?[\p{Any}|\d]*)\s*,\s*(-?[\p{Any}|\d]*)\s*(\]|\))\s*"
    re_matched = match(pattern, domain)
    if re_matched[1] == "["
        if re_matched[4] == "]"
            left = tryparse(Float64, re_matched[2])
            right = tryparse(Float64, re_matched[3])
            domain = ClosedClosedInterval(left, right)
        elseif re_matched[4] == ")"
            left = tryparse(Float64, re_matched[2])
            right = re_matched[3] == "∞" ? oo : tryparse(Float64, re_matched[3])
            domain = ClosedOpenInterval(left, right)
        end
    elseif re_matched[1] == "("
        if re_matched[4] == "]"
            left = re_matched[2] == "-∞" ? -oo : tryparse(Float64, re_matched[2])
            right = tryparse(Float64, re_matched[3])
            domain = OpenClosedInterval(left, right)
        elseif re_matched[4] == ")"
            left = re_matched[2] == "-∞" ? -oo : tryparse(Float64, re_matched[2])
            right = re_matched[3] == "∞" ? oo : tryparse(Float64, re_matched[3])
            domain = OpenOpenInterval(left, right)
        end
    end
end

function extrema(f, domain::String = "(-∞, ∞)")
    f′ = diff(f)
    f′′ = diff(f′)
    crpt_num, crpt_den = criticalpoints(f′)
    second_derivative_at_crpt = []
    for x in crpt_num
        #push!(second_derivative_at_crpt, convert(Float32, f′′(x)))
        push!(second_derivative_at_crpt, (x,f′′(x)))
    end
    #add endpts if they exist
    interval = convert_to_interval(domain)
    left_end, right_end = end_behavior_aux(f, interval)
    max = []
    min = []
    #assuming the derivative is not = 0 at the endpoints, can add code later
    if f′(interval.left) > 0
        push!(min, (interval.left, left_end))
    elseif f′(interval.left) < 0
        push!(max, (interval.left, left_end))
    elseif f′(interval.left) == 0
        throw(DomainError(0, "This can't handle a slope of 0 at the endpoint yet."))
    end
    for (x, x′′) in second_derivative_at_crpt
        if x > interval.left && x < interval.right
            if x′′ < 0
                push!(max, (x,f(x)))
            elseif x′′ > 0
                push!(min, (x,f(x)))
            end
        end
    end
    if f′(interval.right) < 0
        push!(min, (interval.right, right_end))
    elseif f′(interval.right) > 0
        push!(max, (interval.right, right_end))
    elseif f′(interval.right) == 0
        throw(DomainError(0, "This can't handle a slope of 0 at the endpoints yet."))
    end
    min, max
end

function end_behavior_aux(f, interval::ClosedClosedInterval)
    f(interval.left), f(interval.right)
end
function end_behavior_aux(f, interval::ClosedOpenInterval)
    f(interval.left), missing
end
function end_behavior_aux(f, interval::OpenClosedInterval)
    missing, f(interval.right)
end
function end_behavior_aux(f, interval::OpenOpenInterval)
    missing, missing
end

function inflection_points(f, domain::String = "(-∞, ∞)")
    interval = convert_to_interval(domain)
    f′ = diff(f)
    f′′ = diff(f′)
    f′′′ = diff(f′′)
    crpt_num, crpt_den = criticalpoints(f′′)

    infpt = []
    for x in crpt_num
        if x > interval.left && x < interval.right
            if f′′′ ≠ 0
                push!(infpt, (x,f(x)))
            end
        end
    end

    infpt
end

function extrema_approx(ext, str = false)
    strg = ""
    for m in ext
        if !ismissing(m[2])
            if !str
                @printf("(%.3f, %.3f)", convert(Float64, m[1]), convert(Float64, m[2]))
            else
                strg *= @sprintf("(%.3f, %.3f)", convert(Float64, m[1]), convert(Float64, m[2]))
            end
        end
    end
    if str
        strg
    end
end

function end_behavior(f, interval::OpenOpenInterval)
    left_end_val = limit(f(x), x, interval.left)
    if interval.left == -oo
        if left_end_val == -oo
            left_end = "As x → -∞, " * "f(x) → -∞"
        elseif left_end_val == oo
            left_end = "As x → -∞, " * "f(x) → ∞"
        else
            left_end = "As x → -∞, " * "f(x) → $(left_end_val)"
        end
    else
        left_end = "As x → $(interval.left), " * "f(x) → $(left_end_val)"
    end

    right_end_val = limit(f(x), x, interval.right)
    if interval.right == oo
        if right_end_val == -oo
            right_end = "As x → ∞, " * "f(x) → -∞"
        elseif right_end_val == oo
            right_end = "As x → ∞, " * "f(x) → ∞"
        else
            right_end = "As x → ∞, " * "f(x) → $(right_end_val)"
        end
    else
        right_end = "As x → $(interval.right), " * "f(x) → $(right_end_val)"
    end
    left_end, right_end
end
function end_behavior(f, interval::OpenClosedInterval)
    if interval.left == -oo
        left_end_val = limit(f(x), x, -oo)
        if left_end_val == -oo
            left_end = "As x → -∞, " * "f(x) → -∞"
        elseif left_end_val == oo
            left_end = "As x → -∞, " * "f(x) → ∞"
        else
            left_end = "As x → -∞, " * "f(x) → $(left_end_val)"
        end
    else
        left_end = "At endpoint x = $(interval.left), " * "f($(interval.left)) = $(f(interval.left))"
    end

    right_end = "At right endpoint x = $(interval.right), " * "f($(interval.right)) = $(f(interval.right))"

    left_end, right_end
end
function end_behavior(f, interval::ClosedOpenInterval)
    left_end = "At left endpoint x = $(interval.left), " * "f($(interval.left)) = $(f(interval.left))"

    right_end_val = limit(f(x), x, interval.right)
    if interval.right == oo
        if right_end_val == -oo
            right_end = "As x → ∞, " * "f(x) → -∞"
        elseif right_end_val == oo
            right_end = "As x → ∞, " * "f(x) → ∞"
        else
            right_end = "As x → ∞, " * "f(x) → $(right_end_val)"
        end
    else
        right_end = "As x → $(interval.right), " * "f(x) → $(right_end_val)"
    end
    left_end, right_end
end
function end_behavior(f, interval::ClosedClosedInterval)
    left_end = "At left endpoint x = $(interval.left), " * "f($(interval.left)) = $(f(interval.left))"
    right_end = "At right endpoint x = $(interval.right), " * "f($(interval.right)) = $(f(interval.right))"

    left_end, right_end
end

function clean_figs_folder(arg::String)
    if arg == "all"
        file_list = readdir("./figs")
        for f in file_list
            f = "./figs/"*f
            rm(f)
        end
    else
        throw(DomainError(arg, """This has to be either "all" or the name of the function without the "(x)"\n
                                    e.g. clean_figs_folder(f)."""))
    end
end
function clean_figs_folder(f::Function)
    f_name = string(f(x))
    #println("f_name = $f_name")
    f_name = replace(f_name, "*" => "")
    #println("f_name = $f_name")
    file_list = readdir("./figs")
    file_del = ["sc"*f_name*".png" "scp"*f_name*".png" "scpp"*f_name*".png"]
    for f in file_list
        if f in file_del
            rm("./figs/"*f)
        end
    end
end

"""
    graph_f_and_derivative(values, f, horiz_ticks, label = ""; imageFormat = :svg, format=:single, size=(900,800))

Takes a function f and outputs the:\n
It's graph and the graph of the derivative of f′.

The inputs are:\n
values: values to be plotted (left:inc:right),\n
f: is the function to plot,\n
horiz_ticks: horizontal tick values (left:right),\n
label = "" (name of the function);\n
imageFormat = :svg (options :svg, :png, etc),\n
format=:single (:single = for both graphs on a single coordinate system
                :dual = side by side of graphs of f and f′),\n
size=(300,300), but for :dual if size = (x, y) then it will be (2x, y)
"""
function graph_f_and_derivative(values, f, horiz_ticks, label = ""; imageFormat = :svg, format=:single, size=(300,300))
    #imageFormat can be :svg, :png, ...
    f′ = diff(f(x))

    p = plot(values, f, label=label, legend = false, framestyle = :origin, xticks=horiz_ticks, fmt = imageFormat, color = :black)
    p = yaxis!(label)
    if format == :single
        if f.is_polynomial() && degree(f) == 1
            p = hline!([f′], line = :dash)
        else
            p = plot!(values, f′, linestyle = :dash)
        end
        p
    elseif format == :dual
        fst, lst = split(label, "(")
        label = fst*"′("*lst
        if f.is_polynomial() && degree(f) == 1
            #println("f.is_polynomial = $(f.is_polynomial()) and degree(f) = $degree(f)")
            fcn_values = fill(f′, size(values)[1])
            q = plot(values, fcn_values, legend = false, framestyle = :origin, xticks=horiz_ticks, fmt = imageFormat, color = :red, linestyle = :dash)
            q = yaxis!(label)
        else
            #println("f.is_polynomial = $(f.is_polynomial()) and degree(f) = $degree(f)")
            q = plot(values, f′, legend = false, framestyle = :origin, xticks=horiz_ticks, fmt = imageFormat, color = :red, linestyle = :dash)
            q = yaxis!(label)
        end

        size_dual =(size[1]*2, size[2])
        plot(p, q, layout = (1,2), size=size_dual)
    end
end
