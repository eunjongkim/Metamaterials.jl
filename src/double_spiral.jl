

"""
    double_spiral(width, height, conductor_width,
        conductor_gap, clockwise, α0, in_length, out_length)
Create a positive trace of rectangular double spiral coil of specified lengths,
conductor width, and conductor gap. The keyword arguments `clockwise` and
`α0` defines the direction of winding and the angle, respectively.
"""
function double_spiral(width, height, conductor_width,
    conductor_gap, clockwise, α0, in_length, out_length)

    minimal_length = 3 * conductor_width + 2 * conductor_gap
    if (height < minimal_length) | (width < minimal_length)
        error("""
        Cannot make a rectangular double coil with given parameters.
        Both the width and the height of the coil must be larger than
        3*conductor_width + 2*conductor_gap
        """)
    end
    # turn_direction:
    # -1 if?lum clockwise-winding coil, +1 if counter-clockwise-winding coil
    turn_direction = -2 * clockwise + 1

    point0 = -turn_direction * Point(conductor_width/2 * sin(α0),
        - conductor_width/2 * cos(α0))
    p = Path(point0, α0=α0)

    path_style = Paths.Trace(conductor_width)
    corner_style = Paths.SimpleTraceCorner()

    # first length and second length
    l1 = height - conductor_width / 2
    l2 = width - conductor_gap - 2 * conductor_width
    # create a list of lengths of straight lines
    l_list = [l1, l2]
    l = l1 - (1.5 * conductor_width + conductor_gap)
    while l >= (conductor_width + conductor_gap)
        push!(l_list, l)
        l = l_list[end-1] - 2 * (conductor_width + conductor_gap)
    end

    # length of the input line
    l_list[1] += in_length
    # scale the innermost lengths to fit nice and symmetrically
    l_list[end] = l_list[end] / 2 + (conductor_width + conductor_gap) / 2
    # middle section
    l_mid = l_list[end-1] - conductor_width - conductor_gap

    # wind in
    for n in 1:length(l_list)
        straight!(p, l_list[n], path_style)
        corner!(p, turn_direction * π/2, corner_style)
    end
    # connect (middle segment)
    if l_mid > 0 * l_mid
        straight!(p, l_mid, path_style)
    end

    l_list[1] -= in_length
    # length of the output line
    l_list[1] += out_length
    # wind out
    for n in length(l_list):-1:1
        corner!(p, -turn_direction * π/2, corner_style)
        straight!(p, l_list[n], path_style)
    end
    info("Double-spiral coil pathlength : $(pathlength(p) - in_length - out_length)")
    return p
end


"""
    double_spiral!{T,S}(cell::Cell{T,S}, width, height, conductor_width,
        conductor_gap, meta=GDSMeta(0,0); clockwise=true, α0=π/2)
Create a positive trace of rectangular double spiral coil of specified lengths,
conductor width, and conductor gap. The keyword arguments `clockwise` and
`α0` defines the direction of winding and the angle, respectively.
"""
function double_spiral!(c::Cell{T,S}, width, height, conductor_width,
    conductor_gap, meta::Meta=GDSMeta(0,0); clockwise=true, α0=π/2,
    in_length=zero(T), out_length=zero(T)) where {T,S}

    p = double_spiral(width, height, conductor_width,
        conductor_gap, clockwise, α0, in_length, out_length)
    render!(c, p, meta)
    return c
end
