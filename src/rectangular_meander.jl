"""
    rectangular_meander!(c::Cell{T,S}, length, nlegs::Int,
        conductor_width, conductor_spacing, meta::Meta=GDSMeta(0,0);
        first_length=length/2, last_length=length/2, start_upward=true)
Create a positive trace of rectangular meander specified by its
vertical length (`length`), number of legs (`nlegs::Int`), conductor
width and spacing (`conductor_width` and `conductor_spacing`). The
first and last length are initialized as half the `length`, and can
be specified with keyword arguments `first_length` and `last_length`.
The meander starts upward unless the keyword argument
`start_upward' is enforced to `false`.
"""
function rectangular_meander!(c::Cell{T,S}, length, nlegs::Int,
    conductor_width, conductor_spacing, meta::Meta=GDSMeta(0,0);
    first_length=length/2, last_length=length/2,
    start_upward=true) where {T<:Coordinate, S<:GDSMeta}
    if nlegs < 2
        error("The number of legs must be 2 or larger.")
    end

    # starting direction:
    # +1 if start_upward=true, -1 if start_upward=false
    dir = 2 * start_upward - 1

    p = Path(Point(conductor_width/2, zero(T)), α0=dir*π/2)

    path_style = Paths.Trace(conductor_width)
    corner_style = Paths.SimpleTraceCorner()

    straight!(p, first_length + conductor_width/2, path_style)
    corner!(p, -dir * π/2, corner_style)
    straight!(p, (conductor_spacing + conductor_width)/2, path_style)
    for n in 1:(nlegs - 2)
        straight!(p, (conductor_spacing + conductor_width)/2, path_style)
        corner!(p, dir * (-1)^n * π/2, corner_style)
        straight!(p, conductor_width + length, path_style)
        corner!(p, dir * (-1)^n * -π/2, corner_style)
        straight!(p, (conductor_spacing + conductor_width)/2, path_style)
    end
    straight!(p, (conductor_spacing + conductor_width)/2, path_style)
    corner!(p, dir * (-1)^(nlegs-1) * π/2, corner_style)
    straight!(p, conductor_width/2 + last_length, path_style)

    render!(c, p, meta)
end
