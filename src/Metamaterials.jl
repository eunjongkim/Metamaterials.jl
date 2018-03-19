module Metamaterials

using Devices, Clipper
import Devices: Meta, GDSMeta, AbstractPolygon, Coordinate

export double_spiral!, rectangular_meander!, interdigit_capacitor!,
    lumped_elem_resonator!, lumped_elem_resonator_with_Ck!, negative,
    cross_polygon

include("double_spiral.jl")
include("interdigit_capacitor.jl")
include("rectangular_meander.jl")

"""
    lumped_elem_resonator!{T,S}(c::Cell{T,S}, Lr_width, Lr_height,
        Cr_overlap, conductor_width, conductor_gap, meta::Meta=GDSMeta(0,0);
        in_length=zero(T), out_length=zero(T))
Create a positive trace of symmetric lumped-element resonator with
specified double spiral coil dimensions (`Lr_width`, `Lr_height`), capacitor
overlap (`Cr_overlap`). The conductor properties are defined with
`conductor_width` and `conductor_gap`. The length of the conductor that goes
in and out of the resonator can be specified by keyword arguments
`in_length` and `out_length`.
"""
function lumped_elem_resonator!(c::Cell{T,S}, Lr_width, Lr_height, Cr_overlap,
    conductor_width, conductor_gap, meta::Meta=GDSMeta(0,0);
    in_length=zero(T), out_length=zero(T)) where {T,S}

    # width and height of the resonator as a whole
    width, height = Lr_width, 2 * Lr_height + Cr_overlap + conductor_gap

    # constants for rectangular coils
    connect_length = Cr_overlap / 2 + conductor_gap

    # constants for interdigitial capacitor
    terminal_width = zero(T)
    terminal_length = Lr_width - conductor_gap - conductor_width

    # cells
    upper_coil = Cell{T,S}(uniquename("upper_coil_Lr"))
    double_spiral!(upper_coil, Lr_width, Lr_height, conductor_width,
        conductor_gap, meta, in_length=connect_length, out_length=out_length)

    lower_coil = Cell{T,S}(uniquename("lower_coil_Lr"))
    double_spiral!(lower_coil, Lr_width, Lr_height, conductor_width,
        conductor_gap, meta, clockwise=false, in_length=connect_length,
        out_length=in_length)

    cap = Cell{T,S}(uniquename("capacitor_Cr"))
    interdigit_capacitor!(cap, terminal_width, terminal_length, Cr_overlap,
        conductor_width, conductor_gap, conductor_gap, meta)

    # rendering cells
    pointoffset = Point(-width, Lr_height + connect_length + in_length)

    trans_upper_coil = Translation(pointoffset...)
    upper_coil_polys = [trans_upper_coil(p) for
        p in polygon.(elements(upper_coil))]

    trans_lower_coil = Translation(pointoffset...) ∘ Rotation(π)
    lower_coil_polys = [trans_lower_coil(p) for
        p in polygon.(elements(lower_coil))]

    capoffset = pointoffset + Point(conductor_width + conductor_gap,
        connect_length)
    trans_cap = Translation(capoffset...) ∘ Rotation(-π/2)
    cap_polys = [trans_cap(p) for p in polygon.(elements(cap))]

    for p in upper_coil_polys
        render!(c, p, meta)
    end
    for p in cap_polys
        render!(c, p, meta)
    end
    for p in lower_coil_polys
        render!(c, p, meta)
    end
    return c
end

    """ lumped_elem_resonator_with_Ck!{T,S}(c::Cell{T,S}, Lr_width, Lr_height,
            Cr_overlap, Ck_overlap, Ck_terminal_width, conductor_width,
            conductor_gap, meta::Meta=GDSMeta(0,0); ground_length=zero(T),
            Ck_terminal_width2=Ck_terminal_width)
    Create a positive trace of symmetric lumped-element resonator with coupling
    capacitor. The coil dimensions (`Lr_width`, `Lr_height`) and capacitor
    overlap (`Cr_overlap`) defines the resonator. The terminal width
    (`Ck_terminal_width`) and the overlap (`Ck_overlap`) characterizes the
    coupling capacitor. The conductor properties are defined with
    `conductor_width` and `conductor_gap`. The width of the terminal from the
    origin and the length of the conductor into the ground can be specified with
    keyword arguments `Ck_terminal_width2` and `ground_length`, respectively.
    """
    function lumped_elem_resonator_with_Ck!(c::Cell{T,S}, Lr_width, Lr_height,
        Cr_overlap, Ck_overlap, Ck_terminal_width, conductor_width,
        conductor_gap, meta::Meta=GDSMeta(0,0); ground_length=zero(T),
        Ck_terminal_width2=Ck_terminal_width) where {T,S}

        # resonator
        res = Cell{T,S}(uniquename("lumped_elem_res_LrCr"))
        lumped_elem_resonator!(res, Lr_width, Lr_height, Cr_overlap,
            conductor_width, conductor_gap, meta,
            in_length=Ck_terminal_width+conductor_gap,
            out_length=ground_length)

        # capacitor
        cap = Cell{T,S}(uniquename("coupling_cap_Ck"))
        Ck_terminal_length = Lr_width - 2 * conductor_width
        interdigit_capacitor!(cap, Ck_terminal_width, Ck_terminal_length,
            Ck_overlap, conductor_width, conductor_gap, conductor_gap, meta,
            terminal_width2=Ck_terminal_width2)

        # offset point
        p_x = -Ck_terminal_length/2
        p_y = Ck_overlap+2*conductor_gap+Ck_terminal_width+Ck_terminal_width2
        pointoffset = Point(p_x, p_y)

        Ck_ref = CellReference(cap, pointoffset, rot=-π/2)
        res_ref = CellReference(res,
        Point(Lr_width/2, Ck_overlap + 2 * conductor_gap + Ck_terminal_width2))

        push!(c.refs, Ck_ref)
        push!(c.refs, res_ref)
        return c
    end

    """
        negative{T,S}(c::Cell{T,S}, mask::AbstractPolygon,
            meta::Meta=GDSMeta(0,0))
    Take whatever is drawn in the cell off from the mask polygon. i.e. create a
    negative cell drawing within the extent of mask. The references of the
    returned cell is discarded and only the flattened cell remains.
    """
    function negative(c::Cell{T,S}, mask::AbstractPolygon{U},
        meta::Meta=GDSMeta(0,0)) where {T,S,U<:Coordinate}
        c = flatten!(c)
        poly = polygon.(elements(c))

        neg = clip(Clipper.ClipTypeDifference, [mask], poly)[1]
        cell = Cell{T,S}(uniquename("negative_$(name(c))"))
        render!(cell, neg, meta)
        return cell
    end
    function negative(c::Cell{T,S}, mask::AbstractVector{U},
        meta::Meta=GDSMeta(0,0)) where {T,S,U<:AbstractPolygon}
        c = flatten!(c)
        poly = polygon.(elements(c))

        neg = clip(Clipper.ClipTypeDifference, mask, poly)[1]
        cell = Cell{T,S}(uniquename("negative_$(name(c))"))
        render!(cell, neg, meta)
        return cell
    end


end # module
