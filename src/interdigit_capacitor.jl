# code for drawing interdigit capacitor
"""
    interdigit_capacitor!{T,S}(c::Cell{T,S}, terminal_width, terminal_length,
        overlap, finger_width, finger_spacing, end_gap,
        meta::Meta=GDSMeta(0,0); terminal_width2=terminal_width)
Create a positive trace of an interdigital capacitor of a specified
terminal width, terminal length, finger width, finger spacing, and end gap.
The width of the second terminal can be asymetrically defined with the
keyword argument `terminal_width2`.
"""
function interdigit_capacitor!(c::Cell{T,S}, terminal_width, terminal_length,
    overlap, finger_width, finger_spacing, end_gap, meta::Meta=GDSMeta(0,0);
    terminal_width2=terminal_width) where {T<:Coordinate,S<:GDSMeta}

    # terminal
    r1 = Rectangle(terminal_width, terminal_length)
    r2 = Rectangle(terminal_width2, terminal_length)

    # fingers
    fingers = Cell{T,S}(uniquename("fingers"))
    finger_length = overlap + end_gap

    remainder = rem(terminal_length, 2 * (finger_width + finger_spacing))

    # criterion to determine npairs and skiplast
    if remainder < finger_width
        npairs = Int(fld(terminal_length, 2 * (finger_width + finger_spacing)))
        skiplast = 0
    elseif remainder < 2 * finger_width + finger_spacing
        npairs = Int(fld(terminal_length,
            2 * (finger_width + finger_spacing)) + 1)
        skiplast = 1
    else
        npairs = Int(fld(terminal_length,
            2 * (finger_width + finger_spacing)) + 1)
        skiplast = 0
    end
    interdigit!(fingers, finger_width, finger_length, finger_spacing,
        end_gap, npairs, skiplast, meta)

    h = height(bounds(fingers))
    fingers_offset = (terminal_length - h)/2

    trans = Translation(terminal_width, fingers_offset)
    polys = [trans(p) for p in polygon.(elements(fingers))]
    for p in polys
        render!(c, p, meta)
    end
    render!(c, r1, meta)
    render!(c, r2 + Point(terminal_width + end_gap + finger_length,
        zero(T)), meta)
    return c
end
