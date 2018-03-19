# module containing patterns used for Metamaterial waveguide structure
module MetamaterialWaveguide
    using Devices, Clipper
    import Devices: Meta, GDSMeta
    export cross_polygon, mwg_unit_cell!, mwg_last_cell!, mwg_input_cap!, mwg_cell_with_claw!

    """
        cross_polygon(hline_width, hline_height, vline_width,
        vline_height; center=true)
    Create a cross-shaped polygon of specified widths and heights of horizontal
    line and vertical line. The keyword argument `center=true` determines whether
    the returned polygon is centered at the origin or at distance hline_width/2
    along the x-axis away from the origin.
    """
    function cross_polygon(hline_width, hline_height, vline_width,
        vline_height; center=true)

        if hline_width < vline_width
            error("The width of the horizontal line must be larger than that of
                the vertical line")
        elseif hline_height > vline_height
            error("The height of the vertical line must be larger than that of
                the horizontal line")
        end

        offset =  (~center) * Point(hline_width/2, typeof(hline_width)(0))
        hline = centered(Rectangle(hline_width, hline_height)) + offset
        vline = centered(Rectangle(vline_width, vline_height)) + offset

        poly = clip(Clipper.ClipTypeUnion, hline, vline)[1]
        return poly
    end

    """
        mwg_unit_cell!{T,S}(c::Cell{T,S}, Lr_width, Lr_height, Cr_overlap,
        res_ground_length, res_conductor_width, res_conductor_gap,
        Ck_terminal_width, Ck_overlap, wg_meander_length, wg_meander_nlegs::Int,
        wg_ground_gap, wg_conductor_width, wg_conductor_gap, centerline_width,
        centerline_length, meta=GDSMeta(0,0))
    Create a unit cell of the metamaterial waveguide.
    - `c` : a cell to render metamaterial waveguide in
    - `Lr_width` : width of the rectangular coil Lr
    - `Lr_height` : height of the rectangular coil Lr
    - `Cr_overlap` : overlap length of the interdigital capacitor Cr forming
        the resonator
    - `res_ground_length` : length of the line from the top of the resonator
        to the ground. This also sets the gap between the resonator
        and the ground.
    - `res_conductor_width` : width of the conductor forming a resonator
    - `res_conductor_gap` : gap between conductor traces forming a resonator
    - `Ck_terminal_width` : terminal width of the coupling capacitor Ck between
        the resonator and the centerline
    - `Ck_overlap` : overlap length of the coupling capacitor Ck
    - `wg_meander_length` : the length of meandering of waveguide. The trace of
        the meander goes back and forth vertically with this distance.
    - `wg_meander_nlegs` : number of legs in the meander forming the waveguide
    - `wg_ground_gap` : the distance between the meander trace of the waveguide
        and the ground
    - `wg_conductor_width` : width of the conductor forming the meander
    - `wg_conductor_gap` : gap between conductor traces forming the meander
    - `centerline_width` : width of the centerline that connects the resonator
        with waveguide.
    - `centerline_length` : length of the centerline. This has to be equal to or
        greater than `Lr_width`
    """
    function mwg_unit_cell!(c::Cell{T,S}, Lr_width, Lr_height, Cr_overlap,
        res_ground_length, res_conductor_width, res_conductor_gap, Ck_terminal_width,
        Ck_overlap, wg_meander_length, wg_meander_nlegs::Int, wg_ground_gap,
        wg_conductor_width, wg_conductor_gap, centerline_width, centerline_length,
        meta=GDSMeta(0,0)) where {T,S}
        if Lr_width > centerline_length
            error("The length of the centerline must be larger than the width of the coil")
        end
        wg_length = (wg_meander_nlegs - 1) * (wg_conductor_width + wg_conductor_gap) + wg_conductor_width / 2
        res_height = 2 * (Lr_height + res_conductor_gap) + Cr_overlap

        # define a mask for unit cell
        wg_mask_width = 2 * wg_length + centerline_length
        wg_mask_height = wg_meander_length + 2 * (wg_conductor_width + wg_ground_gap)

        res_mask_width = Lr_width + 2 * res_ground_length
        res_mask_height = centerline_width + 2 * (res_ground_length + res_height + Ck_terminal_width + Ck_overlap + 3 * res_conductor_gap)

        # global mask = res_mask ∪ wg_mask
        # mask = Cell{T,S}(uniquename("mask"))
        # mask_poly = cross_polygon(wg_mask_width, wg_mask_height, res_mask_width,
            # res_mask_height; center=false)
        # render!(mask, mask_poly)

        # lumped-element resonator with coupling capacitor
        res_with_Ck = Cell{T,S}(uniquename("lumped_elem_resonator_with_Ck"))
        lumped_elem_resonator_with_Ck!(res_with_Ck, Lr_width, Lr_height, Cr_overlap, Ck_overlap,
            Ck_terminal_width, res_conductor_width, res_conductor_gap, meta,
            ground_length=res_ground_length, Ck_terminal_width2=zero(T))

        res = Cell{T,S}(uniquename("resonator"))
        push!(res.refs, CellReference(res_with_Ck, Point(wg_mask_width/2, centerline_width/2)))
        push!(res.refs, CellReference(res_with_Ck,
            Point(wg_mask_width/2, -centerline_width/2), xrefl=true))

        # meander forming waveguide part of MWG
        meander = Cell{T,S}(uniquename("meander"))
        rectangular_meander!(meander, wg_meander_length, wg_meander_nlegs,
            wg_conductor_width, wg_conductor_gap, meta;
            first_length=wg_meander_length + wg_conductor_width,
            last_length=(wg_meander_length + centerline_width)/2,
            start_upward=true)

        wg = Cell{T,S}(uniquename("waveguide"))
        push!(wg.refs, CellReference(meander, Point(-wg_conductor_width/2,
            -wg_conductor_width-wg_meander_length/2)))
        push!(wg.refs, CellReference(meander,
            Point(wg_conductor_width/2+ wg_mask_width,
            wg_conductor_width + wg_meander_length/2), rot=π))

        # centerline
        centerline = Cell{T,S}(uniquename("centerline"))
        render!(centerline, centered(Rectangle(centerline_length, centerline_width)) +
            Point(wg_mask_width/2, zero(T)), meta)

        # # take the negative trace of each element
        # res_neg = negative(res, mask)
        # wg_neg = negative(wg, mask)
        # centerline_neg = negative(centerline, mask)

        # take the intersection res_neg ∩ wg_neg ∩ centerline_neg
        # unit_cell_poly = clip(Clipper.ClipTypeIntersection,
        #     polygon.(elements(res_neg)), polygon.(elements(wg_neg)))
        # unit_cell_poly = clip(Clipper.ClipTypeIntersection,
        #     unit_cell_poly, polygon.(elements(centerline_neg)))
        #
        # for p in unit_cell_poly
        #     render!(c, p, meta)
        # end
        for part in [res, wg, centerline]
            push!(c.refs, CellReference(part, Point(zero(T), zero(T))))
        end
        return c
    end

    """
    Metamaterial waveguide cell with claw coming out from the resonator node
    """
    function mwg_cell_with_claw!(c::Cell{T,S}, Lr_width, Lr_height, Cr_overlap,
        res_ground_length, res_conductor_width, res_conductor_gap, Ck_terminal_width,
        Ck_overlap, wg_meander_length, wg_meander_nlegs::Int, wg_ground_gap,
        wg_conductor_width, wg_conductor_gap, centerline_width, centerline_length,
        res_claw_gap, res_claw_conductor_width, claw_width, claw_length,
        claw_ground_gap, claw_dist, meta::Meta=GDSMeta(0,0); rot=0) where {T,S}

        wg_length = ((wg_meander_nlegs - 1) * (wg_conductor_width + wg_conductor_gap)
            + wg_conductor_width/2)
        res_height = 2 * (Lr_height + res_conductor_gap) + Cr_overlap

        # define a mask for unit cell
        wg_mask_width = 2 * wg_length + centerline_length
        wg_mask_height = wg_meander_length + 2 * (wg_conductor_width + wg_ground_gap)

        res_mask_width = Lr_width + 2 * res_ground_length
        res_mask_height = centerline_width + 2 * (res_ground_length + res_height
            + Ck_terminal_width + Ck_overlap + 3 * res_conductor_gap)

        claw_mask_left = wg_mask_width/2 + Lr_width/2
        claw_mask_right = (claw_mask_left + res_claw_gap + res_claw_conductor_width
            + claw_length + claw_ground_gap)
        claw_mask_vdist = claw_dist/2 - claw_ground_gap
        claw_mask_hdist = res_claw_gap + res_claw_conductor_width + claw_ground_gap

        wg_mask = Rectangle(Point(zero(T), -wg_mask_height/2),
            Point(wg_mask_width/2, wg_mask_height/2))
        res_mask = centered(Rectangle(res_mask_width, res_mask_height)) + Point(wg_mask_width/2, zero(T))

        claw_mask_bulk = Rectangle(Point(claw_mask_left, -res_mask_height/2),
            Point(claw_mask_right, res_mask_height/2))
        claw_mask_rm = Rectangle(Point(claw_mask_left + claw_mask_hdist,
            -claw_mask_vdist), Point(claw_mask_right, claw_mask_vdist))
        claw_mask = clip(Clipper.ClipTypeDifference, claw_mask_bulk, claw_mask_rm)[1]

        # global mask = res_mask ∪ wg_mask ∪ claw_mask
        # mask = Cell{T,S}(uniquename("mask"))
        # mask_poly = clip(Clipper.ClipTypeUnion, wg_mask, res_mask)[1]
        # mask_poly = clip(Clipper.ClipTypeUnion, mask_poly, claw_mask)[1]
        # render!(mask, mask_poly, meta_mask)

        # lumped-element resonator with coupling capacitor
        res_with_Ck = Cell{T,S}(uniquename("lumped_elem_resonator_with_Ck"))
        lumped_elem_resonator_with_Ck!(res_with_Ck, Lr_width, Lr_height, Cr_overlap, Ck_overlap,
        Ck_terminal_width, res_conductor_width, res_conductor_gap, meta,
        ground_length=res_ground_length, Ck_terminal_width2=zero(T))

        res = Cell{T,S}(uniquename("resonator"))
        push!(res.refs, CellReference(res_with_Ck,
            Point(wg_mask_width/2, centerline_width/2)))
        push!(res.refs, CellReference(res_with_Ck,
            Point(wg_mask_width/2, -centerline_width/2), xrefl=true))

            # meander forming waveguide part of MWG
        meander = Cell{T,S}(uniquename("meander"))
        rectangular_meander!(meander, wg_meander_length, wg_meander_nlegs,
            wg_conductor_width, wg_conductor_gap, meta;
            first_length=wg_meander_length + wg_conductor_width,
            last_length=(wg_meander_length + centerline_width)/2,
            start_upward=true)

        wg = Cell{T,S}(uniquename("waveguide"))
        push!(wg.refs, CellReference(meander, Point(-wg_conductor_width/2,
            -wg_conductor_width-wg_meander_length/2)))

        # centerline
        centerline = Cell{T,S}(uniquename("centerline"))
        render!(centerline, centered(Rectangle(centerline_length, centerline_width)) +
            Point(wg_mask_width/2, zero(T)), meta)

        # claw
        claw = Cell{T,S}(uniquename("claw"))
        res_claw_conn_vdist = (res_mask_height/2 - Lr_height - res_ground_length
            - 2 * res_conductor_gap - Cr_overlap)
        render!(claw, Rectangle(Point(claw_mask_left, res_claw_conn_vdist - 2*res_conductor_width),
            Point(claw_mask_left + res_claw_gap, res_claw_conn_vdist)), meta)
        render!(claw, Rectangle(Point(claw_mask_left, -res_claw_conn_vdist + 2*res_conductor_width),
            Point(claw_mask_left + res_claw_gap, -res_claw_conn_vdist)), meta)

            claw_left = claw_mask_left + res_claw_gap + res_claw_conductor_width
            claw_right = claw_left + claw_length

        render!(claw, Rectangle(Point(claw_left - res_claw_conductor_width,
            -res_claw_conn_vdist), Point(claw_left, res_claw_conn_vdist)), meta)

        render!(claw, Rectangle(Point(claw_left, claw_mask_vdist + claw_ground_gap),
            Point(claw_right, claw_mask_vdist + claw_width + claw_ground_gap)), meta)
        render!(claw, Rectangle(Point(claw_left, -claw_mask_vdist - claw_ground_gap),
            Point(claw_right, -claw_mask_vdist - claw_width - claw_ground_gap)), meta)

        # This part discarded due to instability of Clipper
        # take the negative trace of each element
        # res_neg = negative(res, mask)
        # wg_neg = negative(wg, mask)
        # centerline_neg = negative(centerline, mask)
        # claw_neg = negative(claw, mask)
        # take the intersection res_neg ∩ wg_neg ∩ centerline_neg
        # cell_poly = clip(Clipper.ClipTypeIntersection,
        #     polygon.(elements(res_neg)), polygon.(elements(wg_neg)))
        # cell_poly = clip(Clipper.ClipTypeIntersection,
        #     cell_poly, polygon.(elements(centerline_neg)))
        # cell_poly = clip(Clipper.ClipTypeIntersection,
        #     cell_poly, polygon.(elements(claw_neg)))
        # for p in cell_poly
            # render!(c, p, meta)
        # end

        for part in [res, wg, centerline, claw]# ,mask]
            push!(c.refs, CellReference(part, Point(zero(T), zero(T)), rot=rot))
        end

        # render!(c, mask, meta_mask)
        return c
    end

    #
    # """
    # First cell for the metamaterial waveguide including the input capacitor
    # """
    # function mwg_first_cell!(c::Cell{T,S}, Lr_width, Lr_height, Cr_overlap_first,
    #     res_ground_length, res_conductor_width, res_conductor_gap, Ck_terminal_width,
    #     Ck_overlap, wg_meander_length, wg_meander_nlegs::Int, wg_ground_gap,
    #     wg_conductor_width, wg_conductor_gap, centerline_width, centerline_length,
    #     res_claw_gap, res_claw_conductor_width_first, claw_width_first,
    #     claw_length_first, claw_ground_gap_first, input_cap_width, input_cap_gap, ground_width_first,
    #     meta::Meta=GDSMeta(0,0)) where {T,S}
    #     tmp = Cell{T,S}(uniquename("last_cell"))
    #     mwg_last_cell!(tmp, Lr_width, Lr_height, Cr_overlap_first,
    #         res_ground_length, res_conductor_width, res_conductor_gap,
    #         Ck_terminal_width, Ck_overlap, wg_meander_length, wg_meander_nlegs,
    #         wg_ground_gap, wg_conductor_width, wg_conductor_gap,
    #         centerline_width, centerline_length, res_claw_gap,
    #         res_claw_conductor_width_first, claw_width_first, claw_length_first,
    #         claw_ground_gap_first, input_cap_width, input_cap_gap,
    #         ground_width_first, meta)
    #     for p in polygon.(tmp.elements)
    #         render!(c, Rotation(π)(p)+Point(width(bounds(tmp)), zero(T)), meta)
    #     end
    #     return c
    # end

    # """
    # Input capacitor for the metamaterial waveguide
    # """
    # function mwg_input_cap!(c::Cell{T,S}, res_mask_width, res_mask_height,
    #     wg_mask_width, wg_mask_height, wg_meander_length, wg_meander_nlegs::Int,
    #     wg_ground_gap, wg_conductor_width, wg_conductor_gap, centerline_width,
    #     input_cap_terminal_width, input_cap_terminal_width2, input_cap_overlap,
    #     input_cap_finger_width, input_cap_finger_spacing, input_cap_end_gap,
    #     cpw_gap, cpw_trace, meta=GDSMeta(0,0)) where {T,S}
    #
    #     wg_length = (wg_meander_nlegs - 1) * (wg_conductor_width + wg_conductor_gap) + wg_conductor_width/2
    #
    #     input_cap_terminal_length = 2 * wg_conductor_width + wg_meander_length
    #
    #     # define a mask for input capacitor
    #     wg_mask = (centered(Rectangle(wg_mask_width/2, wg_mask_height))
    #         + Point((wg_mask_width/2+res_mask_width)/2, zero(T)))
    #     res_mask = (centered(Rectangle(res_mask_width, res_mask_height))
    #         + Point(res_mask_width/2, zero(T)))
    #
    #     # global mask = res_mask ∪ wg_mask
    #     mask = clip(Clipper.ClipTypeUnion, wg_mask, res_mask)[1]
    #
    #     # meander forming waveguide part of MWG
    #     meander = Cell{T,S}(uniquename("meander"))
    #     rectangular_meander!(meander, wg_meander_length, wg_meander_nlegs,
    #         wg_conductor_width, wg_conductor_gap;
    #         first_length=wg_meander_length + wg_conductor_width,
    #         last_length=(wg_meander_length + centerline_width)/2,
    #         start_upward=true)
    #     wg = Cell{T,S}(uniquename("waveguide"))
    #     push!(wg.refs, CellReference(meander,
    #         Point(wg_conductor_width/2 + wg_mask_width/2 + res_mask_width/2,
    #         wg_conductor_width + wg_meander_length/2), rot=π))
    #
    #     # input capacitor
    #     cap = Cell{T,S}("capacitor")
    #     interdigit = Cell{T,S}("interdigit")
    #     interdigit_capacitor!(interdigit, input_cap_terminal_width, input_cap_terminal_length,
    #         input_cap_overlap, input_cap_finger_width, input_cap_finger_spacing, input_cap_end_gap,
    #         terminal_width2=input_cap_terminal_width2)
    #     push!(cap.refs, CellReference(interdigit, Point(zero(T), -input_cap_terminal_length/2)))
    #
    #     # centerline
    #     conn_line = Cell{T,S}(uniquename("connecting_line"))
    #     render!(conn_line, Rectangle(Point(width(bounds(interdigit)), -centerline_width/2),
    #         Point(bounds(wg).ll.x, centerline_width/2)), meta)
    #     push!(cap.refs, CellReference(conn_line, Point(zero(T), zero(T))))
    #
    #     input_cap_neg = negative(cap, mask)
    #     input_wg_neg = negative(wg, mask)
    #
    #     # take the intersection res_neg ∩ wg_neg ∩ centerline_neg
    #     input_cap_poly = clip(Clipper.ClipTypeIntersection,
    #         polygon.(elements(input_cap_neg)), polygon.(elements(input_wg_neg)))
    #
    #     input_path = Path(Point(typeof(input_cap_terminal_length)(0),
    #         typeof(input_cap_terminal_length)(0)), α0 = 0)
    #     input_taper_len = (wg_mask_width - res_mask_width)/2
    #     input_trace = input_cap_terminal_length
    #     input_gap = input_trace * cpw_gap/cpw_trace
    #
    #     straight!(input_path, input_taper_len, Paths.CPW(
    #         t->(cpw_trace * (input_taper_len - t)/input_taper_len + t/input_taper_len * input_trace),
    #         t->(cpw_gap * (input_taper_len - t)/input_taper_len + t/input_taper_len * input_gap)))
    #
    #     render!(c, input_path, meta)
    #     trans = Translation(Point(input_taper_len, zero(T)))
    #     for p in input_cap_poly
    #         render!(c, trans(p), meta)
    #     end
    #     return c
    # end
end

# package code goes here
