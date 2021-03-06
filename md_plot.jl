folder = "result_4_10_2/"
# folder = ""
# idx = 0
T = 350
@unpack rp, sp, Vp = rd;
# for idx in LinRange(25, T, Integer(T/25))
for idx in [400]
    tag = string(folder, "malpasset_h");
    # tag = "malpasset_btm"
    filename = string(tag, string(N),"_", string(Integer(round(idx))), ".csv");
    h = readdlm(filename, ',', Float64);
    # h = h - h0
    fig1 = Makie.Figure();
    fig2 = Makie.Figure();
    ax1 = Makie.Axis(fig1[1, 1],
                    aspect = DataAspect(),
                    show_axis=false, resolution = (2500,2500));
    ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
                    elevation = .25*pi, azimuth = -.25*pi,
                    show_axis=false, resolution = (2500, 2500));
    plot_data = Vp*Pq*(btm_q)

    plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
            color=vec(plot_data),
            # shading = false,
            colormap =:gray1,
            # colormap =:viridis,
            );
    w_idx = [];
    for e = 1:size(h,2)
        if maximum(h[:,e]) > 0.0002
            append!(w_idx, e)
        end
    end
    w_idx = append!(w_idx, ocean_idx);
    w_idx = collect(w_idx)

    plot_data = Vp*Pq*(h[:,w_idx])
    # plot_data = Vp*Pq*(h0[:,ocean_idx])
    # plot_data = Vp*Pq*(btm_q[:,sb_idx])

    plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x[:, w_idx], Vp*y[:,w_idx]), set_z_coordinate_to_zero=true),
            color=vec(plot_data),
            # shading = false,
            colormap =:blues,
            # colormap =:viridis,
            );

    # plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x[:, sb_idx], Vp*y[:,sb_idx]), set_z_coordinate_to_zero=true),
    #         color=vec(plot_data),
    #         # shading = false,
    #         colormap =:blues,
    #         # colormap =:viridis,
    #         );

    # xs = xf[sb_b_idx]; ys = yf[sb_b_idx]
    # plt1 = Makie.scatter!(ax1, xs, ys, color = :yellow, markersize = 2)
    # plt = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
    #         color=vec(plot_data), shading = false, colormap = :blues);
    #
    # ax = [ax1, ax2]
    Makie.hidespines!(ax1)
    Makie.hidedecorations!(ax1)
    # save("h1.0_mf_90.png", fig, px_per_unit = 2)
    Makie.tightlimits!(ax1)
    # Makie.Colorbar(fig1[1], plt1)
    # trim!(fig.layout);
    rowsize!(fig1.layout, 1, ax1.scene.px_area[].widths[2]);
    # rowsize!(fig.layout[1,3], 1, ax1.scene.px_area[].widths[2]*1.2)
    # rowsize!(fig.layout, 1, Aspect(3, 1))
    # fig1

    # plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
    #     color=vec(plot_data),
    #     # shading = false,
    #     colormap = :blues
    #     # colormap = :viridis
    #     );
    # Makie.hidespines!(ax2)
    # Makie.hidedecorations!(ax2)
    # fig2

    # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_90.png");
    filename = string( tag,"_", string(Integer(round(idx))), "_btm_90.png");
    save(filename, fig1, px_per_unit = 2)
    # # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_45.png");
    # filename = string(tag,"_", string(Integer(round(idx))),  "_45.png");
    # save(filename, fig2, px_per_unit = 2)
end