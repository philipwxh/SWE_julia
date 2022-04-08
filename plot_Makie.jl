using Triangulate
# rp, sp = NodesAndModes.Tri.equi_nodes(20)

function compute_triangle_area(tri)
   A,B,C = tri
   return .5*(A[1]*(B[2] - C[2]) + B[1]*(C[2]-A[2]) + C[1]*(A[2]-B[2]))
end

function plotting_triangulation(rst_plot,tol=50*eps())

   # on-the-fly triangulation of plotting nodes on the reference element
   triin = Triangulate.TriangulateIO()
   triin.pointlist = permutedims(hcat(rst_plot...))
   triout,_ = triangulate("Q", triin)
   t = triout.trianglelist

   # filter out sliver triangles
   has_volume = fill(true,size(t,2))
   for i in axes(t,2)
       ids = @view t[:,i]
       x_points = @view triout.pointlist[1,ids]
       y_points = @view triout.pointlist[2,ids]
       area = compute_triangle_area(zip(x_points,y_points))
       if abs(area) < tol
           has_volume[i] = false
       end
   end
   return t[:,findall(has_volume)]
end

using CairoMakie
using GeometryBasics
using ColorSchemes

# call `mesh(uplot, rst_plot, xyz_plot)`
function build_plot_data(uplot, rst_plot, xyz_plot;
                        set_z_coordinate_to_zero=false)

   t = permutedims(plotting_triangulation(rst_plot))
   makie_triangles = Makie.to_triangles(t)

   num_elements = size(first(xyz_plot),2)
   trimesh = Vector{GeometryBasics.Mesh{3,Float32}}(undef,num_elements)
   coordinates = zeros(length(first(rst_plot)),3)
   for e = 1:num_elements
       for i = 1:2
           coordinates[:,i] .= view(xyz_plot[i],:,e)
       end
       if !set_z_coordinate_to_zero
           coordinates[:,3] .= view(uplot,:,e)
       end
       trimesh[e] = GeometryBasics.normal_mesh(Makie.to_vertices(coordinates),makie_triangles) # speed this up?
   end
   return merge([trimesh...])
end

function makie_plot_from_csv(filename, save_img=false)
    h  = readdlm(filename, ',', Float64);
    fig1 = Makie.Figure();
    fig2 = Makie.Figure();
    ax1 = Makie.Axis(fig1[1, 1],
                    aspect = DataAspect(),
                    show_axis=false, resolution = (2500,2500));
    ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
                    elevation = .25*pi, azimuth = -.25*pi,
                    show_axis=false, resolution = (2500, 2500));
    plot_data = Vp*Pq*(h)

    plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
            color=vec(plot_data),
            # shading = false,
            colormap =:blues,
            # colormap =:viridis,
            );
    # plt = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
    #         color=vec(plot_data), shading = false, colormap = :blues);
    #
    # ax = [ax1, ax2]
    Makie.hidespines!(ax1)
    Makie.hidedecorations!(ax1)
    # save("h1.0_mf_90.png", fig, px_per_unit = 2)
    Makie.tightlimits!(ax1)
    Makie.Colorbar(fig1[1,2], plt1)
    # trim!(fig.layout);
    rowsize!(fig1.layout, 1, ax1.scene.px_area[].widths[2]);
    # rowsize!(fig.layout[1,3], 1, ax1.scene.px_area[].widths[2]*1.2)
    # rowsize!(fig.layout, 1, Aspect(3, 1))
    fig1

    plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
        color=vec(plot_data),
        # shading = false,
        colormap = :blues
        # colormap = :viridis
        );
    Makie.hidespines!(ax2)
    Makie.hidedecorations!(ax2)
    fig2
    if save_img
        # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_90.png");
        filename = string("db_", tag, "_h_", string(K1D),"_", string(idx-1), "_90.png");
        save(filename, fig1, px_per_unit = 2)
        # # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_45.png");
        filename = string("db_", tag, "_h_", string(K1D),"_", string(idx-1), "_45.png");
        save(filename, fig2, px_per_unit = 2)
    end
end