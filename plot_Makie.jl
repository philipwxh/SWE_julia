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
