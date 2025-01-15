using FerriteGridUtil
using Ferrite, FerriteGmsh, SparseArrays
using OrderedCollections
# using GLMakie
using CairoMakie
# import GLMakie: Makie, GeometryBasics, Figure, Axis3, Axis, Colorbar, Legend, NoShading, wireframe!, mesh!

nel = 10
shape = Quadrilateral
refshape = RefQuadrilateral
grid = generate_grid(shape, (nel, nel), Tensors.Vec((0.0, 0.0)), Tensors.Vec((1.0, 1.0)))
addcellset!(grid, "right-half", x -> x[1] >= 0.5 - 1e-3)
# addfacetset!(grid, "top", x -> x[2] ≈ 1.0) # facets for which x[2] ≈ 1.0 for all nodes
# addfacetset!(grid, "left", x -> abs(x[1]) < 1e-6)
# addfacetset!(grid, "bottom", x -> abs(x[2]) < 1e-6);

dim = 2
order = 1 # linear interpolation
ip = Lagrange{refshape, order}()^dim; # vector valued interpolation

qr = QuadratureRule{refshape}(1) # 1 quadrature point
qr_face = FacetQuadratureRule{refshape}(1);

cv = CellValues(qr, ip)
fv = FacetValues(qr_face, ip);

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0, 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "left"),   (x, t) -> 0.0, 1))
close!(ch);

traction(x) = Tensors.Vec(0.0, 20e3 * x[1]);

function assemble_external_forces!(f_ext, dh, facetset, fv, prescribed_traction)
    # Create a temporary array for the facet's local contributions to the external force vector
    fe_ext = zeros(getnbasefunctions(fv))
    for facet in FacetIterator(dh, facetset)
        # Update the facetvalues to the correct facet number
        reinit!(fv, facet)
        # Reset the temporary array for the next facet
        fill!(fe_ext, 0.0)
        # Access the cell's coordinates
        cell_coordinates = getcoordinates(facet)
        for qp in 1:getnquadpoints(fv)
            # Calculate the global coordinate of the quadrature point.
            x = spatial_coordinate(fv, qp, cell_coordinates)
            tₚ = prescribed_traction(x)
            # Get the integration weight for the current quadrature point.
            dΓ = getdetJdV(fv, qp)
            for i in 1:getnbasefunctions(fv)
                Nᵢ = shape_value(fv, qp, i)
                fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
            end
        end
        # Add the local contributions to the correct indices in the global external force vector
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end

Emod = 200e3 # Young's modulus [MPa]
ν = 0.3      # Poisson's ratio [-]

Gmod = Emod / (2(1 + ν))  # Shear modulus
Kmod = Emod / (3(1 - 2ν)) # Bulk modulus

C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2,2}));

function assemble_cell!(ke, cv, C)
    for q_point in 1:getnquadpoints(cv)
        # Get the integration weight for the quadrature point
        dΩ = getdetJdV(cv, q_point)
        for i in 1:getnbasefunctions(cv)
            # Gradient of the test function
            ∇Nᵢ = shape_gradient(cv, q_point, i)
            for j in 1:getnbasefunctions(cv)
                # Symmetric gradient of the trial function
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cv, q_point, j)
                ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
            end
        end
    end
    return ke
end

function assemble_global!(K, dh, cv, C)
    # Allocate the element stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    ke = zeros(n_basefuncs, n_basefuncs)
    # Create an assembler
    assembler = start_assemble(K)
    # Loop over all cells
    for cell in CellIterator(dh)
        # Update the shape function gradients based on the cell coordinates
        reinit!(cv, cell)
        # Reset the element stiffness matrix
        fill!(ke, 0.0)
        # Compute element contribution
        assemble_cell!(ke, cv, C)
        # Assemble ke into K
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end

K = allocate_matrix(dh)
assemble_global!(K, dh, cv, C);

f_ext = zeros(ndofs(dh))
assemble_external_forces!(f_ext, dh, getfacetset(grid, "top"), fv, traction);

apply!(K, f_ext, ch)
u = K \ f_ext;

# function calculate_stresses(grid, dh, cv, u, C)
#     qp_stresses = [
#         [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
#         for _ in 1:getncells(grid)]
#     avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:3)...)
#     for cell in CellIterator(dh)
#         reinit!(cv, cell)
#         cell_stresses = qp_stresses[cellid(cell)]
#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             cell_stresses[q_point] = C ⊡ ε
#         end
#         σ_avg = sum(cell_stresses) / getnquadpoints(cv)
#         avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
#         avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
#         avg_cell_stresses[3][cellid(cell)] = σ_avg[1, 2]
#     end
#     return qp_stresses, avg_cell_stresses
# end

function compute_von_Mises_stresses(dh, cv, u, C)
	σvM = zeros(getncells(dh.grid))
	for cell in CellIterator(dh.grid)
		reinit!(cv, cell)
		ue = @view u[celldofs(dh, cellid(cell))]
		ϵ = function_gradient(cv, 1, ue)
		σ = C ⊡ ϵ
		σvM[cellid(cell)] = sqrt(1.5)*norm(dev(σ))
	end
	return σvM
end


########################################################################
function plot_this(values::Vector, mesh::Makie.GeometryBasics.Mesh)
    fig = Figure()
    colour=:coolwarm
    # colour=:jet
    limits = (minimum(values), maximum(values))
    ax = Axis(fig[1,1]; aspect=1, title=" ")
    mesh!(ax, mesh; color = values, colormap=colour, colorrange=limits)
    wireframe!(ax, mesh; color=:black, linewidth=0.25)
    Colorbar(fig[1, 2]; vertical=true, colormap=colour, limits=limits)
    return fig
end

σvM_cell = compute_von_Mises_stresses(dh, cv, u, C)
# cellset = getcellset(grid, "right-half")
# cellset = OrderedSet(1:getncells(grid))

matrix_field = collect(Float64(i)*j for i in OrderedSet(1:getncells(grid)), j in 1:4)
dmesh, dmakienodes, dmakiecells = convert_to_makie_mesh(grid; cellset=cellset, disconnectcells=true)
# mesh, makienodes, makiecells = convert_to_makie_mesh(grid; disconnectcells=false)

σvM_cell_discon = FerriteGridUtil.disconnect_field(matrix_field, grid; cellset="right-half")

plot_this(σvM_cell_discon, dmesh)