using FerriteGridUtil
using Ferrite, FerriteGmsh, SparseArrays
using OrderedCollections
import GLMakie: Makie, GeometryBasics, Figure, Axis3, Colorbar, Legend, NoShading, wireframe!, mesh!

using Downloads: download
logo_mesh = "logo.geo"
asset_url = "https://raw.githubusercontent.com/Ferrite-FEM/Ferrite.jl/gh-pages/assets/"
isfile(logo_mesh) || download(string(asset_url, logo_mesh), logo_mesh)

grid = togrid(logo_mesh);

addfacetset!(grid, "top",    x -> x[2] ≈ 1.0) # facets for which x[2] ≈ 1.0 for all nodes
addfacetset!(grid, "left",   x -> abs(x[1]) < 1e-6)
addfacetset!(grid, "bottom", x -> abs(x[2]) < 1e-6);

dim = 2
order = 1 # linear interpolation
ip = Lagrange{RefTriangle, order}()^dim; # vector valued interpolation

qr = QuadratureRule{RefTriangle}(1) # 1 quadrature point
qr_face = FacetQuadratureRule{RefTriangle}(1);

cv = CellValues(qr, ip)
fv = FacetValues(qr_face, ip);

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0, 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "left"),   (x, t) -> 0.0, 1))
close!(ch);

traction(x) = Vec(0.0, 20e3 * x[1]);

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

function compute_von_Mises_stresses(grid, dh, cv, u, C)
	σvm = zeros(getncells(grid))
	for cell in CellIterator(grid)
		reinit!(cv, cell)
		ue = @view u[celldofs(dh, cellid(cell))]
		ε = function_symmetric_gradient(cv, 1, ue)
		σ = C ⊡ ε
		σvm[cellid(cell)] = sqrt(1.5)*norm(dev(σ))
	end
    return σvm
end

σvm = compute_von_Mises_stresses(grid, dh, cv, u, C)
u_nodes = evaluate_at_grid_nodes(dh, u, :u)
mesh = convert_to_makie_mesh(grid; displ=u_nodes)#, disconnectcells=false)

# fig = Figure()
# ax  = Axis3(fig[1,1]; aspect=:data, title="Grid", xlabel="x [m]", ylabel="y [m]", zlabel="z[m]")
# mesh!(ax, mesh; color=:gray, shading=NoShading)
# wireframe!(ax, mesh; color=:black, linewidth=0.2)
# fig

fig = Makie.Figure()
limits = (minimum(norm.(u_nodes)), maximum(norm.(u_nodes)))
ax = Makie.Axis(fig[1,1]; aspect=1, title="displacement")
Makie.mesh!(ax, mesh; color=norm.(u_nodes))
# Makie.wireframe!(ax, mesh; linewidth=0.25)
Makie.Colorbar(fig[1, 2]; vertical=true, colormap=:jet, limits=limits)
fig