
const Depth = 5300
const Nz    = 250

function linear_z_faces(Nz, Depth; first_Δ = 5)

    m = (Depth - (Nz * first_Δ)) / sum(1:Nz-1)
    Δ = m .* (0:Nz-1) .+ first_Δ

    z_faces = zeros(Nz+1)
    for k in Nz:-1:1
        z_faces[k] = z_faces[k+1] - Δ[Nz - k + 1]
    end

    return z_faces
end