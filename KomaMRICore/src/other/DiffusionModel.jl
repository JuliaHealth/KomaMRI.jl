## UNDER CONSTRUCTION!
"""
Matrix A such as A ⋅ b = n × b
"""
cross(n) = begin
   nx,ny,nz = n
    [0 -nz ny;
     nz 0 -nx;
    -ny nx 0]
end

"""
Rodrigues' formula: Rotation matrix that when applied rotates with respect to "n" in an
angle θ anti clock-wise
"""
Un(θ,n) = [1 0 0; 0 1 0; 0 0 1] * cos(θ) + sin(θ) * cross(n) + (1-cos(θ)) * (n * n')
