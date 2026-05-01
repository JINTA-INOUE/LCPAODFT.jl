function EulerAngle_Spin(Re11, Re22, Re12, Im12)

    r, theta, phi = xyz_to_spherical(2*Re12, -2*Im12, Re11-Re22)

    return 0.5*(Re11+Re22+r), 0.5*(Re11+Re22-r), theta, phi
end