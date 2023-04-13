


subroutine generate_perterb_factors(ens_size, pert_factors)

    Use, Intrinsic :: IEEE_ARITHMETIC
    
    Implicit None

    Integer              :: ens_size   ! = 20
    Real, Intent(Out)    :: pert_factors(8, ens_size)
    Real, parameter      :: temp_corr_hr = 24.0
    Real, parameter      :: spat_corr_km = 50.0

    var_names = ["precipitation_bilinear", "precipitation_conserve", "temperature", &
                 "solar_radiation", "longwave_radiation", &
                 "specific_humidity", "wind_speed", "surface_pressure"]
    var_units = ["mm/s", "mm/s", "K", "W/m2", "W/m2", "kg/kg", "m/s", "Pa"]
    ! only for precip, temp, SWRad, LWRad
    std_dev = [0.5, 2.0, 0.3, 50.0, 0.1, 0.1, 0.2]
    mean_var = [1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0]
    multi_bool = [True, False, True, False, True, True, True]

    ! # only for          precip, SWRad, LWRad
    ! # pert_corr_matrix = [[1.0, -0.8, 0.5],
    ! #                     [-0.8, 1.0, -0.5],
    ! #                     [0.5, -0.5, 1.0]]
    pert_corr_matrix = [[1.0, -0.8, 0.5, 0.0, 0.0, 0.0, 0.0],
                        [-0.8, 1.0, -0.5, 0.0, 0.0, 0.0, 0.0],
                        [0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

    corr_pert = scipy.stats.multivariate_normal.rvs(mean=np.zeros(3), cov=pert_corr_matrix, size=ens_size)

    print(corr_pert)

    pert_factors = np.full((8, ens_size), 1.0)
    pert_factors[0, :] = 1 + std_dev[0] * corr_pert[:, 0]
    for indx in range(7):
        if multi_bool[indx]:
            pert_factors[indx + 1, :] = 1 + std_dev[indx] * corr_pert[:, indx]
        else:
            pert_factors[indx + 1, :] = std_dev[indx] * corr_pert[:, indx]

    ! return pert_factors

end subroutine generate_perterb_factors