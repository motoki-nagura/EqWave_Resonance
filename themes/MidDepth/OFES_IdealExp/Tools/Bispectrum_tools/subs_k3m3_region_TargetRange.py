def km_range_180_90_60(period, i_cluster):
    match (period, i_cluster):
        case ("180", 0):  k_range, m_range = (-0.45e-5,-0.05e-5),( 0.10e-2, 0.70e-2)
        case ( "90", 0):  k_range, m_range = (-0.85e-5,-0.40e-5),(-0.55e-2,-0.05e-2)
        case ( "60", 0):  k_range, m_range = (-1.20e-5,-0.60e-5),(-0.20e-2, 0.35e-2)

        case ("180", 1):  k_range, m_range = (-0.45e-5,-0.10e-5),(-0.65e-2,-0.32e-2)
        case ( "90", 1):  k_range, m_range = (-1.28e-5,-0.92e-5),( 0.45e-2, 0.80e-2)
        case ( "60", 1):  k_range, m_range = (-1.45e-5,-1.15e-5),( 0.05e-2, 0.40e-2)

        case ("180", 2):  k_range, m_range = (-0.07e-5, 0.07e-5),(-0.10e-2, 0.22e-2)
        case ( "90", 2):  k_range, m_range = (-0.81e-5,-0.48e-5),(-0.10e-2, 0.25e-2)
        case ( "60", 2):  k_range, m_range = (-0.85e-5,-0.55e-5),(-0.10e-2, 0.40e-2)

        case (_,_):       k_range, m_range = None, None

    return k_range, m_range


def km_range_180_180_90(period, xlab, i_cluster):
    match (period, xlab, i_cluster):
        case ("180", "$k_1$", 0):  k_range, m_range = (-0.45e-5,-0.10e-5),( 0.15e-2, 0.70e-2)
        case ("180", "$k_2$", 0):  k_range, m_range = (-0.65e-5,-0.25e-5),(-0.67e-2,-0.03e-2)
        case ( "90", "$k_3$", 0):  k_range, m_range = (-1.20e-5,-0.55e-5),( 0.00e-2, 0.80e-2)

        case ("180", "$k_1$", 1):  k_range, m_range = (-0.35e-5, 0.05e-5),(-0.40e-2, 0.10e-2)
        case ("180", "$k_2$", 1):  k_range, m_range = (-0.54e-5,-0.11e-5),( 0.18e-2, 0.68e-2)
        case ( "90", "$k_3$", 1):  k_range, m_range = (-0.75e-5,-0.25e-5),( 0.00e-2, 0.60e-2)

        case ("180", "$k_1$", 2):  k_range, m_range = (-0.35e-5,-0.10e-5),( 0.05e-2, 0.25e-2)
        case ("180", "$k_2$", 2):  k_range, m_range = (-1.20e-5,-0.80e-5),( 0.15e-2, 0.65e-2)
        case ( "90", "$k_3$", 2):  k_range, m_range = (-0.75e-5,-0.25e-5),( 0.00e-2, 0.60e-2)

        case (_,_,_):              k_range, m_range = None, None

    return k_range, m_range
