import mod_param
import mod_scalars
import mod_boundary
import mod_forces
import mod_grid
import mod_mixing


def mod_arrays():
#  This routine routine allocates and initializes model state arrays


    # LBi =BOUNDS.LBi
    # UBi =BOUNDS.UBi
    # LBj =BOUNDS.LBj
    # UBj =BOUNDS.UBj
    # LBij=BOUNDS.LBij
    # UBij=BOUNDS.UBij


    # Intialize variables
    # ----------------------------------------------------------------------
    initialize_boundary(model)

    if Lclimatology:
        initialize_clima()

    initialize_forces(model)
    initialize_grid  (model)
    initialize_mixing(model)
    initialize_ocean (model)

