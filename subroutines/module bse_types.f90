module bse_types

    type bse_coeff
        integer                 :: nc ! number of conduction bands
        integer                 :: nv ! number of valence bands
        integer                 :: nkpts ! number of kpoints
        integer                 :: ic !ic index
        integer                 :: iv !iv index
        integer                 :: ik  !ik index
        integer                 :: nT ! number of transitions
        INTEGER                 :: lmbd ! lambda 
        complex, allocatable    :: A_table(:,:,:,:,:) ! A lambda(nff) c v k q
    end type bse_coeff
end module