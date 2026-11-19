! SALSA++
! Copyright (c) 2026 Wageningen University & Research
! SPDX-License-Identifier: Apache-2.0
!
! Stand-in for the UCLALES LES module `grid` (GPL, not vendored). It provides only the variables
! that the vendored SALSA files use, written from their use sites. The driver sets the values.
module grid
    implicit none
    save

    ! Grid size including ghost cells. run_SALSA loops over i,j = 3..n-2 and k = nzp-1..2,
    ! so 5 x 5 x 3 gives exactly one interior point.
    integer :: nxp = 5
    integer :: nyp = 5
    integer :: nzp = 3

    integer :: level = 4            ! Thermodynamic level: 4 = SALSA without ice, 5 = with ice.
    real :: dtlt = 1.               ! Time step (s).
    logical :: lemission = .false.  ! Aerosol emissions (not used in the 0D reference).
end module grid
