! SALSA++
! Copyright (c) 2026 Wageningen University & Research
! SPDX-License-Identifier: Apache-2.0
!
! Stand-in for the UCLALES LES module `util` (GPL, not vendored). It provides only the procedures
! that the vendored SALSA files use, with interfaces taken from their call sites.
! getMassIndex and calc_correlation are implemented; the procedures that the phase-1 code path
! does not reach stop with an error if they are called.
module util
    implicit none
    private
    public :: getMassIndex, calc_correlation, getBinMassArray, getBinTotalMass, Extend_last

contains

    ! Index of compound nm of bin nb in a mass array that stores all nbtot bins sequentially per compound.
    pure integer function getMassIndex(nbtot, nb, nm)
        integer, intent(in) :: nbtot, nb, nm
        getMassIndex = (nm-1)*nbtot + nb
    end function getMassIndex

    ! Used by findDry4Wet (constrain_SALSA.f90) to compare the dry composition of an evaporating drop
    ! with that of the candidate a and b aerosol bins. It reproduces upstream util.f90 (commit 8c77958),
    ! which is meant to be the Pearson correlation coefficient
    !     r = (n*sum(x*y) - sum(x)*sum(y)) / (sqrt(n*sum(x**2) - sum(x)**2) * sqrt(n*sum(y**2) - sum(y)**2)).
    ! Two upstream quirks are kept on purpose, so that the reference behaves as upstream:
    ! - the cross term holds only the last product x(n)*y(n), not sum(x*y);
    ! - r = 0 if either variance term is below epsilon(1.) (absolute). For mass mixing ratios (kg/kg)
    !   this is nearly always the case, so findDry4Wet then picks the a bin.
    real function calc_correlation(x, y, n)
        integer, intent(in) :: n
        real, intent(in) :: x(n), y(n)
        real :: sx, sy, sxx, syy, sxy

        sx = sum(x)
        sy = sum(y)
        sxx = sum(x**2)
        syy = sum(y**2)
        sxy = x(n)*y(n)

        if (n*sxx - sx**2 < epsilon(1.) .or. n*syy - sy**2 < epsilon(1.)) then
            calc_correlation = 0.
        else
            calc_correlation = (n*sxy - sx*sy) / (sqrt(n*sxx - sx**2)*sqrt(n*syy - sy**2))
        end if
    end function calc_correlation

    ! Only used by LES output diagnostics and ice (level 5).
    subroutine getBinMassArray(nb, ns, ibin, mass, binmass)
        integer, intent(in) :: nb, ns, ibin
        real, intent(in) :: mass(:)
        real, intent(out) :: binmass(:)
        binmass = 0.
        error stop "util stand-in: getBinMassArray is not implemented (not needed for phase 1)"
    end subroutine getBinMassArray

    ! Only used by LES output diagnostics and ice (level 5).
    subroutine getBinTotalMass(nb, iend, ibin, mass, total)
        integer, intent(in) :: nb, iend, ibin
        real, intent(in) :: mass(:)
        real, intent(out) :: total
        total = 0.
        error stop "util stand-in: getBinTotalMass is not implemented (not needed for phase 1)"
    end subroutine getBinTotalMass

    ! Imported by mo_diag_state.f90, but never called.
    subroutine Extend_last()
        error stop "util stand-in: Extend_last is not implemented (never called)"
    end subroutine Extend_last

end module util
