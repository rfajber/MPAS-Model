! Copyright (c) 2023 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_string_utils
!
!> \brief Collection of functions used for string manipulation
!> \author Matthew Dimond
!> \date   25 July 2023
!> \details
!>  This module provides functions and subroutines used for string
!>  manipulations and utilities.
!
!-----------------------------------------------------------------------
module mpas_string_utils

    contains

    !-----------------------------------------------------------------------
    !  routine mpas_split_string
    !
    !> \brief This routine splits a string on a specified delimiting character
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/23/2014
    !> \details This routine splits the given "string" on the delimiter
    !>          character, and returns an array of pointers to the substrings
    !>          between the delimiting characters. Strings are trimmed before
    !>          splitting such that all trailing whitespace is ignored.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_split_string(string, delimiter, subStrings)

        implicit none

        ! Arguments
        character(len=*), intent(in) :: string
        character, intent(in) :: delimiter
        character(len=*), pointer, dimension(:) :: subStrings

        ! Local variables
        character(len=len_trim(string)) :: trimString
        integer :: i, start, index

        trimString = trim(string)
        index = 1

        do i = 1, len(trimString)
            if (trimString(i:i) == delimiter) then
                index = index + 1
            end if
        end do

        allocate(subStrings(1:index))

        start = 1
        index = 1
        do i = 1, len(trimString)
            if (trimString(i:i) == delimiter) then
                subStrings(index) = trimString(start:i-1)
                index = index + 1
                start = i + 1
            end if
        end do
        subStrings(index) = trimString(start:len(trimString))

    end subroutine mpas_split_string

    !-----------------------------------------------------------------------
    !  function mpas_strip_whitespace
    !
    !> \brief Returns string with all whitespace characters removed
    !> \details This function removes every whitespace character from "string",
    !>          including characters that the Fortran intrinsic trim() does not
    !>          treat as blanks: horizontal tab, line feed, vertical tab, form
    !>          feed, carriage return, and NUL. Leading and embedded whitespace
    !>          is removed as well as trailing whitespace. The result is blank
    !>          padded to the length of "string".
    !>
    !>          This is intended for values that are parsed character by
    !>          character rather than with a list-directed read -- date and time
    !>          strings, in particular -- where a stray carriage return from a
    !>          file that has passed through a DOS editor, or NUL padding in a
    !>          character variable written by a netCDF library, would otherwise
    !>          be carried into the parser.
    !
    !-----------------------------------------------------------------------
    function mpas_strip_whitespace(string) result(stringOut)

        implicit none

        ! Arguments
        character(len=*), intent(in) :: string

        ! Result
        character(len=len(string)) :: stringOut

        ! Local variables
        integer :: i, n, ic

        stringOut = ' '
        n = 0

        do i = 1, len(string)
            ic = iachar(string(i:i))

            ! Space, or one of NUL, HT, LF, VT, FF, CR
            if (ic == 32 .or. ic == 0 .or. (ic >= 9 .and. ic <= 13)) cycle

            n = n + 1
            stringOut(n:n) = string(i:i)
        end do

    end function mpas_strip_whitespace

    !-----------------------------------------------------------------------
    !  routine mpas_string_replace
    !
    !> \brief Returns string with charToReplace replaced with targetChar
    !> \author Matthew Dimond
    !> \date   07/26/2023
    !> \details This function replaces all characters matching charToReplace in
    !>          "string" with the char "targetChar" after trimming "string"
    !
    !-----------------------------------------------------------------------
    function mpas_string_replace(string, charToReplace, targetChar) result(stringOut)

        implicit none

        ! Arguments
        character(len=*), intent(in) :: string 
        character, intent(in) :: targetChar, charToReplace

        ! Local variables
        integer :: i

        ! Result
        character(len=len_trim(string)) :: stringOut

        stringOut = trim(string)

        do i = 1, len_trim(string)
            if (string(i:i) == charToReplace) then
                stringOut(i:i) = targetChar
            end if
        end do

    end function mpas_string_replace

end module mpas_string_utils

