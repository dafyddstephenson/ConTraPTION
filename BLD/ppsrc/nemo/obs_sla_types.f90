MODULE obs_sla_types
   !!======================================================================
   !!                       ***  MODULE obs_sla_type  ***
   !! Observation operators : Satellite identifiers for SLA data
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_sla_types.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obssla_types.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   INTEGER, PARAMETER :: imaxmissions=8
   CHARACTER(len=3) :: cmissions(0:imaxmissions) = &
      & (/ 'XXX', 'E1 ', 'E2 ', 'TP ', 'TPM', 'G2 ', 'J1 ', 'EN ', 'J2 ' /)

END MODULE obs_sla_types
